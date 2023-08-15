#include <iostream>     /* cout */
#include <vector>
#include <math.h>       /* atan */
#include "ng_func.h"
#include <fstream>
#include <chrono>

#include <Eigen/Dense>
#include "petsc.h"

// Get the name of the variable
#define PRINTER(name) #name

static char help[] = "2D Neuron Growth Solver\n";

int main(int argc, char **argv)
{
	int n_process = atoi(argv[1]);

	int rank, nProcs;
	PetscErrorCode ierr;
	/// start up petsc
	ierr = PetscInitialize(&argc, &argv, (char*)0, help); if (ierr) return ierr;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	MPI_Comm_size(PETSC_COMM_WORLD, &nProcs);

	PetscPrintf(PETSC_COMM_WORLD, "PETSc Check Done!\n");
	ierr = PetscFinalize(); CHKERRQ(ierr);

	// Start Simulation Model
	std::cout << "**********************************************************" << std::endl;
	std::cout << "2D Phase-field Neuron Growth solver using IGA-Collocation" << std::endl;

	// auto start = std::chrono::system_clock::now();
	// auto stop = std::chrono::system_clock::now();
	// auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
	// std::cout << duration.count() << " milliseconds" << std::endl;

	Eigen::setNbThreads(n_process);

	// variable and png save frequency
	int var_save_invl = 100;
	int png_save_invl = 100;
	int png_plot_invl = 100;

	// variable setup
	int numNeuron = 1;
	int gc_sz = 2;
	float kappa = 2.0;
	int M_axon = 100;
	int M_neurites = 60;
	int delta = 0.20;

	// Simulation Parameter Initialization
	// time stepping variables
	float dtime = 1e-2;
	int end_iter = 1; // for debuging
	// int end_iter = 100000; // large enough number

	// setup domain size
	int Nx, Ny;
	neuron_domain_setup(numNeuron,Nx,Ny);

	// B-spline curve order (U,V direction)
	int p = 3;
	int q = 3;
	// vector<float> knotvectorU(Nx+1+6), knotvectorV(Nx+1+6);
	Eigen::VectorXd knotvectorU = Eigen::VectorXd::Zero(Nx+1+6);
	Eigen::VectorXd knotvectorV = Eigen::VectorXd::Zero(Nx+1+6);

	gen_knotvector(knotvectorU,0,Nx,Nx+1,q);
	gen_knotvector(knotvectorV,0,Nx,Nx+1,q);

	// setting lenu lenv this way for easier access to ghost nodes later on
	int lenu = knotvectorU.rows()-2*(p-1);
	int lenv = knotvectorV.rows()-2*(p-1);
	std::cout << "lenu: " << lenu << " | lenv: " << lenv << std::endl;

	// neuron growth variables
	int aniso = 6;
	// float kappa = 2;end // kappa= 2;
	float alpha = 0.9;
	float pix=4.0*atan(1.0);
	float alphOverPix = alpha/pix;
	int gamma = 10;
	float tau = 0.3;
	float M_phi = 60;
	float M_theta = 0.5*M_phi;
	float s_coeff = 0.007;
	// float delta = 0.1;
	float epsilonb = 0.04;

	// Tubulin parameters
	int r = 5;
	float g = 0.1;
	float alpha_t = 0.001;
	float beta_t = 0.001;
	int Diff = 4;
	float source_coeff = 15;

	// tolerance for NR iterations in phi equation
	float residual, tol;
	tol = 1e-4;

	// Initializing phi and concentration based on neuron seed position
	int seed_radius = 20; // Seed size
	Nx = 20; Ny = 20;
	seed_radius = 5;

	Eigen::VectorXd seed_x = Eigen::VectorXd::Zero(numNeuron);
	Eigen::VectorXd seed_y = Eigen::VectorXd::Zero(numNeuron);
	Eigen::MatrixXd phi = Eigen::MatrixXd::Zero(lenu, lenv);
	Eigen::MatrixXd conct = Eigen::MatrixXd::Zero(lenu, lenv);
	Eigen::MatrixXd theta = Eigen::MatrixXd::Random(lenu, lenv);
	Eigen::MatrixXd tempr = Eigen::MatrixXd::Zero(lenu, lenv);

	initialize_neurite_growth(phi, conct, seed_x, seed_y, seed_radius, lenu, lenv, numNeuron);

	// initializing initial phi,theta,tempr for boundary condition (Dirichlet)
	Eigen::MatrixXd phi_initial = phi;
	Eigen::MatrixXd theta_initial = theta;
	Eigen::MatrixXd tempr_initial = tempr;
	// std::cout << "test2.5" << std::endl;
	// for (int i = 2; i < lenu-2; i++)
	// {
	//     for (int j = 2; j < lenv-2; j++)
	//     {
	//         std::cout << i << "|" << j << std::endl;
	//         phi_initial(i,j) = 0;
	//         theta_initial(i,j) = 0;
	//         tempr_initial(i,j) = 0;
	//     }
	// }

	phi_initial = phi_initial.reshaped(lenu*lenv,1);
	theta_initial  = theta_initial.reshaped(lenu*lenv,1);
	tempr_initial  = tempr_initial.reshaped(lenu*lenv,1);

	phi = phi.reshaped(lenu*lenv,1);
	conct = phi.reshaped(lenu*lenv,1);
	tempr = tempr.reshaped(lenu*lenv,1);

	// Expanding domain parameters
	int BC_clearance = 25;
	int expd_sz = 10;

	// Iterating Variable Initialization
	// constructing collocation basis
	int order_deriv = 2;    // highest order of derivatives to calculate
	// std::vector<Eigen::MatrixXd> cm;
	// for (int i = 0; i < 7; ++i) // constructing 7 vars: NuNv, N1uNv, NuN1v, N1uN1v, N2uNv, NuN2v, N2uN2v
	// {
	//     cm.push_back(Eigen::MatrixXd::Zero(lenu*lenv, lenu*lenv));
	// }
	std::vector<Eigen::MatrixXd> cm = vectorMatInit(7, lenu*lenv, lenu*lenv);
	collocationDers(knotvectorU, p, knotvectorV, q, order_deriv, cm);

	Eigen::MatrixXd tmp = cm[0].block(0,0,50,50);
	printArray2TXT(tmp, "./cm0_initial.txt");

	Eigen::MatrixXd lap = cm[5]+cm[6];

	std::cout << "phi shape: " << phi.rows() << "|" << phi.cols() << std::endl;
	phi = linSol(cm[0], phi);

	// std::cout << "Testing second approach" << std::endl;
	// start = std::chrono::system_clock::now();
	// phi = cm[0].lu().solve(phi);
	// stop = std::chrono::system_clock::now();
	// duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
	// std::cout << duration.count() << " milliseconds" << std::endl;

	std::cout << "conct shape: " << conct.rows() << "|" << conct.cols() << std::endl;
	conct = linSol(cm[0], conct);

	// initializing theta and temperature
	std::cout << "theta shape: " << theta.rows() << "|" << theta.cols() << std::endl;   
	theta = linSol(cm[0], theta.reshaped(lenu*lenv,1));

	Eigen::MatrixXd NNtheta = cm[0]*theta; //NNtheta
	Eigen::MatrixXd N1Ntheta = cm[1]*theta;
	Eigen::MatrixXd NN1theta = cm[2]*theta;

	// theta does not evolve over time, only need to compute initially or 
	// expanding domain magnitude of theta gradient
	Eigen::MatrixXd mag_grad_theta = ((cm[1]*theta).cwiseProduct(cm[1]*theta)+\
		(cm[2]*theta).cwiseProduct(cm[2]*theta)).cwiseSqrt();
	Eigen::MatrixXd C0 = matInit(lenu,lenv,0.5);
	C0 = mag_grad_theta.cwiseProduct(matInit(lenu*lenv,1,6*s_coeff));

	// binary ID for boundary location (define 4 edges)
	// id = 1 means there is bc
	Eigen::MatrixXd bcid = Eigen::MatrixXd::Zero(lenu, lenv);
	std::cout << "bcid shape: " << bcid.rows() << "|" << bcid.cols() << std::endl;       
	for (int i = 0; i < lenu; i++)
	{
		bcid(0,i) = 1;
		bcid(lenu-1,0) = 1;
		bcid(i,0) = 1;
		bcid(i,lenv-1) = 1;
	}
	bcid = bcid.reshaped(lenu*lenv,1);

	Eigen::MatrixXd dist = Eigen::MatrixXd::Zero(lenu, lenv);

	std::cout << "Iterating Variable Initialization - Done!" << std::endl;
	std::cout << "**********************************************************" << std::endl;

	std::cout << "Starting Neuron Growth Model transient iterations..." << std::endl;

	Eigen::MatrixXd phiK, R, dR, dp;

	Eigen::ArrayXd E, rMat, gMat, nnT, theta_ori, term_change, aOpMat, gammaMat, C1,\
		NNa, N1Na, NN1a, NNaap, N1Naap, NN1aap, t5, t6, oneMat, twoMat, threeMat,\
		NNpk, N1Npk, NN1pk, N1N1pk, LAPpk, NNtempr, NNct, MphiMat, tauMat;
	aOpMat = matInit(lenu*lenv, 1, alphOverPix).array();
	gammaMat = matInit(lenu*lenv, 1, gamma).array();
	rMat = matInit(lenu*lenv, 1, r).array();
	gMat = matInit(lenu*lenv, 1, g).array();
	theta_ori = matInit(lenu*lenv, 1, 1).array();
	oneMat = matInit(lenu*lenv, 1, 1).array();
	twoMat = matInit(lenu*lenv, 1, 2).array();
	threeMat = matInit(lenu*lenv, 1, 3).array();
	MphiMat = matInit(lenu*lenv, 1, M_phi).array();
	tauMat = matInit(lenu*lenv, 1, tau).array();

	int ind_check;

	for (int iter = 0; iter < end_iter; iter++)
	{
		std::vector<Eigen::MatrixXd> ep_aap = vectorMatInit(5, lenu*lenv, lenu*lenv);
		getEpsilonAndAap(ep_aap, epsilonb, delta, phi, NNtheta, cm, lenu, lenv);

		NNtempr = (cm[0]*tempr).array();
		NNct = (cm[0]*conct).array();

		if (iter<=100)
		{      
			E = aOpMat*(getAtan(gammaMat*(oneMat-NNtempr), lenu, lenv).array());
		} else {
			nnT = theta_ori.reshaped(lenu*lenv,1);
			// adjust tip r g value
			updateRgSg(rMat, gMat, nnT, lenu, lenv); // rMat(nnT==1) = 50; gMat(nnT==1) = 0;
			term_change = regular_Heiviside_fun(rMat*(NNct) - gMat, lenu, lenv);

			E = aOpMat*(getAtan(gammaMat*(term_change*(oneMat-NNtempr)), lenu, lenv).array());
		}
		
		std::cout << "Done E" << std::endl;

		// Phi (Implicit Nonlinear NR method)
		// NR method initial guess (guess current phi)
		phiK = phi;
		// initial residual for NR method
		residual = 2*tol;

		// splitted C0 from C1 because E mag_grad_theta dimension mismatch
		// during domain expansion. Compute here to fix conflict
		C1 = E-C0.array();

		// NR method calculation
		ind_check = 0;
		NNa = (cm[0]*ep_aap[0]).array();
		N1Na = (cm[1]*ep_aap[0]).array();
		NN1a = (cm[2]*ep_aap[0]).array();
		NNaap = (cm[0]*ep_aap[2]).array();
		N1Naap = (cm[1]*ep_aap[2]).array();
		NN1aap = (cm[2]*ep_aap[2]).array();

		std::vector<Eigen::ArrayXd> out = arrayMatInit(4, lenu*lenv, lenu*lenv);
		float dt_t = 0;
		std::cout << "NR iter" << std::endl;
		while (residual >= tol)
		{
			NNpk = (cm[0]*phiK).array();
			N1Npk = (cm[1]*phiK).array();
			NN1pk = (cm[2]*phiK).array();
			N1N1pk = (cm[3]*phiK).array();
			LAPpk = (lap*phiK).array();
			
			std::cout << "test1" << std::endl;

			// term a2
			out[0]  = twoMat*(NNa*N1Na*N1Npk)+(NNa*NNa*LAPpk)+twoMat*(NNa*NN1a*NN1pk);
			std::cout << "test1 - a2" << std::endl;

			// termadx
			out[1] = N1Naap*NN1pk+NNaap*N1N1pk;
			std::cout << "test1 - adx" << std::endl;

			// termady
			out[2] = NN1aap*N1Npk+NNaap*N1N1pk;
			std::cout << "test1 - ady" << std::endl;

			// termNL
			out[3] = -NNpk*NNpk*NNpk+(oneMat-C1)*NNpk*NNpk+C1*NNpk;
			std::cout << "test1 - out3" << std::endl;

			if (dt_t==0) // these terms only needs to be calculated once
			{
				std::cout << "test1 - dtt" << std::endl;

				// terma2_deriv
				t5 =  ((twoMat*NNa*N1Na+N1Naap).matrix()*cm[1]).array()+\
					((NNa*NNa).matrix()*lap).array()+((twoMat*NNa*NN1a-N1Naap).matrix()*cm[2]).array();
				std::cout << "test1 - t5" << std::endl;

			}
			std::cout << "test1 - NL" << std::endl;

		    	// termNL_deriv
			t6 = (- threeMat*NNpk*NNpk+twoMat*(oneMat-C1)*NNpk+C1).matrix()*cm[0];
			std::cout << "test1 - NL_deriv" << std::endl;

			R = ((MphiMat/tauMat*(out[0]+out[1]+out[2]+out[3]))*dtime-NNpk).matrix()+cm[0]*phi;;
			dR = ((MphiMat/tauMat*(t5+t6))*dtime).matrix()-cm[0];

			// check residual and update guess
			R = R - dR*phi_initial;
			stiffMatSetupBCID(dR, R, bcid, phi_initial, lenu, lenv);
			dp = linSol(dR,-R);
			phiK = phiK + dp;
			
		//     max_phi_R = full(max(abs(R)));
		//     if (ind_check >= 100 || max(abs(R))>1e20)
		//         error('Phi NR method NOT converging!-Max residual: %.2d\n',...
		//             max_phi_R);
		//     end

			ind_check += 1;
			dt_t += dtime;
			std::cout << ind_check << std::endl;
		}
	}

	std::cout << "**********************************************************" << std::endl;
	std::cout << "All simulations complete!" << std::endl;
	std::cout << "**********************************************************" << std::endl;

	return 0;
}