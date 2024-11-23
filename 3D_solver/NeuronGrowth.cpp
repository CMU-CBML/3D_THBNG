#include "NeuronGrowth.h"
#include "utils.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <stack>	// for tic toc
#include <ctime>	// for tic toc
#include <time.h>	// for srand
#include <queue>	// for geodesic distance

#include <cfloat>	// for FLT_MAX
#include <limits>	// for numeric_limits

#include <tuple>	// for find vertex with indices and distances

using namespace std;
typedef unsigned int uint;
const float PI = 3.1415926;

std::stack<clock_t> tictoc_stack;
float t_phi, t_syn, t_tub, t_collect, t_write(0), t_total(0);

const int INF = 1e9;

void tic() 
{
	tictoc_stack.push(clock());
}
void toc(float &t) 
{
	t = ((float)(clock() - tictoc_stack.top())) / CLOCKS_PER_SEC;
	tictoc_stack.pop();
}

float MatrixDet(float dxdt[3][3])
{
	float det = dxdt[0][0] * dxdt[1][1] * dxdt[2][2] + dxdt[0][1] * dxdt[1][2] * dxdt[2][0] + dxdt[0][2] * dxdt[2][1] * dxdt[1][0] -
		(dxdt[0][2] * dxdt[1][1] * dxdt[2][0] + dxdt[0][0] * dxdt[1][2] * dxdt[2][1] + dxdt[1][0] * dxdt[0][1] * dxdt[2][2]);

	return det;
}

double atan2_3D(float x, float y, float z) {
    return atan2(y, x) + atan2(z, sqrt(x * x + y * y));
}

void Matrix3DInverse(float dxdt[3][3], float dtdx[3][3])
{
	float det = MatrixDet(dxdt);
	dtdx[0][0] = 1 / det*(dxdt[1][1] * dxdt[2][2] - dxdt[1][2] * dxdt[2][1]);
	dtdx[0][1] = 1 / det*(dxdt[2][1] * dxdt[0][2] - dxdt[0][1] * dxdt[2][2]);
	dtdx[0][2] = 1 / det*(dxdt[0][1] * dxdt[1][2] - dxdt[1][1] * dxdt[0][2]);
	dtdx[1][0] = 1 / det*(dxdt[2][0] * dxdt[1][2] - dxdt[1][0] * dxdt[2][2]);
	dtdx[1][1] = 1 / det*(dxdt[0][0] * dxdt[2][2] - dxdt[0][2] * dxdt[2][0]);
	dtdx[1][2] = 1 / det*(dxdt[1][0] * dxdt[0][2] - dxdt[0][0] * dxdt[1][2]);
	dtdx[2][0] = 1 / det*(dxdt[1][0] * dxdt[2][1] - dxdt[1][1] * dxdt[2][0]);
	dtdx[2][1] = 1 / det*(dxdt[0][1] * dxdt[2][0] - dxdt[0][0] * dxdt[2][1]);
	dtdx[2][2] = 1 / det*(dxdt[0][0] * dxdt[1][1] - dxdt[0][1] * dxdt[1][0]);
}

NeuronGrowth::NeuronGrowth(){
	comm = PETSC_COMM_WORLD; //MPI_COMM_WORLD;
	mpiErr = MPI_Comm_rank(comm, &comRank);
	mpiErr = MPI_Comm_size(comm, &comSize);
	nProcess = comSize;
	judge_phi = 0;
	judge_syn = 0;
	judge_tub = 0;
	sum_grad_phi0_local = 0;
	sum_grad_phi0_global = 0;

	// // integer variable setup
	expandCK_invl		= 1; 		// var_save_invl
	// integer variable setup
	var_save_invl		= 25; 		// var_save_invl
	numNeuron 		= 1;	     	// numNeuron
	// M_axon			= 100;		// M_axon
	// M_neurites 		= 60;  		// M_neurites
	aniso 			= 6;   		// aniso
	gamma 			= 10;  		// gamma
	// seed_radius 		= 2;		// seed radius
	// seed_radius 		= 1;		// seed radius
	// seed_radius 		= 8;		// seed radius
	seed_radius 		= 4;		// seed radius

	// variable setup
	// kappa			= 2.0;		// kappa;
	kappa			= 1.8;		// kappa;
	dt			= 1e-2;		// time step
	Dc			= 3;		// syn D
	// Dc			= 0.1;		// syn D
	alpha			= 0.9;		// alpha
	alphaOverPi		= alpha / PI; 	// alphOverPix
	M_phi			= 10;		// M_phi
	s_coeff			= 0.007;	// s_coeff
	// delta			= 0.20;		// delta
	delta			= 0.50;		// delta
	// epsilonb		= 0.04;		// epsilonb
	epsilonb		= 0.01;		// epsilonb
	r			= 5;		// r
	g			= 0.1;		// g
	alphaT 			= 0.001;	// alpha_t
	betaT			= 0.001;	// beta_t
	Diff			= 0.1;		// Diff
	Diff			= 4;		// Diff
	// Diff			= 60;		// Diff
	source_coeff		= 15;		// source_coeff
	// source_coeff		= 0.012;		// source_coeff
	// source_coeff		= 2;		// source_coeff
}

void NeuronGrowth::AssignProcessor(vector<vector<int>> &ele_proc)
{
	ele_process.clear();
	for (int i = 0; i < ele_proc[comRank].size(); i++)
		ele_process.push_back(ele_proc[comRank][i]);
}

void NeuronGrowth::SetVariables(string fn_par)
{
	string fname(fn_par), stmp;
	stringstream ss;
	ifstream fin;
	fin.open(fname);
	if (fin.is_open()) {
		fin >> stmp >> var_save_invl;
		fin >> stmp >> expandCK_invl;
		fin >> stmp >> gc_sz;
		fin >> stmp >> aniso;
		fin >> stmp >> gamma;
		fin >> stmp >> seed_radius;
		fin >> stmp >> kappa;
		fin >> stmp >> dt;
		fin >> stmp >> Dc;
		fin >> stmp >> alpha;
		alphaOverPi = alpha / PI; 	// alphOverPi
		fin >> stmp >> M_phi;
		fin >> stmp >> s_coeff;
		fin >> stmp >> delta;
		fin >> stmp >> epsilonb;
		fin >> stmp >> r;
		fin >> stmp >> g;
		fin >> stmp >> alphaT;
		fin >> stmp >> betaT;
		fin >> stmp >> Diff;
		fin >> stmp >> source_coeff;

		fin.close();
		PetscPrintf(PETSC_COMM_WORLD, "Parameter Loaded!\n");
	} else {
		PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s!\n", fname.c_str());
	}
}

void NeuronGrowth::InitializeProblemNG(const int n_bz, vector<Vertex3D>& cpts, vector<Vertex3D> prev_cpts, vector<vector<float>> &NGvars, vector<array<float, 3>> &seed)
{
	MPI_Barrier(PETSC_COMM_WORLD);
	PetscPrintf(PETSC_COMM_WORLD, "Initializing-----------------------------------------------------------------\n");

	/*Initialize parameters*/
	int cpt_sz = cpts.size();
	// GaussInfo(3);
	n_bzmesh = n_bz;
	N_0.resize(cpt_sz);

	// std::cout << cpts.size() << std::endl << std::endl << std::endl;
	float max_px, min_px, max_py, min_py, max_pz, min_pz;
	for (int i = 0; i<cpts.size(); i++) {
		if (i==0) {
			max_x = cpts[i].coor[0];
			min_x = cpts[i].coor[0];
			max_y = cpts[i].coor[1];
			min_y = cpts[i].coor[1];
			max_z = cpts[i].coor[2];
			min_z = cpts[i].coor[2];
		} else {
			max_x = max(cpts[i].coor[0],max_x);
			min_x = min(cpts[i].coor[0],min_x);
			max_y = max(cpts[i].coor[1],max_y);
			min_y = min(cpts[i].coor[1],min_y);
			max_z = max(cpts[i].coor[2],max_z);
			min_z = min(cpts[i].coor[2],min_z);
		}
		// std::cout << cpts[i].coor[0] << " " << cpts[i].coor[1] << " " << cpts[i].coor[2] << std::endl;
	}

	for (int i = 0; i<prev_cpts.size(); i++) {
		if (i==0) {
			max_px = prev_cpts[i].coor[0];
			min_px = prev_cpts[i].coor[0];
			max_py = prev_cpts[i].coor[1];
			min_py = prev_cpts[i].coor[1];
			max_pz = prev_cpts[i].coor[2];
			min_pz = prev_cpts[i].coor[2];
		} else {
			max_px = max(prev_cpts[i].coor[0],max_px);
			min_px = min(prev_cpts[i].coor[0],min_px);
			max_py = max(prev_cpts[i].coor[1],max_py);
			min_py = min(prev_cpts[i].coor[1],min_py);
			max_pz = max(prev_cpts[i].coor[2],max_pz);
			min_pz = min(prev_cpts[i].coor[2],min_pz);
		}
		// std::cout << cpts[i].coor[0] << " " << cpts[i].coor[1] << " " << cpts[i].coor[2] << std::endl;
	}

	PetscPrintf(PETSC_COMM_WORLD, "Checking maximum x, y, and z value-----------------------------------------------\n");
	PetscPrintf(PETSC_COMM_WORLD, "max x: %f min x: %f | max y: %f min y %f | max z: %f min z %f\n", max_x, min_x, max_y, min_y, max_z, min_z);

	float r, x, y, z; // radius, x, y
	int nx = sqrt(cpts.size()); // only used at the beginning (square mesh)
	float dx = max_x/(float)nx; // only used at the beginning (square mesh)
	srand(time(NULL)); // random seed

	float maxDistI(0);
	// std::cout << dx << std::endl;
	if (n == 0) { // Initialize
		for (int i = 0; i<cpts.size(); i++) {
			phi.push_back(0.f);
			tub.push_back(0.f);
			theta.push_back((float)(rand()%100)/(float)100); // orientation theta
			// polar.push_back((float)(rand()%100)/(float)100); // orientation theta
			// azimuth.push_back((float)(rand()%100)/(float)100); // orientation theta
			syn.push_back(0.f); // synaptogenesis
			Mphi.push_back(M_phi);
			distI.push_back(0);
		}
		for (int i = 0; i<cpts.size(); i++) {
			x = cpts[i].coor[0];
			y = cpts[i].coor[1];
			z = cpts[i].coor[2];

			// calculate new boundary label
			if (x==min_x || x==max_x || y==min_y || y==max_y || z==min_z || z==max_z) {
				cpts[i].label = 1;
			} else {
				cpts[i].label = 0;
			}   

			for (int j = 0; j<seed.size(); j++) { // multiple neurons
				r = sqrt( (x-(float)seed[j][0])*(x-(float)seed[j][0]) 
					+(y-(float)seed[j][1])*(y-(float)seed[j][1])
					+(z-(float)seed[j][2])*(z-(float)seed[j][2]) ); // radius of initial soma
				if (r <= (seed_radius)) { // initial soma
					phi[i] = 1.0f;
					tub[i] = (0.5+0.5*tanh((sqrt(seed_radius)-r)/2)); // based eqn in literature
				} 
			} 

			auto closestVertices = FindClosestVerticesWithIndicesAndDistances(cpts, cpts[i], 6);
			for (const auto& item : closestVertices) {
				// float distance = std::get<2>(item);
				distI[i] += std::get<2>(item);;
			}

			maxDistI = max(distI[i], maxDistI);

		}
		phi_0 = phi; // initial phi
		tub_0 = tub; // initial tub
		
	} else { // Collect from NGvars - continue simulation
		phi.clear(); 	phi.resize(cpts.size());
		syn.clear(); 	syn.resize(cpts.size());
		tub.clear(); 	tub.resize(cpts.size());
		theta.clear();	theta.resize(cpts.size());
		phi_0.clear(); 	phi_0.resize(cpts.size());
		tub_0.clear(); 	tub_0.resize(cpts.size());
		distI.clear(); 	distI.resize(cpts.size());

		for (int i = 0; i < cpts.size(); i++) {
			if ((cpts[i].coor[0] > min_px) && (cpts[i].coor[0] < max_px)
				&& (cpts[i].coor[1] > min_py) && (cpts[i].coor[1] < max_py)
				&& (cpts[i].coor[2] > min_pz) && (cpts[i].coor[2] < max_pz)) {
				
				auto closestVertices = FindClosestVerticesWithIndicesAndDistances(prev_cpts, cpts[i], 6);

				int index;
				if (SearchVertex(prev_cpts, cpts[i].coor[0], cpts[i].coor[1], cpts[i].coor[2], index)) {
					// Exact match found, no need for interpolation
					phi[i] = NGvars[0][index];
					syn[i] = NGvars[1][index];
					tub[i] = NGvars[2][index];
					theta[i] = NGvars[3][index];
					phi_0[i] = NGvars[4][index];
					tub_0[i] = NGvars[5][index];

					for (const auto& item : closestVertices) {
						int index = std::get<1>(item);
						distI[i] += std::get<2>(item);
					}

				} 
				else {
					// auto closestVertices = FindClosestVerticesWithIndices(prev_cpts, cpts[i], 8);

					// for (const auto& vertex : closestVertices) {
					// 	phi[i] += NGvars[0][vertex.second];
					// 	syn[i] += NGvars[1][vertex.second];
					// 	tub[i] += NGvars[2][vertex.second];
					// 	theta[i] += NGvars[3][vertex.second];
					// 	phi_0[i] += NGvars[4][vertex.second];
					// 	tub_0[i] += NGvars[5][vertex.second];
					// }

					// auto closestVertices = FindClosestVerticesWithIndicesAndDistances(prev_cpts, cpts[i], 8);

					for (const auto& item : closestVertices) {
						// const Vertex3D& vertex = std::get<0>(item);
						int index = std::get<1>(item);
						// float distance = std::get<2>(item);
						phi[i] += NGvars[0][index];
						syn[i] += NGvars[1][index];
						tub[i] += NGvars[2][index];
						theta[i] += NGvars[3][index];
						phi_0[i] += NGvars[4][index];
						tub_0[i] += NGvars[5][index];
						distI[i] += std::get<2>(item);
					}

					phi[i] = phi[i] / closestVertices.size();
					syn[i] = syn[i] / closestVertices.size();
					tub[i] = tub[i] / closestVertices.size();
					theta[i] = theta[i] / closestVertices.size();
					phi_0[i] = phi_0[i] / closestVertices.size();
					tub_0[i] = tub_0[i] / closestVertices.size();
				}
			} else {
				phi[i] = 0;
				syn[i] = 0;
				tub[i] = 0;
				theta[i] = (float)(rand()%100)/(float)100;
				phi_0[i] = 0;
				tub_0[i] = 0;

				auto closestVertices = FindClosestVerticesWithIndicesAndDistances(cpts, cpts[i], 6);
				for (const auto& item : closestVertices) {
					// int index = std::get<1>(item);
					// float distance = std::get<2>(item);
					distI[i] += std::get<2>(item);
				}
			}

			maxDistI = max(distI[i], maxDistI);
			
		}
		// Mphi.push_back(M_phi);
	}

	for (int i = 0; i < distI.size(); i++) {
		distI[i] = distI[i]/maxDistI;
	}

	this->cpts = cpts; // for passing into formFunction

	// // writing initialized variables for debugging purposes
	// bool visualization = true;
	// PrintVec2TXT(phi, "/home/kuanrenqian/Research/3D_NeuronGrowth/io/2DNG/var_check/phi_inital.txt", visualization);
	// PrintVec2TXT(tub, "/home/kuanrenqian/Research/3D_NeuronGrowth/io/2DNG/var_check/tub_inital.txt", visualization);
	// PrintVec2TXT(theta, "/home/kuanrenqian/Research/3D_NeuronGrowth/io/2DNG/var_check/theta_inital.txt", visualization);
	// PrintVec2TXT(syn, "/home/kuanrenqian/Research/3D_NeuronGrowth/io/2DNG/var_check/syn_inital.txt", visualization);

	/*Initialize petsc vector, matrix*/
	PetscInt mat_dim = cpt_sz;

	// 250 = 2*5^3. (2 var, 5 basis, 3D) || 25 = 1*5^2. (1 var, 5 basis, 2D)	
	ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, mat_dim, &temp_phi);
	
	ierr = MatCreate(PETSC_COMM_WORLD, &J);
	ierr = MatSetSizes(J, PETSC_DECIDE, PETSC_DECIDE, mat_dim, mat_dim);
	ierr = MatSetType(J, MATMPIAIJ);
	ierr = MatMPIAIJSetPreallocation(J, 250, NULL, 250, NULL);
	ierr = MatSetOption(J, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
	ierr = MatSetUp(J);

	ierr = MatCreate(PETSC_COMM_WORLD, &GK_syn);
	ierr = MatSetSizes(GK_syn, PETSC_DECIDE, PETSC_DECIDE, mat_dim, mat_dim);
	ierr = MatSetType(GK_syn, MATMPIAIJ);
	ierr = MatMPIAIJSetPreallocation(GK_syn, 250, NULL, 250, NULL);
	ierr = MatSetOption(GK_syn, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
	ierr = MatSetUp(GK_syn);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, mat_dim, &GR_syn);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, mat_dim, &temp_syn);

	ierr = MatCreate(PETSC_COMM_WORLD, &GK_tub);
	ierr = MatSetSizes(GK_tub, PETSC_DECIDE, PETSC_DECIDE, mat_dim, mat_dim);
	ierr = MatSetType(GK_tub, MATMPIAIJ);
	ierr = MatMPIAIJSetPreallocation(GK_tub, 250, NULL, 250, NULL);
	ierr = MatSetOption(GK_tub, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
	ierr = MatSetUp(GK_tub);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, mat_dim, &GR_tub);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, mat_dim, &temp_tub);

	vars.reserve(21*sizeof(float));
	vars.resize(21);

	pre_EMatrixSolve.clear();
	pre_EVectorSolve.clear();
}


void NeuronGrowth::CheckVar(string fn, vector<Vertex3D> cpts, vector<float> input)
{
	string fname(fn + "_" + to_string(n) + ".vtk");
	ofstream fout;
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << cpts.size() << " float\n";
		for (uint i = 0; i<cpts.size(); i++)
		{
			fout << cpts[i].coor[0] << " " << cpts[i].coor[1] << " " << cpts[i].coor[2] << "\n";
		}

		fout << "POINT_DATA " << input.size() << "\nSCALARS pact float 1\nLOOKUP_TABLE default\n";
		for (uint i = 0; i<input.size(); i++)
		{
			fout << input[i] << "\n";
		}


		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}
}

void NeuronGrowth::ToPETScVec(vector<float> input, Vec& petscVec)
{
	PetscInt localStart, localEnd;
	VecGetOwnershipRange(petscVec, &localStart, &localEnd);
	for (int i = localStart; i<localEnd; i++) {
		PetscScalar value = input[i]; // Adjust the index to match the localValues
		VecSetValues(petscVec, 1, &i, &value, INSERT_VALUES);
	}
	VecAssemblyBegin(petscVec);
	VecAssemblyEnd(petscVec);
}

void NeuronGrowth::ReadBezierElementProcess(string fn)
{
	string stmp;
	int npts, neles, nfunctions, itmp, itmp1;
	int add(0);

	string fname_cmat = fn + "cmat.txt";
	ifstream fin_cmat;
	fin_cmat.open(fname_cmat);

	int numPTs(0); 

	if (fin_cmat.is_open()) {
		fin_cmat >> neles;
		bzmesh_process.resize(ele_process.size());
		for (int i = 0; i<neles; i++) {
			// if (i == ele_process[add]) {
			if ((i == ele_process[add]) && (add < ele_process.size())) {
				fin_cmat >> itmp >> nfunctions >> bzmesh_process[add].type;
				bzmesh_process[add].cmat.resize(nfunctions);
				bzmesh_process[add].IEN.resize(nfunctions);
				for (int j = 0; j < nfunctions; j++) {
					fin_cmat >> bzmesh_process[add].IEN[j];
				}

				for (int j = 0; j < nfunctions; j++) {
					for (int k = 0; k < 64; k++) {
						fin_cmat >> bzmesh_process[add].cmat[j][k];
					}
				}
				add++;
			}
			else {
				fin_cmat >> stmp >> nfunctions >> itmp;
				for (int j = 0; j < nfunctions; j++)
					fin_cmat >> stmp;
				for (int j = 0; j < nfunctions; j++)
					for (int k = 0; k < 64; k++)
						fin_cmat >> stmp;
			}
		}
		fin_cmat.close();
		PetscPrintf(PETSC_COMM_WORLD, "Bezier Matrices Loaded!\n");
	}
	else {
		PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s!\n", fname_cmat.c_str());
	}

	string fname_bzpt = fn + "bzpt.txt";
	ifstream fin_bzpt;
	fin_bzpt.open(fname_bzpt);
	add = 0;
	if (fin_bzpt.is_open()) {
		fin_bzpt >> npts;
		getline(fin_bzpt, stmp);
		for (int e = 0; e < neles; e++) {
			// if (e == ele_process[add]) {
			if ((e == ele_process[add]) && (add < ele_process.size())) {
				bzmesh_process[add].pts.resize(bzpt_num);
				for (int i = 0; i < bzpt_num; i++) {
					fin_bzpt >> bzmesh_process[add].pts[i][0] >> bzmesh_process[add].pts[i][1] >> bzmesh_process[add].pts[i][2];
				}
				add++;
			}
			else {
				for (int i = 0; i < bzpt_num; i++)
					fin_bzpt >> stmp >> stmp >> stmp;
			}
		}
		fin_bzpt.close();
		PetscPrintf(PETSC_COMM_WORLD, "Bezier Points Loaded!\n");
	}
	else {
		PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s!\n", fname_bzpt.c_str());
	}
}

// void NeuronGrowth::ReadBezierElementProcess(string fn)
// {
// 	string stmp;
// 	int npts, neles, nfunctions, itmp, itmp1;
// 	int add(0);

// 	string fname_cmat = fn + "cmat.txt";
// 	ifstream fin_cmat;
// 	fin_cmat.open(fname_cmat);

// 	int numPTs(0); 

// 	if (fin_cmat.is_open()) {
// 		fin_cmat >> neles;
// 		bzmesh_process.resize(ele_process.size());
// 		for (int i = 0; i<neles; i++) {
// 			if (i == ele_process[add]) {
// 				fin_cmat >> itmp >> nfunctions >> bzmesh_process[add].type;
// 				bzmesh_process[add].cmat.resize(nfunctions);
// 				bzmesh_process[add].IEN.resize(nfunctions);
// 				for (int j = 0; j < nfunctions; j++) {
// 					fin_cmat >> bzmesh_process[add].IEN[j];
// 				}

// 				for (int j = 0; j < nfunctions; j++) {
// 					for (int k = 0; k < 64; k++) {
// 						fin_cmat >> bzmesh_process[add].cmat[j][k];
// 					}
// 				}
// 				add++;
// 			}
// 			else {
// 				fin_cmat >> stmp >> nfunctions >> itmp;
// 				for (int j = 0; j < nfunctions; j++)
// 					fin_cmat >> stmp;
// 				for (int j = 0; j < nfunctions; j++)
// 					for (int k = 0; k < 64; k++)
// 						fin_cmat >> stmp;
// 			}
// 		}
// 		fin_cmat.close();
// 		PetscPrintf(PETSC_COMM_WORLD, "Bezier Matrices Loaded!\n");
// 	}
// 	else {
// 		PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s!\n", fname_cmat.c_str());
// 	}

// 	string fname_bzpt = fn + "bzpt.txt";
// 	ifstream fin_bzpt;
// 	fin_bzpt.open(fname_bzpt);
// 	add = 0;
// 	if (fin_bzpt.is_open()) {
// 		fin_bzpt >> npts;
// 		getline(fin_bzpt, stmp);
// 		for (int e = 0; e < neles; e++) {
// 			if (e == ele_process[add]) {
// 				bzmesh_process[add].pts.resize(bzpt_num);
// 				for (int i = 0; i < bzpt_num; i++) {
// 					fin_bzpt >> bzmesh_process[add].pts[i][0] >> bzmesh_process[add].pts[i][1] >> bzmesh_process[add].pts[i][2];
// 				}
// 				add++;
// 			}
// 			else {
// 				for (int i = 0; i < bzpt_num; i++)
// 					fin_bzpt >> stmp >> stmp >> stmp;
// 			}
// 		}
// 		fin_bzpt.close();
// 		PetscPrintf(PETSC_COMM_WORLD, "Bezier Points Loaded!\n");
// 	}
// 	else {
// 		PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s!\n", fname_bzpt.c_str());
// 	}
// }

void NeuronGrowth::GaussInfo(int ng)
{
	Gpt.clear();
	wght.clear();
	switch (ng)
	{
		case 2:
		{
			Gpt.resize(ng);
			wght.resize(ng);
			Gpt[0] = 0.2113248654051871;
			Gpt[1] = 0.7886751345948129;
			wght[0] = 1.;
			wght[1] = 1.;
			break;
		}
		case 3:
		{
			Gpt.resize(ng);
			wght.resize(ng);
			Gpt[0] = 0.1127016653792583;
			Gpt[1] = 0.5;
			Gpt[2] = 0.8872983346207417;
			wght[0] = 0.5555555555555556;
			wght[1] = 0.8888888888888889;
			wght[2] = 0.5555555555555556;
			break;
		}
		case 4:
		{
			Gpt.resize(ng);
			wght.resize(ng);
			Gpt[0] = 0.06943184420297371;
			Gpt[1] = 0.33000947820757187;
			Gpt[2] = 0.6699905217924281;
			Gpt[3] = 0.9305681557970262;
			wght[0] = 0.3478548451374539;
			wght[1] = 0.6521451548625461;
			wght[2] = 0.6521451548625461;
			wght[3] = 0.3478548451374539;
			break;
		}
		case 5:
		{
			Gpt.resize(ng);
			wght.resize(ng);
			Gpt[0] = 0.046910077030668;
			Gpt[1] = 0.2307653449471585;
			Gpt[2] = 0.5;
			Gpt[3] = 0.7692346550528415;
			Gpt[4] = 0.953089922969332;
			wght[0] = 0.2369268850561891;
			wght[1] = 0.4786286704993665;
			wght[2] = 0.5688888888888889;
			wght[3] = 0.4786286704993665;
			wght[4] = 0.2369268850561891;
			break;
		}
		default:
		{
			Gpt.resize(2);
			wght.resize(2);
			Gpt[0] = 0.2113248654051871;
			Gpt[1] = 0.7886751345948129;
			wght[0] = 1.;
			wght[1] = 1.;
			break;
		}
	}
}

void NeuronGrowth::BasisFunction(float u, float v, float w, const vector<array<float, 3>>& pt, const vector<array<float, 64>> &cmat, vector<float> &Nx, vector<array<float, 3>> &dNdx, float dudx[3][3], float& detJ)
{
	float Nu[4] = { (1. - u)*(1. - u)*(1. - u), 3.*(1. - u)*(1. - u)*u, 3.*(1. - u)*u*u, u*u*u };
	float Nv[4] = { (1. - v)*(1. - v)*(1. - v), 3.*(1. - v)*(1. - v)*v, 3.*(1. - v)*v*v, v*v*v };
	float Nw[4] = { (1. - w)*(1. - w)*(1. - w), 3.*(1. - w)*(1. - w)*w, 3.*(1. - w)*w*w, w*w*w };
	float dNdu[4] = { -3.*(1. - u)*(1. - u), 3. - 12.*u + 9.*u*u, 3.*(2. - 3.*u)*u, 3.*u*u };
	float dNdv[4] = { -3.*(1. - v)*(1. - v), 3. - 12.*v + 9.*v*v, 3.*(2. - 3.*v)*v, 3.*v*v };
	float dNdw[4] = { -3.*(1. - w)*(1. - w), 3. - 12.*w + 9.*w*w, 3.*(2. - 3.*w)*w, 3.*w*w };
	float dNdt[bzpt_num][3];
	float Nx_bz[bzpt_num];
	float dNdx_bz[bzpt_num][3];

	Nx.clear();
	dNdx.clear();
	Nx.resize(cmat.size(), 0);
	dNdx.resize(cmat.size(), { 0, 0, 0 });

	int i, j, k, a, b, c, loc;
	loc = 0;
	for (i = 0; i<4; i++){
		for (j = 0; j<4; j++){
			for (k = 0; k < 4; k++)	{
				Nx_bz[loc] = Nu[k] * Nv[j] * Nw[i];
				dNdt[loc][0] = dNdu[k] * Nv[j] * Nw[i];
				dNdt[loc][1] = Nu[k] * dNdv[j] * Nw[i];
				dNdt[loc][2] = Nu[k] * Nv[j] * dNdw[i];
				loc++;
			}
		}
	}

	float dxdt[3][3] = { { 0 } };
	for (loc = 0; loc < bzpt_num; loc++)
		for (a = 0; a<3; a++)
			for (b = 0; b<3; b++)
				dxdt[a][b] += pt[loc][a] * dNdt[loc][b];

	float dtdx[3][3] = { { 0 } };
	Matrix3DInverse(dxdt, dtdx);

	//1st derivatives
	for (i = 0; i<bzpt_num; i++){
		dNdx_bz[i][0] = dNdt[i][0] * dtdx[0][0] + dNdt[i][1] * dtdx[1][0] + dNdt[i][2] * dtdx[2][0];
		dNdx_bz[i][1] = dNdt[i][0] * dtdx[0][1] + dNdt[i][1] * dtdx[1][1] + dNdt[i][2] * dtdx[2][1];
		dNdx_bz[i][2] = dNdt[i][0] * dtdx[0][2] + dNdt[i][1] * dtdx[1][2] + dNdt[i][2] * dtdx[2][2];
	}
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			dudx[i][j] = dtdx[i][j];

	detJ = MatrixDet(dxdt);
	detJ = 0.125*detJ;

	for (i = 0; i < cmat.size(); i++){
		for (j = 0; j < bzpt_num; j++){
			Nx[i] += cmat[i][j] * Nx_bz[j];
			for (int m = 0; m < 3; m++)	{
				dNdx[i][m] += cmat[i][j] * dNdx_bz[j][m];
			}
		}
	}
}

void NeuronGrowth::ApplyBoundaryCondition(const float bc_value, int pt_num, int variable_num, vector<vector<float>> &EMatrixSolve, vector<float> &EVectorSolve)
{
	int j, k;
	int nen = EVectorSolve.size();
	for (j = 0; j < nen; j++)
	{
		EVectorSolve[j] -= bc_value * EMatrixSolve[j][pt_num + variable_num * nen];
	}
	for (j = 0; j < nen; j++)
	{
		EMatrixSolve[j][pt_num + variable_num * nen] = 0.0;
		EMatrixSolve[pt_num + variable_num * nen][j] = 0.0;
	}
	EMatrixSolve[pt_num + variable_num * nen][pt_num + variable_num * nen] = 1.0;
	EVectorSolve[pt_num + variable_num * nen] = bc_value;
}

void NeuronGrowth::MatrixAssembly(vector<vector<float>> &EMatrixSolve, const vector<int> &IEN, Mat &GK)
{
	int i, j, A, B, m, n;
	int row_start, row_end, row_now;
	int add = 0;
	int nen = IEN.size();

	PetscInt *nodeList = new PetscInt[nen * 1];
	PetscReal *tmpGK = new PetscReal[nen * 1 * nen * 1];

	for (m = 0; m < IEN.size(); m++)
	{
		A = IEN[m];
		nodeList[1 * m] = 1 * A;
		for (n = 0; n < IEN.size(); n++)
		{
				tmpGK[add] = EMatrixSolve[m][n];
				add++;
		}
	}
	MatSetValues(GK, nen * 1, nodeList, nen * 1, nodeList, tmpGK, ADD_VALUES);
	delete nodeList;
	delete tmpGK;
}

void NeuronGrowth::ResidualAssembly(vector<float> &EVectorSolve, const vector<int> &IEN, Vec &GR)
{
	int i, m, A;
	int add = 0;
	int nen = IEN.size();

	PetscInt *nodeList = new PetscInt[nen];
	PetscReal *tmpGR = new PetscReal[nen];

	for (i = 0; i < IEN.size(); i++)
	{
		A = IEN[i];
		nodeList[i * 1 + 0] = A * 1 + 0;
		tmpGR[add] = EVectorSolve[i];
		add += 1;
	}
	VecSetValues(GR, nen * 1, nodeList, tmpGR, ADD_VALUES);
	delete nodeList;
	delete tmpGR;
}

void NeuronGrowth::MatrixAssembly_insert(vector<vector<float>> &EMatrixSolve, const vector<int> &IEN, Mat &GK)
{
	int i, j, A, B, m, n;
	int row_start, row_end, row_now;
	int add = 0;
	int nen = IEN.size();

	PetscInt *nodeList = new PetscInt[nen * 1];
	PetscReal *tmpGK = new PetscReal[nen * 1 * nen * 1];

	for (m = 0; m < IEN.size(); m++)
	{
		A = IEN[m];
		nodeList[1 * m] = 1 * A;
		for (n = 0; n < IEN.size(); n++)
		{
				tmpGK[add] = EMatrixSolve[m][n];
				add++;
		}
	}
	MatSetValues(GK, nen * 1, nodeList, nen * 1, nodeList, tmpGK, INSERT_VALUES);
	delete nodeList;
	delete tmpGK;
}

void NeuronGrowth::ResidualAssembly_insert(vector<float> &EVectorSolve, const vector<int> &IEN, Vec &GR)
{
	int i, m, A;
	int add = 0;
	int nen = IEN.size();

	PetscInt *nodeList = new PetscInt[nen];
	PetscReal *tmpGR = new PetscReal[nen];

	for (i = 0; i < IEN.size(); i++)
	{
		A = IEN[i];
		nodeList[i * 1 + 0] = A * 1 + 0;
		tmpGR[add] = EVectorSolve[i];
		add += 1;
	}
	VecSetValues(GR, nen * 1, nodeList, tmpGR, INSERT_VALUES);
	delete nodeList;
	delete tmpGR;
}

void NeuronGrowth::MatrixAssembly_2var(vector<vector<float>>& EMatrixSolve, const vector<int>& IEN, Mat& GK)
{
	int i, j, A, B, m, n;
	int row_start, row_end, row_now;
	int add = 0;
	int nen = IEN.size();


	PetscInt *nodeList = new PetscInt[nen * 2];	
	PetscReal *tmpGK = new PetscReal[nen * 2 * nen * 2];
	
	for (m = 0; m<IEN.size(); m++){
		A = IEN[m];
		for (i = 0; i < 2; i++)	{
			nodeList[2 * m + i] = 2 * A + i;
			for (n = 0; n<IEN.size(); n++){
				for (j = 0; j < 2; j++)	{
					tmpGK[add] = EMatrixSolve[m + i*nen][n + j*nen];
					add++;
				}   
			}
		}
	}	
	MatSetValues(GK, nen * 2, nodeList, nen * 2, nodeList, tmpGK, ADD_VALUES);
	delete nodeList;
	delete tmpGK;
}

void NeuronGrowth::ResidualAssembly_2var(vector<float>& EVectorSolve, const vector<int>& IEN, Vec& GR)
{
	int i, m, A;
	int add = 0;
	int nen = IEN.size();

	PetscInt *nodeList = new PetscInt[nen * 2];
	PetscReal *tmpGR = new PetscReal[nen * 2];
	
	for (i = 0; i<IEN.size(); i++){
		A = IEN[i];
		nodeList[i * 2 + 0] = A * 2 + 0;
		nodeList[i * 2 + 1] = A * 2 + 1;
		tmpGR[add] = EVectorSolve[i];
		tmpGR[add + 1] = EVectorSolve[i + nen];
		add+=2;
	}
	VecSetValues(GR, nen * 2, nodeList, tmpGR, ADD_VALUES);
	delete nodeList;
	delete tmpGR;
}

void NeuronGrowth::VisualizeVTK_ControlMesh(const vector<Vertex3D> &spt, const vector<Element3D> &mesh, int step, string fn, vector<float> var, string varName)
{
	string fname;
	stringstream ss;
	ss << std::setw(6) << std::setfill('0') << step;
	ofstream fout;
	unsigned int i;
	fname = fn + "/controlmesh_" + varName + ss.str() + ".vtk";
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << spt.size() << " float\n";
		for (i = 0; i < spt.size(); i++)
		{
			fout << spt[i].coor[0] << " " << spt[i].coor[1] << " " << spt[i].coor[2] << "\n";
		}
		fout << "\nCELLS " << mesh.size() << " " << 9 * mesh.size() << '\n';
		for (i = 0; i < mesh.size(); i++)
		{
			fout << "8 " << mesh[i].IEN[0] << " " << mesh[i].IEN[1] << " " << mesh[i].IEN[2] << " " << mesh[i].IEN[3]
			     << " " << mesh[i].IEN[4] << " " << mesh[i].IEN[5] << " " << mesh[i].IEN[6] << " " << mesh[i].IEN[7] << '\n';
		}
		fout << "\nCELL_TYPES " << mesh.size() << '\n';
		for (i = 0; i < mesh.size(); i++)
		{
			fout << "12\n";
		}
		fout << "POINT_DATA " << var.size() << "\nSCALARS AllParticles float 1\nLOOKUP_TABLE default\n";
		for (uint i = 0; i < var.size(); i++)
		{
			fout << var[i] /* +N_plus[i]+N_minus[i] */ << "\n";
			// fout << spt[i].label /* +N_plus[i]+N_minus[i] */ << "\n";
		}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}
}

void NeuronGrowth::ConcentrationCal_Coupling_Bezier(float u, float v, float w, const Element3D& bzel, float pt[3], float& disp, float dudx[3], float& detJ)
{
	float dUdx[3][3];
	vector<float> Nx(bzel.IEN.size());
	vector<array<float, 3>> dNdx(bzel.IEN.size());
	bzel.Para2Phys(u, v, w, pt);
	BasisFunction(u, v, w, bzel.pts, bzel.cmat, Nx, dNdx, dUdx, detJ);
	disp = 0.;
	dudx[0] = 0.; dudx[1] = 0.; dudx[2] = 0.;
	for (uint i = 0; i < bzel.IEN.size(); i++) {
		disp += Nx[i] * (N_0[bzel.IEN[i]]);
		dudx[0] += dNdx[i][0] * (N_0[bzel.IEN[i]]);
		dudx[1] += dNdx[i][1] * (N_0[bzel.IEN[i]]);
		dudx[2] += dNdx[i][2] * (N_0[bzel.IEN[i]]);
	}
}

void NeuronGrowth::VisualizeVTK_PhysicalDomain(int step, string var, string fn)
{
	vector<array<float, 3>> spt_all;//sample points
	vector<float> sresult_all;
	vector<array<int, 8>> sele_all;
	float detJ;
	int num_bzmesh_ele = bzmesh_process.size();
	float spt_proc[num_bzmesh_ele * 24];
	float sresult_proc[num_bzmesh_ele * 8];
	int sele_proc[num_bzmesh_ele * 8];
	for (unsigned int e = 0; e<num_bzmesh_ele; e++)
	{
		int ns(2);
		vector<float> su(ns);
		for (int i = 0; i < ns; i++)
		{
			su[i] = float(i) / (float(ns) - 1.);
		}

		int loc(0);
		for (int a = 0; a<ns; a++)
		{
			for (int b = 0; b<ns; b++)
			{
				for (int c = 0; c < ns; c++)
				{
					float pt1[3], dudx[3];
					float result;
					ConcentrationCal_Coupling_Bezier(su[c], su[b], su[a], bzmesh_process[e] , pt1, result, dudx, detJ);
					spt_proc[24 * e + loc * 3 + 0] = pt1[0];
					spt_proc[24 * e + loc * 3 + 1] = pt1[1];
					spt_proc[24 * e + loc * 3 + 2] = pt1[2];
					sresult_proc[8 * e + loc] = result;
					loc++;
				}
			}
		}
		int nns[2] = { ns*ns*ns, ns*ns };
		for (int a = 0; a<ns - 1; a++)
		{
			for (int b = 0; b<ns - 1; b++)
			{
				for (int c = 0; c < ns - 1; c++)
				{
					sele_proc[8 * e + 0] = 8 * e + a*nns[1] + b*ns + c;
					sele_proc[8 * e + 1] = 8 * e + a*nns[1] + b*ns + c + 1;
					sele_proc[8 * e + 2] = 8 * e + a*nns[1] + (b + 1)*ns + c + 1;
					sele_proc[8 * e + 3] = 8 * e + a*nns[1] + (b + 1)*ns + c;
					sele_proc[8 * e + 4] = 8 * e + (a + 1)*nns[1] + b*ns + c;
					sele_proc[8 * e + 5] = 8 * e + (a + 1)*nns[1] + b*ns + c + 1;
					sele_proc[8 * e + 6] = 8 * e + (a + 1)*nns[1] + (b + 1)*ns + c + 1;
					sele_proc[8 * e + 7] = 8 * e + (a + 1)*nns[1] + (b + 1)*ns + c;
				}
			}
		}
	}


	float *spts = NULL;
	float *sresults = NULL;
	int *seles = NULL;
	int *displs_spts = NULL;
	int *displs_sresults = NULL;
	int *displs_seles = NULL;
	int *num_bzmesh_eles = NULL;
	int *recvcounts_spts = NULL;
	int *recvcounts_sresults = NULL;
	int *recvcounts_seles = NULL;

	if (comRank == 0)
	{
		num_bzmesh_eles = (int*)malloc(sizeof(int)*nProcess);
		recvcounts_spts = (int*)malloc(sizeof(int)*nProcess);
		recvcounts_sresults = (int*)malloc(sizeof(int)*nProcess);
		recvcounts_seles = (int*)malloc(sizeof(int)*nProcess);
	}
	MPI_Gather(&num_bzmesh_ele, 1, MPI_INT, num_bzmesh_eles, 1, MPI_INT, 0, PETSC_COMM_WORLD);
	MPI_Barrier(comm);

	if (comRank == 0)
	{
		spts = (float*)malloc(sizeof(float) * 24 * n_bzmesh);
		sresults = (float*)malloc(sizeof(float) * 8* n_bzmesh);
		seles = (int*)malloc(sizeof(int) * 8 * n_bzmesh);

		displs_spts = (int*)malloc(nProcess * sizeof(int));
		displs_sresults = (int*)malloc(nProcess * sizeof(int));
		displs_seles = (int*)malloc(nProcess * sizeof(int));
		displs_spts[0] = 0;
		displs_sresults[0] = 0;
		displs_seles[0] = 0;

		for (int i = 1; i<nProcess; i++) {
			displs_spts[i] = displs_spts[i - 1] + num_bzmesh_eles[i - 1] * 24;
			displs_sresults[i] = displs_sresults[i - 1] + num_bzmesh_eles[i - 1] * 8;
			displs_seles[i] = displs_seles[i - 1] + num_bzmesh_eles[i - 1] * 8;
		}

		for (int i = 0; i < nProcess; i++)
		{
			recvcounts_spts[i] = num_bzmesh_eles[i] * 24;
			recvcounts_sresults[i] = num_bzmesh_eles[i] * 8;
			recvcounts_seles[i] = num_bzmesh_eles[i] * 8;
		}
	}

	MPI_Gatherv(spt_proc, num_bzmesh_ele * 8 * 3, MPI_FLOAT, spts, recvcounts_spts, displs_spts, MPI_FLOAT, 0, PETSC_COMM_WORLD);
	MPI_Gatherv(sresult_proc, num_bzmesh_ele * 8 , MPI_FLOAT, sresults, recvcounts_sresults, displs_sresults, MPI_FLOAT, 0, PETSC_COMM_WORLD);
	MPI_Gatherv(sele_proc, num_bzmesh_ele * 8, MPI_INT, seles, recvcounts_seles, displs_seles, MPI_INT, 0, PETSC_COMM_WORLD);


	if (comRank == 0)
	{
		for (int i = 0; i < n_bzmesh; i++)
		{
			for (int j = 0; j < 8; j++)
			{
				array<float, 3> pt = { spts[i * 24 + j * 3 + 0], spts[i * 24 + j * 3 + 1], spts[i * 24 + j * 3 + 2] };
				spt_all.push_back(pt);
				sresult_all.push_back(sresults[i * 8 + j]);

			}
		}
		int sum_ele = 0;
		int pstart = 0;
		for (int i = 0; i < nProcess; i++)
		{
			for (int e = 0; e < num_bzmesh_eles[i]; e++)
			{
				array<int, 8> el;
				el[0] = pstart + seles[8 * sum_ele + 0];
				el[1] = pstart + seles[8 * sum_ele + 1];
				el[2] = pstart + seles[8 * sum_ele + 2];
				el[3] = pstart + seles[8 * sum_ele + 3];
				el[4] = pstart + seles[8 * sum_ele + 4];
				el[5] = pstart + seles[8 * sum_ele + 5];
				el[6] = pstart + seles[8 * sum_ele + 6];
				el[7] = pstart + seles[8 * sum_ele + 7];
				sele_all.push_back(el);
				sum_ele++;
			}
			pstart = pstart + num_bzmesh_eles[i] * 8;
		}
		cout << "Visualizing in Physical Domain...\n";
		WriteVTK(spt_all, sresult_all, sele_all, step, fn);
	}
}

void NeuronGrowth::WriteVTK(const vector<array<float, 3>> spt, const vector<float> sdisp, const vector<array<int, 8>> sele, int step, string fn)
{
	stringstream ss;
	ss << step;
	string fname = fn + "/physics_allparticle_" + ss.str() + ".vtk";
	ofstream fout;
	fout.open(fname.c_str());
	unsigned int i;
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nHex test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << spt.size() << " float\n";
		for (i = 0; i<spt.size(); i++)
		{
			fout << spt[i][0] << " " << spt[i][1] << " " << spt[i][2] << "\n";
		}
		fout << "\nCELLS " << sele.size() << " " << 9 * sele.size() << '\n';
		for (i = 0; i<sele.size(); i++)
		{
			fout << "8 " << sele[i][0] << " " << sele[i][1] << " " << sele[i][2] << " " << sele[i][3]
				<< " " << sele[i][4] << " " << sele[i][5] << " " << sele[i][6] << " " << sele[i][7] << '\n';
		}
		fout << "\nCELL_TYPES " << sele.size() << '\n';
		for (i = 0; i<sele.size(); i++)
		{
			fout << "12\n";
		}
		fout << "POINT_DATA " << sdisp.size() << "\nSCALARS AllParticle float 1\nLOOKUP_TABLE default\n";
		for (uint i = 0; i<sdisp.size(); i++)
		{
			fout << sdisp[i] << "\n";
		}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}
}

void NeuronGrowth::CalculateVarsForOutput(vector<array<float, 3>> &spt_all, vector<float> &sresult_all, vector<array<int, 8>> &sele_all) 
{
	float detJ;
	int num_bzmesh_ele = bzmesh_process.size();
	float spt_proc[num_bzmesh_ele * 24];
	float sresult_proc[num_bzmesh_ele * 8];
	int sele_proc[num_bzmesh_ele * 8];
	for (unsigned int e = 0; e<num_bzmesh_ele; e++)
	{
		int ns(2);
		vector<float> su(ns);
		for (int i = 0; i < ns; i++)
		{
			su[i] = float(i) / (float(ns) - 1.);
		}

		int loc(0);
		for (int a = 0; a<ns; a++)
		{
			for (int b = 0; b<ns; b++)
			{
				for (int c = 0; c < ns; c++)
				{
					float pt1[3], dudx[3];
					float result;
					ConcentrationCal_Coupling_Bezier(su[c], su[b], su[a], bzmesh_process[e] , pt1, result, dudx, detJ);
					spt_proc[24 * e + loc * 3 + 0] = pt1[0];
					spt_proc[24 * e + loc * 3 + 1] = pt1[1];
					spt_proc[24 * e + loc * 3 + 2] = pt1[2];
					sresult_proc[8 * e + loc] = result;
					loc++;
				}
			}
		}
		int nns[2] = { ns*ns*ns, ns*ns };
		for (int a = 0; a<ns - 1; a++)
		{
			for (int b = 0; b<ns - 1; b++)
			{
				for (int c = 0; c < ns - 1; c++)
				{
					sele_proc[8 * e + 0] = 8 * e + a*nns[1] + b*ns + c;
					sele_proc[8 * e + 1] = 8 * e + a*nns[1] + b*ns + c + 1;
					sele_proc[8 * e + 2] = 8 * e + a*nns[1] + (b + 1)*ns + c + 1;
					sele_proc[8 * e + 3] = 8 * e + a*nns[1] + (b + 1)*ns + c;
					sele_proc[8 * e + 4] = 8 * e + (a + 1)*nns[1] + b*ns + c;
					sele_proc[8 * e + 5] = 8 * e + (a + 1)*nns[1] + b*ns + c + 1;
					sele_proc[8 * e + 6] = 8 * e + (a + 1)*nns[1] + (b + 1)*ns + c + 1;
					sele_proc[8 * e + 7] = 8 * e + (a + 1)*nns[1] + (b + 1)*ns + c;
				}
			}
		}
	}


	float *spts = NULL;
	float *sresults = NULL;
	int *seles = NULL;
	int *displs_spts = NULL;
	int *displs_sresults = NULL;
	int *displs_seles = NULL;
	int *num_bzmesh_eles = NULL;
	int *recvcounts_spts = NULL;
	int *recvcounts_sresults = NULL;
	int *recvcounts_seles = NULL;

	if (comRank == 0)
	{
		num_bzmesh_eles = (int*)malloc(sizeof(int)*nProcess);
		recvcounts_spts = (int*)malloc(sizeof(int)*nProcess);
		recvcounts_sresults = (int*)malloc(sizeof(int)*nProcess);
		recvcounts_seles = (int*)malloc(sizeof(int)*nProcess);
	}
	MPI_Gather(&num_bzmesh_ele, 1, MPI_INT, num_bzmesh_eles, 1, MPI_INT, 0, PETSC_COMM_WORLD);
	MPI_Barrier(comm);

	if (comRank == 0)
	{
		spts = (float*)malloc(sizeof(float) * 24 * n_bzmesh);
		sresults = (float*)malloc(sizeof(float) * 8* n_bzmesh);
		seles = (int*)malloc(sizeof(int) * 8 * n_bzmesh);

		displs_spts = (int*)malloc(nProcess * sizeof(int));
		displs_sresults = (int*)malloc(nProcess * sizeof(int));
		displs_seles = (int*)malloc(nProcess * sizeof(int));
		displs_spts[0] = 0;
		displs_sresults[0] = 0;
		displs_seles[0] = 0;

		for (int i = 1; i<nProcess; i++) {
			displs_spts[i] = displs_spts[i - 1] + num_bzmesh_eles[i - 1] * 24;
			displs_sresults[i] = displs_sresults[i - 1] + num_bzmesh_eles[i - 1] * 8;
			displs_seles[i] = displs_seles[i - 1] + num_bzmesh_eles[i - 1] * 8;
		}

		for (int i = 0; i < nProcess; i++)
		{
			recvcounts_spts[i] = num_bzmesh_eles[i] * 24;
			recvcounts_sresults[i] = num_bzmesh_eles[i] * 8;
			recvcounts_seles[i] = num_bzmesh_eles[i] * 8;
		}
	}

	MPI_Gatherv(spt_proc, num_bzmesh_ele * 8 * 3, MPI_FLOAT, spts, recvcounts_spts, displs_spts, MPI_FLOAT, 0, PETSC_COMM_WORLD);
	MPI_Gatherv(sresult_proc, num_bzmesh_ele * 8 , MPI_FLOAT, sresults, recvcounts_sresults, displs_sresults, MPI_FLOAT, 0, PETSC_COMM_WORLD);
	MPI_Gatherv(sele_proc, num_bzmesh_ele * 8, MPI_INT, seles, recvcounts_seles, displs_seles, MPI_INT, 0, PETSC_COMM_WORLD);


	if (comRank == 0)
	{
		for (int i = 0; i < n_bzmesh; i++)
		{
			for (int j = 0; j < 8; j++)
			{
				array<float, 3> pt = { spts[i * 24 + j * 3 + 0], spts[i * 24 + j * 3 + 1], spts[i * 24 + j * 3 + 2] };
				spt_all.push_back(pt);
				sresult_all.push_back(sresults[i * 8 + j]);

			}
		}
		int sum_ele = 0;
		int pstart = 0;
		for (int i = 0; i < nProcess; i++)
		{
			for (int e = 0; e < num_bzmesh_eles[i]; e++)
			{
				array<int, 8> el;
				el[0] = pstart + seles[8 * sum_ele + 0];
				el[1] = pstart + seles[8 * sum_ele + 1];
				el[2] = pstart + seles[8 * sum_ele + 2];
				el[3] = pstart + seles[8 * sum_ele + 3];
				el[4] = pstart + seles[8 * sum_ele + 4];
				el[5] = pstart + seles[8 * sum_ele + 5];
				el[6] = pstart + seles[8 * sum_ele + 6];
				el[7] = pstart + seles[8 * sum_ele + 7];
				sele_all.push_back(el);
				sum_ele++;
			}
			pstart = pstart + num_bzmesh_eles[i] * 8;
		}
	}
}

/**
 * Visualizes the physical domain of neuron growth for all variables using VTK.
 *
 * @param step The current simulation step (used for time-stepping control).
 * @param fn The filename prefix for the output VTK files.
 */
void NeuronGrowth::VisualizeVTK_PhysicalDomain_All(int step, string fn) {
	// Initialize containers for sample points, simulation results, and element connectivity for all variables
	vector<vector<array<float, 3>>> spt_all_4var(4); // Sample points for all 4 variables
	vector<vector<float>> sresult_all_4var(4); // Simulation results for all 4 variables
	vector<vector<array<int, 8>>> sele_all_4var(4); // Element connectivity for all 4 variables

	// Set current variable to phi and calculate its output variables
	N_0 = phi;
	CalculateVarsForOutput(spt_all_4var[0], sresult_all_4var[0], sele_all_4var[0]);

	// Set current variable to syn and calculate its output variables
	N_0 = syn;
	CalculateVarsForOutput(spt_all_4var[1], sresult_all_4var[1], sele_all_4var[1]);

	// Set current variable to tub and calculate its output variables
	N_0 = tub;
	CalculateVarsForOutput(spt_all_4var[2], sresult_all_4var[2], sele_all_4var[2]);

	// Set current variable to tips and calculate its output variables
	N_0 = tips;
	CalculateVarsForOutput(spt_all_4var[3], sresult_all_4var[3], sele_all_4var[3]);

	// If running on the master process (comRank == 0), write the VTK file for visualization
	if (comRank == 0) {
		WriteVTK_ALL(spt_all_4var[0], sresult_all_4var, sele_all_4var[0], step, fn);
	}
}

/**
 * Writes the simulation data to a VTK file for visualization.
 *
 * @param spt Sample points' coordinates.
 * @param sdisp Scalar displacements for each variable at the sample points.
 * @param sele Connectivity information for the cells (elements).
 * @param step The current simulation step, used for generating the filename.
 * @param fn The base filename to which the step number and extension are appended.
 */
void NeuronGrowth::WriteVTK_ALL(const vector<array<float, 3>> spt, const vector<vector<float>> sdisp, const vector<array<int, 8>> sele, int step, string fn) {
	stringstream ss;
	ss << step;
	string fname = fn + "/physics_allparticle_" + ss.str() + ".vtk";
	ofstream fout(fname.c_str());

	if (fout.is_open()) {
		// Write the VTK file header
		fout << "# vtk DataFile Version 2.0\nHex test\nASCII\nDATASET UNSTRUCTURED_GRID\n";

		// Write the points
		fout << "POINTS " << spt.size() << " float\n";
		for (unsigned int i = 0; i < spt.size(); i++) {
			fout << spt[i][0] << " " << spt[i][1] << " " << spt[i][2] << "\n";
		}

		// Write the cells
		fout << "\nCELLS " << sele.size() << " " << 9 * sele.size() << '\n';
		for (unsigned int i = 0; i < sele.size(); i++) {
			fout << "8 " << sele[i][0] << " " << sele[i][1] << " " << sele[i][2] << " " << sele[i][3]
				<< " " << sele[i][4] << " " << sele[i][5] << " " << sele[i][6] << " " << sele[i][7] << '\n';
		}

		// Write the cell types
		fout << "\nCELL_TYPES " << sele.size() << '\n';
		for (unsigned int i = 0; i < sele.size(); i++) {
			fout << "12\n"; // 12 corresponds to VTK_HEXAHEDRON
		}

		// Write scalar fields
		fout << "POINT_DATA " << sdisp[0].size() << "\n";
		const char* scalarNames[] = {"phi", "synaptogenesis", "tubulin", "tips"};
		for (size_t varIndex = 0; varIndex < sdisp.size(); ++varIndex) {
			fout << "\nSCALARS " << scalarNames[varIndex] << " float 1\nLOOKUP_TABLE default\n";
			for (size_t i = 0; i < sdisp[varIndex].size(); i++) {
				fout << sdisp[varIndex][i] << "\n";
			}
		}

		fout.close();
	} else {
		cout << "Cannot open " << fname << "!\n";
	}
}

void NeuronGrowth::PointFormValue(vector<float> &Nx, const vector<float> &U, float Value)
{
	Value = 0;

	for (int j = 0; j < Nx.size(); j++)
		Value += U[j] * Nx[j];
}

void NeuronGrowth::PointFormGrad(vector<array<float, 3>> &dNdx, const vector<float> &U, float Value[2])
{
	for (int j = 0; j < 2; j++)
		Value[j] = 0.;

	for (int j = 0; j < 2; j++)
		for (int k = 0; k < dNdx.size(); k++)
			Value[j] += U[k] * dNdx[k][j];
}

void NeuronGrowth::PointFormHess(vector<array<array<float, 2>, 2>> &d2Ndx2, const vector<float> &U, float Value[2][2])
{
	for (int j = 0; j < 2; j++)
	{
		for (int k = 0; k < 2; k++)
		{
			Value[j][k] = 0.;
		}
	}
	for (int j = 0; j < 2; j++)
	{
		for (int k = 0; k < 2; k++)
		{
			for (int l = 0; l < d2Ndx2.size(); l++)
			{
				Value[j][k] += U[l] * d2Ndx2[l][j][k];
			}
		}
	}
}

void NeuronGrowth::ElementValue(const vector<float> &Nx, const vector<float> value_node, float &value)
{
	value = 0.;
	for (int i = 0; i < Nx.size(); i++)
		value += value_node[i] * Nx[i];
}

void NeuronGrowth::ElementValueAll(const vector<float> &Nx, const vector<float> elePhiGuess, float &elePG, const vector<float> elePhi, float &eleP, const vector<float> eleSyn, float &eleS, const vector<float> eleTips, float &eleTp, const vector<float> eleTubulin, float &eleTb, const vector<float> eleEpsilon, float &eleEP, const vector<float> eleEpsilonP, float &eleEEP)
{
	elePG = 0.;
	eleP = 0.;
	eleS = 0.;
	eleTp = 0.;
	eleTb = 0.;
	eleEP = 0.;
	eleEEP = 0.;
	for (int i = 0; i < Nx.size(); i++) {
		elePG += elePhiGuess[i] * Nx[i];
		eleP += elePhi[i] * Nx[i];
		eleS += eleSyn[i] * Nx[i];
		eleTp += eleTips[i] * Nx[i];
		eleTb += eleTubulin[i] * Nx[i];
		eleEP += eleEpsilon[i] * Nx[i];
		eleEEP += eleEpsilonP[i] * Nx[i];
	}
}

void NeuronGrowth::ElementDeriv(const int nen, vector<array<float, 3>> &dNdx, const vector<float> value_node, float &dVdx, float &dVdy, float &dVdz)
{
	dVdx = 0;
	dVdy = 0.;
	dVdz = 0.;
	for (int i = 0; i < nen; i++) {
		dVdx += value_node[i] * dNdx[i][0];
		dVdy += value_node[i] * dNdx[i][1];
		dVdz += value_node[i] * dNdx[i][2];
	}
}


void NeuronGrowth::ElementDerivAll(const int nen, vector<array<float, 3>> &dNdx, const vector<float> elePhiGuess, float &dPGdx, float &dPGdy, const vector<float> eleTheta, float &dThedx, float &dThedy, const vector<float> eleEpsilon, float &dAdx, float &dAdy, const vector<float> eleEpsilonP, float &dAPdx, float &dAPdy)
{	
	dPGdx = 0; dPGdy = 0.;
	dThedx = 0; dThedy = 0.;
	dAdx = 0; dAdy = 0.;
	dAPdx = 0; dAPdy = 0.;
	for (int i = 0; i < nen; i++) {
		dPGdx += elePhiGuess[i] * dNdx[i][0];	dPGdy += elePhiGuess[i] * dNdx[i][1];
		dThedx += eleTheta[i] * dNdx[i][0];	dThedy += eleTheta[i] * dNdx[i][1];
		dAdx += eleEpsilon[i] * dNdx[i][0];	dAdy += eleEpsilon[i] * dNdx[i][1];
		dAPdx += eleEpsilonP[i] * dNdx[i][0];	dAPdy += eleEpsilonP[i] * dNdx[i][1];
	}
}

void NeuronGrowth::ElementEvaluationAll_phi(const int nen, const vector<float> &Nx, vector<array<float, 3>> &dNdx, const vector<float> elePhiGuess, float &elePG, const vector<float> elePhi, float &eleP, const vector<float> eleSyn, float &eleS, const vector<float> eleTips, float &eleTp, const vector<float> eleTubulin, float &eleTb, float &dPGdx, float &dPGdy, float &dPGdz)
{
	elePG = 0.;
	eleP = 0.;
	eleS = 0.;
	eleTp = 0.;
	eleTb = 0.;

	dPGdx = 0; dPGdy = 0, dPGdz = 0;

	for (int i = 0; i < nen; i++) {
		elePG += elePhiGuess[i] * Nx[i];
		eleP += elePhi[i] * Nx[i];
		eleS += eleSyn[i] * Nx[i];
		eleTp += eleTips[i] * Nx[i];
		eleTb += eleTubulin[i] * Nx[i];

		dPGdx += elePhiGuess[i] * dNdx[i][0];	dPGdy += elePhiGuess[i] * dNdx[i][1];	dPGdz += elePhiGuess[i] * dNdx[i][2];
	}
}

void NeuronGrowth::ElementEvaluationAll_phi(const int nen, const vector<float> &Nx, vector<array<float, 3>> &dNdx, const vector<float> elePhiGuess, float &elePG, const vector<float> elePhi, float &eleP, const vector<float> eleSyn, float &eleS, const vector<float> eleTips, float &eleTp, const vector<float> eleTubulin, float &eleTb, const vector<float> eleEpsilon, float &eleEP, const vector<float> eleEpsilonP, float &eleEEP, float &dPGdx, float &dPGdy, float &dAdx, float &dAdy, float &dAPdx, float &dAPdy)
{
	elePG = 0.;
	eleP = 0.;
	eleS = 0.;
	eleTp = 0.;
	eleTb = 0.;
	eleEP = 0.;
	eleEEP = 0.;

	dPGdx = 0; dPGdy = 0.;
	dAdx = 0; dAdy = 0.;
	dAPdx = 0; dAPdy = 0.;

	for (int i = 0; i < nen; i++) {
		elePG += elePhiGuess[i] * Nx[i];
		eleP += elePhi[i] * Nx[i];
		eleS += eleSyn[i] * Nx[i];
		eleTp += eleTips[i] * Nx[i];
		eleTb += eleTubulin[i] * Nx[i];
		eleEP += eleEpsilon[i] * Nx[i];
		eleEEP += eleEpsilonP[i] * Nx[i];

		dPGdx += elePhiGuess[i] * dNdx[i][0];	dPGdy += elePhiGuess[i] * dNdx[i][1];
		dAdx += eleEpsilon[i] * dNdx[i][0];	dAdy += eleEpsilon[i] * dNdx[i][1];
		dAPdx += eleEpsilonP[i] * dNdx[i][0];	dAPdy += eleEpsilonP[i] * dNdx[i][1];
	}
}

void NeuronGrowth::ElementEvaluationAll_phi(const int nen, const vector<float> &Nx, vector<array<float, 3>> &dNdx, vector<vector<float>> &eleVal, vector<float> &vars)
{
	// eleVal Array Description
	// Index | Description
	// ------|----------------
	// 0     | elePhiGuess  - guess for elePhi
	// 1     | elePhi       - Phi value
	// 2     | eleSyn       - Syn value
	// 3     | eleTubulin   - Tubulin value
	// 4     | eleTheta     - Theta value
	// 5     | eleTips      - Tips value
	// 6     | eleEpsilon   - Epsilon value
	// 7     | eleEpsilonP  - Epsilon derivative value

	// Vars Array Description
	// Index | Variable | Description
	// ------|----------|-------------------------------------------------
	// 0     | C1       | Coefficient 1
	// 1     | C0       | Coefficient 0
	// 2     | elePG    | Phi guess element
	// 3     | dPGdx    | Derivative of PG with respect to x
	// 4     | dPGdy    | Derivative of PG with respect to y
	// 5     | eleP     | Phi element
	// 6     | eleS     | Syn element
	// 7     | eleTb    | Tubulin element
	// 8     | eleE     | Energy element
	// 9     | eleTp    | Tip element
	// 10    | dThedx   | Derivative of Theta with respect to x
	// 11    | dThedy   | Derivative of Theta with respect to y
	// 12    | eleEP    | Epsilon element
	// 13    | dAdx     | Derivative of A with respect to x
	// 14    | dAdy     | Derivative of A with respect to y
	// 15    | eleEEP   | Epsilon derivative element
	// 16    | dAPdx    | Derivative of AP with respect to x
	// 17    | dAPdy    | Derivative of AP with respect to y
	// 18    | dPGdz    | Derivative of PG with respect to z
	// 19    | dAdz     | Derivative of A with respect to z
	// 20    | dAPdz    | Derivative of AP with respect to z

	vars[2] = 0.;
	vars[5] = 0.;
	vars[6] = 0.;
	vars[9] = 0.;
	vars[7] = 0.;
	vars[12] = 0.;
	vars[15] = 0.;

	vars[3] = 0; vars[4] = 0.; vars[18] = 0.;
	vars[13] = 0; vars[14] = 0.; vars[19] = 0.;
	vars[16] = 0; vars[17] = 0.; vars[20] = 0.;

	for (int i = 0; i < nen; i++) {
		vars[2] += eleVal[0][i] * Nx[i];
		vars[5] += eleVal[1][i] * Nx[i];
		vars[6] += eleVal[2][i] * Nx[i];
		vars[9] += eleVal[5][i] * Nx[i];
		vars[7] += eleVal[3][i] * Nx[i];
		// vars[12] += eleVal[6][i] * Nx[i];
		// vars[15] += eleVal[7][i] * Nx[i];

		vars[3] += eleVal[0][i] * dNdx[i][0];	vars[4] += eleVal[0][i] * dNdx[i][1];	vars[18] += eleVal[0][i] * dNdx[i][2];
		// vars[13] += eleVal[6][i] * dNdx[i][0];	vars[14] += eleVal[6][i] * dNdx[i][1];	vars[19] += eleVal[6][i] * dNdx[i][2];
		// vars[16] += eleVal[7][i] * dNdx[i][0];	vars[17] += eleVal[7][i] * dNdx[i][1];	vars[20] += eleVal[7][i] * dNdx[i][2];
	}
}

void NeuronGrowth::ElementEvaluationAll_syn_tub(const int nen, const vector<float> &Nx, vector<array<float, 3>> &dNdx, const vector<float> elePhiDiff, float &elePf, const vector<float> eleSyn, float &eleS, const vector<float> elePhi, float &eleP, const vector<float> elePhiPrev, float &elePprev, const vector<float> eleConct, float &eleC, float &dPdx, float &dPdy, float &dPdz)
{
	elePf = 0.;
	eleS = 0.;
	eleP = 0.;
	elePprev = 0.;
	eleC = 0.;

	dPdx = 0; dPdy = 0.; dPdy = 0.;

	for (int i = 0; i < nen; i++) {
		elePf += elePhiDiff[i] * Nx[i];
		eleS += eleSyn[i] * Nx[i];

		eleP += elePhi[i] * Nx[i];
		elePprev += elePhiPrev[i] * Nx[i];
		eleC += eleConct[i] * Nx[i];

		dPdx += elePhi[i] * dNdx[i][0];	dPdy += elePhi[i] * dNdx[i][1];	dPdz += elePhi[i] * dNdx[i][2];
	}
}

void NeuronGrowth::ElementEvaluationAll_syn_tub(const int nen, const vector<float> &Nx, vector<array<float, 3>> &dNdx, vector<vector<float>> &eleVal, vector<float> &vars)
{
	// Element Value Definitions (eleVal)
	// eleVal[0]: elePhiDiff   - Difference in Phi
	// eleVal[1]: eleSyn       - Synthesis element
	// eleVal[2]: elePhi       - Current Phi
	// eleVal[3]: elePhiPrev   - Previous Phi
	// eleVal[4]: eleConct     - Conductivity

	// Variable Definitions
	// Index  Variable         Description
	// -----  ---------------  -------------------------------------------
	// 0      elePf            - phi difference
	// 1      eleS             - synaptogenesis
	// 2      eleP             - current phi
	// 3      dPdx             - Partial derivative of phi with respect to x
	// 4      dPdy             - Partial derivative of phi with respect to y
	// 5      elePprev         - previous phi
	// 6      eleC             - tubulin concentration
	// 7      mag_grad_phi0    - Magnitude of gradient of phi0
	// 8      term_diff        - Diffusion term
	// 9      term_alph        - Alpha term (active transport)
	// 10     term_beta        - Beta term (degradation)
	// 11     term_source      - Source term (source generation)
	// 12     dPdz             - Partial derivative of P with respect to z

	vars[0] = 0.;
	vars[1] = 0.;
	vars[2] = 0.;
	vars[5] = 0.;
	vars[6] = 0.;

	vars[3] = 0; vars[4] = 0.; vars[12] = 0;

	for (int i = 0; i < nen; i++) {
		vars[0] += eleVal[0][i] * Nx[i];
		vars[1] += eleVal[1][i] * Nx[i];

		vars[2] += eleVal[2][i] * Nx[i];
		vars[5] += eleVal[3][i] * Nx[i];
		vars[6] += eleVal[4][i] * Nx[i];

		vars[3] += eleVal[2][i] * dNdx[i][0];	vars[4] += eleVal[2][i] * dNdx[i][1]; 	vars[12] += eleVal[2][i] * dNdx[i][2];
	}
}

void NeuronGrowth::prepareBasis() {

	/*Build linear system in each process*/
	float detJ;
	float dudx[3][3];
	
	float dThedx, dThedy, dThedz;
	sum_grad_phi0_local = 0;
	int e;

	// vector<float> elePolar, eleAzimuth;
	// float elePo, eleAz;

	for (e = 0; e < bzmesh_process.size(); e++) {
		int nen = bzmesh_process[e].IEN.size();

		elePhi0.resize(nen);
		eleTheta.resize(nen);
		// elePolar.resize(nen);
		// eleAzimuth.resize(nen);
		
		for (int i = 0; i < nen; i++) {
			elePhi0[i] = phi_0[bzmesh_process[e].IEN[i]];
			eleTheta[i] = theta[bzmesh_process[e].IEN[i]];
			// elePolar[i] = polar[bzmesh_process[e].IEN[i]];
			// eleAzimuth[i] = azimuth[bzmesh_process[e].IEN[i]];
		}

		for (int i = 0; i < Gpt.size(); i++) {
			for (int j = 0; j < Gpt.size(); j++) {
				for (int k = 0; k < Gpt.size(); k++) {
					BasisFunction(Gpt[i], Gpt[j], Gpt[k], bzmesh_process[e].pts, bzmesh_process[e].cmat, Nx, dNdx, dudx, detJ);
					detJ = wght[i] * wght[j] * wght[k] * detJ;

					pre_Nx.push_back(Nx);
					pre_dNdx.push_back(dNdx);
					pre_detJ.push_back(detJ);

					ElementDeriv(nen, dNdx, eleTheta, dThedx, dThedy, dThedz);
					pre_C0.push_back((0.5 - 6 * s_coeff * sqrt(pow(dThedx, 2) + pow(dThedy, 2) + pow(dThedz, 2))));
					// pre_C0.push_back((6 * s_coeff * sqrt(pow(dThedx, 2) + pow(dThedy, 2) + pow(dThedz, 2))));

					// ElementValue(Nx, elePolar, elePo);
					// ElementValue(Nx, eleAzimuth, eleAz);
					// pre_C0_sp.push_back((0.5 + 6 * s_coeff * sqrt(pow(elePo, 2) + pow(eleAz, 2))));

					ElementDeriv(nen, dNdx, elePhi0, dP0dx, dP0dy, dP0dz);
					pre_mag_grad_phi0.push_back(pow(dP0dx, 2) + pow(dP0dy, 2) + pow(dP0dz, 2));
					sum_grad_phi0_local += pow(dP0dx, 2) + pow(dP0dy, 2) + pow(dP0dz, 2);
				}
			}
		}
	}
	MPI_Barrier(PETSC_COMM_WORLD);
}

// void NeuronGrowth::preparePhaseField() 
// {
//     // Clear previous data
//     pre_eleEP.clear();
//     pre_eleEEP.clear();
//     pre_dAdx.clear();
//     pre_dAdy.clear();
//     pre_dAdz.clear();  // Added for 3D
//     pre_eleP.clear();
//     pre_eleTh.clear();
//     pre_eleMp.clear();
//     pre_C1.clear();

//     // Initialize variables
//     float eleEP(0), eleEEP(0), dAdx(0), dAdy(0), dAdz(0), // Added dAdz
//           eleP(0), eleTh(0), eleS(0), eleTb(0), eleTp(0), eleE(0);
//     uint ind(0);

//     // Loop over elements
//     for (size_t e = 0; e < bzmesh_process.size(); e++) {
//         size_t nen = bzmesh_process[e].IEN.size();

//         // Initialize element-wise variables
//         vector<float> elePhi(nen, 0), eleTheta(nen, 0), eleEpsilon(nen, 0), eleEpsilonP(nen, 0);
//         vector<float> eleSyn(nen, 0), eleTubulin(nen, 0), eleTips(nen, 0);

//         // Extract nodal values for the current element
//         for (size_t i = 0; i < nen; i++) {
//             size_t nodeIndex = bzmesh_process[e].IEN[i];
//             elePhi[i]    = phi[nodeIndex];       // elePhi
//             eleTheta[i]  = theta[nodeIndex];     // eleTheta
//             eleSyn[i]    = syn[nodeIndex];       // eleSyn
//             eleTubulin[i]= tub[nodeIndex];       // eleTubulin
//             eleTips[i]   = tips[nodeIndex];      // eleTips
//         }

//         // Loop over Gaussian quadrature points in 3D
//         for (size_t i = 0; i < Gpt.size(); i++) {
//             for (size_t j = 0; j < Gpt.size(); j++) {
//                 for (size_t k = 0; k < Gpt.size(); k++) { // Added k-loop for 3D

//                     // Evaluate orientation and compute element variables
//                     EvaluateOrientation(nen, pre_Nx[ind], pre_dNdx[ind], elePhi, eleTheta, eleEpsilon, eleEpsilonP);
// 					EvaluateOrientation(nen, pre_Nx[ind], pre_dNdx[ind], elePhi, eleTheta, eleAniso, dA_dPdx, dA_dPdy, dA_dPdz)

//                     // Compute element values and derivatives
//                     ElementValue(pre_Nx[ind], eleEpsilon, eleEP);
//                     pre_eleEP.push_back(eleEP);

//                     ElementValue(pre_Nx[ind], eleEpsilonP, eleEEP);
//                     pre_eleEEP.push_back(eleEEP);

//                     ElementDeriv(nen, pre_dNdx[ind], eleEpsilonP, dAdx, dAdy, dAdz); // Modified to include dAdz
//                     pre_dAdx.push_back(dAdx);
//                     pre_dAdy.push_back(dAdy);
//                     pre_dAdz.push_back(dAdz); // Added for 3D

//                     ElementValue(pre_Nx[ind], elePhi, eleP);
//                     pre_eleP.push_back(eleP);

//                     ElementValue(pre_Nx[ind], eleTheta, eleTh);
//                     pre_eleTh.push_back(eleTh);

//                     ElementValue(pre_Nx[ind], eleSyn, eleS);
//                     ElementValue(pre_Nx[ind], eleTubulin, eleTb);
//                     ElementValue(pre_Nx[ind], eleTips, eleTp);

//                     // Adjust assembly and disassembly rates based on detected tips
//                     if (n < 50) { // Assuming 'n' is a time step or iteration count; ensure it's properly defined
//                         eleE = alphaOverPi * atan(gamma * (c_opt - eleS));
//                         pre_eleMp.push_back(M_neurite);
//                     } else {
//                         if (eleTp != 0) {
//                             eleE = alphaOverPi * atan(gamma * Regular_Heiviside_fun(50 * eleTb) * (c_opt - eleS));
//                             if (eleTp < 0) {
//                                 pre_eleMp.push_back(M_axon);
//                             } else {
//                                 pre_eleMp.push_back(M_neurite);
//                             }
//                         } else {
//                             eleE = alphaOverPi * atan(gamma * Regular_Heiviside_fun(r * eleTb - g) * (c_opt - eleS));
//                             pre_eleMp.push_back(5);
//                         }
//                     }

//                     // Calculate C1 variable for the phase-field energy term
//                     float C1 = eleE - pre_C0[ind]; // Update pre_C0[ind] if necessary to include 3D effects
//                     pre_C1.push_back(C1);

//                     // Increment the index for precomputed variables
//                     ind += 1;
//                 }
//             }
//         }
//     }
// }

void NeuronGrowth::prepareTerm_source() {
	int ind(0);
	pre_term_source.clear();
	for (int e = 0; e < bzmesh_process.size(); e++) {
		for (int i = 0; i < Gpt.size(); i++) {
			for (int j = 0; j < Gpt.size(); j++) {
				for (int k = 0; k < Gpt.size(); k++) {
					pre_term_source.push_back(source_coeff * pre_mag_grad_phi0[ind] / sum_grad_phi0_global);
					ind += 1;
				}
			}
		}
	}
}

// void NeuronGrowth::prepareEE()
// {
// 	/*Build linear system in each process*/
// 	PetscInt ind(0), e;
// 	for (e = 0; e < bzmesh_process.size(); e++) {
// 		PetscInt nen = bzmesh_process[e].IEN.size();

// 		for (PetscInt i = 0; i < nen * 1; i++) {
// 			EVectorSolve[i] = 0.0;
// 			for (PetscInt j = 0; j < nen * 1; j++) {
// 				EMatrixSolve[i][j] = 0.0;
// 			}
// 		}
	
// 		for (PetscInt i = 0; i < nen * 1; i++) {
// 			eleVal[0][i] = phi[bzmesh_process[e].IEN[i]];		// elePhiGuess
// 			// eleVal[1][i] = phi[bzmesh_process[e].IEN[i]];		// elePhi
// 			eleVal[0][i] = 0;		// elePhiGuess
// 			eleVal[1][i] = 0;		// elePhi

// 			eleVal[2][i] = syn[bzmesh_process[e].IEN[i]];		// eleSyn
// 			eleVal[3][i] = tub[bzmesh_process[e].IEN[i]];		// eleTubulin
// 			eleVal[4][i] = theta[bzmesh_process[e].IEN[i]];	// eleTheta
// 			eleVal[5][i] = tips[bzmesh_process[e].IEN[i]];	// eleTips
// 			eleVal[6][i] = 0;							// eleEpsilon
// 			eleVal[7][i] = 0;							// eleEpsilonP
// 		}

// 		for (PetscInt i = 0; i < Gpt.size(); i++) {
// 			for (PetscInt j = 0; j < Gpt.size(); j++) {

// 				PetscReal detJ;

// 				vector<float> Nx;
// 				vector<array<float, 3>> dNdx;
// 				Nx.clear();
// 				dNdx.clear();
// 				Nx.resize(nen);
// 				dNdx.resize(nen);

// 				// use pre-calculated values (save computational cost)
// 				Nx = pre_Nx[ind];
// 				dNdx = pre_dNdx[ind];
// 				detJ = pre_detJ[ind];
// 				vars[1] = pre_C0[ind]; // C0
// 				ind += 1;

// 				EvaluateOrientation(nen, Nx, dNdx, eleVal[1], eleVal[4], eleVal[6], eleVal[7]);
// 				ElementEvaluationAll_phi(nen, Nx, dNdx, eleVal, vars);
// 				// Vars definition
// 				//	 0   1    2     3      4      5     6     7      8     9      10      11      12     13    14    15      16     17   	18     19   20  
// 				// float C1, C0, elePG, dPGdx, dPGdy, eleP, eleS, eleTb, eleE, eleTp, dThedx, dThedy, eleEP, dAdx, dAdy, eleEEP, dAPdx, dAPdy, dPGdz, dAdz, dAPdz
				
// 				vars[8] = alphaOverPi*atan(gamma * (1 - vars[6]));

// 				vars[0] = vars[8] - vars[1]; // E - (0.5 + 6 * s_coeff * sqrt(pow(dThedx, 2) + pow(dThedy, 2));
									
// 				for (PetscInt m = 0; m < nen; m++) {
// 						EVectorSolve[m] += (vars[2] * Nx[m] - dt * M_phi *\
// 						((- vars[12] * vars[12] * (vars[3] * dNdx[m][0] + vars[4] * dNdx[m][1])) -\
// 						(- vars[13] * vars[15] * vars[4] * dNdx[m][0]) +\
// 						(- vars[14] * vars[15] * vars[3] * dNdx[m][1]) +\
// 						(- vars[2] * vars[2] * vars[2] + (1 - vars[0]) * vars[2] * vars[2] + vars[0] * vars[2]) * Nx[m])
// 						- vars[5] * Nx[m]) * detJ;
// 					for (PetscInt n = 0; n < nen; n++) {
// 						EMatrixSolve[m][n] += (Nx[m] * Nx[n] - dt * M_phi *
// 							((- vars[12] * vars[12] * (dNdx[m][0] * dNdx[n][0] + dNdx[m][1] * dNdx[n][1])) // terma2
// 							- (- vars[13] * vars[15] * dNdx[m][1] * dNdx[n][0]) // termadx
// 							+ (- vars[14] * vars[15] * dNdx[m][0] * dNdx[n][1]) // termady
// 							+ (- 3 * vars[2] * vars[2] + 2 * (1 - vars[0]) * vars[2] + vars[0] * Nx[m]) * Nx[n]) // termdbl
// 							) * detJ;
// 					}
// 				}
// 			}
// 		}

// 		/*Apply Boundary Condition*/
// 		for (PetscInt i = 0; i < nen; i++) {
// 			PetscInt A = bzmesh_process[e].IEN[i];
// 			if (cpts[A].label == 1) // domain boundary
// 				ApplyBoundaryCondition(0, i, 0, EMatrixSolve, EVectorSolve);
// 		}
		
// 		pre_EMatrixSolve.push_back(EMatrixSolve);
// 		pre_EVectorSolve.push_back(EVectorSolve);
// 	}
// }

void NeuronGrowth::EvaluateEnergy(const int nen, const vector<float> &Nx, const vector<float> eleSyn, vector<float>& E)
{
	for (int i = 0; i < nen; i++) {
		E[i] = alphaOverPi*atan(gamma*(1-eleSyn[i]));
	}
}

float NeuronGrowth::Regular_Heiviside_fun(float x) 
{
	// float epsilon = 0.00001; // the number is not fixed.
	float epsilon = 0.00001; // the number is not fixed.
	return 0.5*(1+(2/PI)*atan(x/epsilon));
}

/**
 * EvaluateOrientation - Calculates the orientation and its derivatives for neuronal growth.
 * 
 * This function evaluates the anisotropic properties and their spatial derivatives
 * based on the phase field (elePhi) and orientation (eleTheta) of elements. It is
 * designed to support the modeling of neuron growth by determining the directional
 * properties of neurite extension.
 * 
 * Parameters:
 * - nen: The number of elements (int) to consider in the computation.
 * - Nx: A vector of shape functions (vector<float>) evaluated at the integration points.
 * - dNdx: A vector of derivatives of shape functions (vector<array<float, 3>>) with respect to x, y, z.
 * - elePhi: A vector<float> representing the phase field values of elements.
 * - eleTheta: A vector<float> representing the orientation angles of elements.
 * - eleAniso: A reference to a float to store the computed anisotropy of the element.
 * - dA_dPdx: A reference to a float to store the derivative of anisotropy with respect to x.
 * - dA_dPdy: A reference to a float to store the derivative of anisotropy with respect to y.
 * - dA_dPdz: A reference to a float to store the derivative of anisotropy with respect to z.
 * 
 * The function iterates over each element, calculating the anisotropy (eleAniso) and its
 * derivatives (dA_dPdx, dA_dPdy, dA_dPdz) based on the phase field and orientation data.
 * It applies a condition to filter elements based on their phase field value, ensuring
 * computations are only performed on relevant elements. The function utilizes fourth and
 * third powers of the phase field derivatives to compute the anisotropy and its derivatives,
 * incorporating stability checks to handle computational anomalies.
 * 
 * Outputs:
 * - The function directly modifies eleAniso, dA_dPdx, dA_dPdy, and dA_dPdz to reflect the
 *   computed values.
 * 
 * Note:
 * - The function prints warnings to standard output if computed values are NaN or unexpectedly large,
 *   setting those values to zero to maintain stability.
 */
void NeuronGrowth::EvaluateOrientation(const int nen, const vector<float> &Nx, const vector<array<float, 3>> &dNdx, const vector<float> elePhi, const vector<float> eleTheta,  float& eleAniso, float& dA_dPdx, float& dA_dPdy, float& dA_dPdz)
{
	dA_dPdx = 0; dA_dPdy = 0; dA_dPdz = 0;
	for (int i = 0; i < nen; i++) {
		// if ((elePhi[i] > 0.05) && (elePhi[i] < 0.95)) {	
			float dPdx4 = pow(elePhi[i] * dNdx[i][0], 4);
			float dPdy4 = pow(elePhi[i] * dNdx[i][1], 4);
			float dPdz4 = pow(elePhi[i] * dNdx[i][2], 4);
			float dPdx3 = pow(elePhi[i] * dNdx[i][0], 3);
			float dPdy3 = pow(elePhi[i] * dNdx[i][1], 3);
			float dPdz3 = pow(elePhi[i] * dNdx[i][2], 3);
			float dPdx2 = pow(elePhi[i] * dNdx[i][0], 2);
			float dPdy2 = pow(elePhi[i] * dNdx[i][1], 2);
			float dPdz2 = pow(elePhi[i] * dNdx[i][2], 2);
			float dPdx1 = elePhi[i] * dNdx[i][0];
			float dPdy1 = elePhi[i] * dNdx[i][1];
			float dPdz1 = elePhi[i] * dNdx[i][2];

			float dP4 = dPdx4 + dPdy4 + dPdz4;
			float md4 = pow((dPdx2 + dPdy2 + dPdz2), 2);

			float stable(1e-2);

			float tmp(0);
			tmp = epsilonb * (1 - 3 * delta) + epsilonb * 4 * delta * dP4 / (md4 + stable);

			if ((isnan(tmp) == 1) || (abs(tmp) > 10)) {
				std::cout << "nan 1: " << tmp << std::endl;
				tmp = 0;
			}
			eleAniso += tmp;

			float C1x = dPdy4 + dPdz4;
			float C2x = dPdy2 + dPdz2;
			tmp = epsilonb * 4 * delta * ( (4 * dPdx3) / (pow(C2x + dPdx2, 2) + stable)
				- (4 * dPdx1 * (C1x + dPdx4)) / (pow((C2x + dPdx2), 3) + stable) );
			if ((isnan(tmp) == 1) || (abs(tmp) > 10)) {
				std::cout << "nan 2: " << tmp << std::endl;
				tmp = 0;
			}
			dA_dPdx += tmp;

			float C1y = dPdx4 + dPdz4;
			float C2y = dPdx2 + dPdz2;
			tmp =  epsilonb * 4 * delta * ( (4 * dPdy3) / (pow(C2y + dPdy2, 2) + stable)
				- (4 * dPdy1 * (C1y + dPdy4)) / (pow((C2y + dPdy2), 3) + stable) );
			if ((isnan(tmp) == 1) || (abs(tmp) > 10)) {
				std::cout << "nan 3: " << tmp << std::endl;
				tmp = 0;
			}
			dA_dPdy += tmp;
		
			float C1z = dPdx4 + dPdy4;
			float C2z = dPdx2 + dPdy2;
			tmp =  epsilonb * 4 * delta * ( (4 * dPdz3) / (pow(C2z + dPdz2, 2) + stable)
				- (4 * dPdz1 * (C1z + dPdz4)) / (pow((C2z + dPdz2), 3) + stable) );
			if ((isnan(tmp) == 1) || (abs(tmp) > 10)) {
				std::cout << "nan 4: " << tmp << std::endl;
				tmp = 0;
			}
			dA_dPdz += tmp;
		// }
	}
}

/**
 * @brief Evaluates the orientation of elements in a spherical coordinate system.
 *
 * This function computes the orientation of neuron growth elements based on their polar and azimuthal angles,
 * as well as their contribution to the overall orientation energy of the system, represented by `eleEpsilon`.
 * The gradients of the energy with respect to polar and azimuthal angles are also calculated to understand
 * the directional sensitivity of the growth process.
 *
 * @param nen The number of elements to evaluate.
 * @param Nx A vector of shape functions evaluated at the current point.
 * @param dNdx A vector of shape function derivatives with respect to the x, y, and z coordinates.
 * @param elePhi A vector containing the values of the field variable phi for each element.
 * @param elePolar A vector containing the polar angles for each element.
 * @param eleAzimuth A vector containing the azimuthal angles for each element.
 * @param eleEpsilon Reference to a float where the calculated orientation energy will be stored.
 * @param dEdp Reference to a float where the gradient of the energy with respect to the polar angle will be stored.
 * @param dEda Reference to a float where the gradient of the energy with respect to the azimuthal angle will be stored.
 *
 * @note The function modifies `eleEpsilon`, `dEdp`, and `dEda` in place to accumulate the orientation energy and its gradients.
 * @note The variable `delta` is assumed to be a predefined constant within the class or global scope, affecting the calculation.
 */
void NeuronGrowth::EvaluateOrientationSpherical(const int nen, const vector<float> &Nx, const vector<array<float, 3>> &dNdx, const vector<float> elePhi, const vector<float> elePolar, const vector<float> eleAzimuth, float& eleEpsilon, float dEdp, float dEda)
{
	for (int i = 0; i < nen; i++) {
		float polar = Nx[i] * elePolar[i];
		float azimuth = Nx[i] * eleAzimuth[i];
		float mode_grad_p = sqrt(pow(elePhi[i] * dNdx[i][0], 2) + pow(elePhi[i] * dNdx[i][1], 2) + pow(elePhi[i] * dNdx[i][2], 2));
		float p_angle = -acos((elePhi[i] * dNdx[i][2]) / mode_grad_p) - polar;
		float a_angle = -atan((elePhi[i] * dNdx[i][1]) / (elePhi[i] * dNdx[i][0])) - azimuth;
		float sin_p = sin(p_angle);
		float cos_p = cos(p_angle);
		float sin_a = sin(a_angle);
		float cos_a = cos(a_angle);

		eleEpsilon += 1 - 3 * delta + 4 * delta * (pow(sin_p, 4) * (pow(cos_a, 4) + pow(sin_a, 4)) + pow(cos_p, 4));
		dEdp += 4 * delta * ((pow(cos_a, 4) + pow(sin_a, 4)) * 4 * pow(sin_p, 3) * cos_p - 4 * pow(cos_p, 3) * sin_p);
		dEda += 4 * delta * (pow(sin_p, 4) * (-4 * pow(cos_a, 3) * sin_a + 4 * pow(sin_a, 3) * cos_a));
	}
}

/**
 * Calculate the sum of the gradient squared of Phi0 across the neuron growth model's Bézier mesh.
 *
 * This function iterates over each element in the neuron growth model's process-specific Bézier mesh,
 * computing the gradient of Phi0 at each element and accumulating the sum of the squares of these gradients.
 * The calculation involves evaluating the basis functions at Gaussian points to compute gradients
 * with respect to the mesh's local coordinates. This sum provides insight into the variation of Phi0
 * across the mesh, which is useful for understanding the growth dynamics in the neuron growth model.
 *
 * @param cpts A constant reference to a vector of Vertex3D objects representing control points in the mesh.
 *             This parameter is not directly used in the function but is part of the function signature for
 *             potential use in extended implementations or overrides.
 *
 * Process:
 * 1. Iterate over each element in the Bézier mesh.
 * 2. For each element, prepare a vector of Phi0 values corresponding to the element's control points.
 * 3. Iterate through a set of Gaussian points to calculate the determinant of the Jacobian (detJ) and gradients of Phi0.
 * 4. For each Gaussian point, calculate the gradients (dP0dx, dP0dy, dP0dz) using the basis function derivatives.
 * 5. Accumulate the sum of the squares of the gradients for all elements in the mesh.
 *
 * Note:
 * - The function uses Gaussian quadrature for integration over each element's domain to compute gradients.
 * - The BasisFunction and ElementDeriv functions are used to compute basis function values and gradients, respectively.
 * - The function updates the member variable `sum_grad_phi0_local` with the accumulated sum of gradient squares.
 */
void NeuronGrowth::CalculateSumGradPhi0(const vector<Vertex3D>& cpts) {
	// Iterate over each element in the process-specific Bézier mesh
	for (int e = 0; e < bzmesh_process.size(); e++) {
		int nen = bzmesh_process[e].IEN.size();
		vector<float> elePhi0(nen);

		// Prepare element-specific phi_0 values
		for (int i = 0; i < nen; i++) {
			elePhi0[i] = phi_0[bzmesh_process[e].IEN[i]];
		}

		// Iterate through Gaussian points to calculate gradients
		for (int i = 0; i < Gpt.size(); i++) {
			for (int j = 0; j < Gpt.size(); j++) {
				for (int k = 0; k < Gpt.size(); k++) {
					float detJ;
					float dudx[3][3];
					vector<float> Nx;
					vector<array<float, 3>> dNdx;

					BasisFunction(Gpt[i], Gpt[j], Gpt[k], bzmesh_process[e].pts, bzmesh_process[e].cmat, Nx, dNdx, dudx, detJ);
					detJ = wght[i] * wght[j] * wght[k] * detJ;

					float dP0dx, dP0dy, dP0dz;
					ElementDeriv(nen, dNdx, elePhi0, dP0dx, dP0dy, dP0dz);

					// Accumulate local gradient sum
					sum_grad_phi0_local += pow(dP0dx, 2) + pow(dP0dy, 2) + pow(dP0dz, 2);
				}
			}
		}
	}
}

/**
 * Builds the linear system for synaptic and tubulin dynamics within each process of the neuron growth model.
 *
 * This function constructs and assembles the linear systems for synaptic transmission (syn) and tubulin (tub)
 * concentration dynamics across the neuron mesh. It utilizes a detailed computational model to simulate
 * the complex interactions and growth processes of neurons, focusing on two key aspects: the synaptic
 * strength (syn) and tubulin concentration (tub), which are crucial for neural connectivity and structure formation.
 *
 * @param cpts The control points (vertices) of the neuron mesh, encapsulating the spatial structure of the neuron.
 *
 * Key steps involved:
 * 1. Iterates over each element (bezier mesh patch) in the process-specific neuron mesh.
 * 2. Initializes and clears matrices and vectors for solving syn and tub dynamics.
 * 3. Resizes solution matrices and vectors based on the number of element nodes (nen).
 * 4. Sets up elemental values (eleVal_st) for cell boundary conditions, synaptic strength, and tubulin concentration.
 * 5. Iterates through Gaussian points (Gpt) to evaluate elemental contributions to the solution vectors and matrices.
 * 6. Applies boundary conditions for syn and tub dynamics based on the neuron domain boundaries and cell properties.
 * 7. Assembles the global residual vectors and stiffness matrices for syn and tub dynamics.
 * 8. Finalizes the assembly of global vectors and matrices for incorporation into the broader neuron growth model.
 *
 * The function leverages adaptive computational techniques to efficiently simulate the distribution and evolution
 * of synaptogenesis and tubulin concentrations, contributing to the overall understanding of neuron development
 * and neurodevelopmental disorders.
 *
 * Note: The function assumes access to global variables and methods for matrix and vector operations, boundary condition
 * application, and element evaluation specific to synaptic and tubulin dynamics.
 */
void NeuronGrowth::BuildLinearSystemProcessNG_syn_tub(const vector<Vertex3D> &cpts)
{
	/*Build linear system in each process*/
	int ind(0);
	for (int e = 0; e < bzmesh_process.size(); e++) {
		int nen = bzmesh_process[e].IEN.size();

		vector<vector<float>> EMatrixSolve_syn;
		vector<float> EVectorSolve_syn;
		EMatrixSolve_syn.clear();
		EVectorSolve_syn.clear();
		EMatrixSolve_syn.resize(nen * 1);
		EVectorSolve_syn.resize(nen * 1);
		
		vector<vector<float>> EMatrixSolve_tub;
		vector<float> EVectorSolve_tub;
		EMatrixSolve_tub.clear();
		EVectorSolve_tub.clear();
		EMatrixSolve_tub.resize(nen);
		EVectorSolve_tub.resize(nen);

		vector<vector<float>> eleVal_st;
		eleVal_st.clear();
		eleVal_st.resize(5);

		for (int i = 0; i < nen * 1; i++) {
			EMatrixSolve_syn[i].resize(nen * 1, 0.0);
			EMatrixSolve_tub[i].resize(nen, 0.0);
			for (int j = 0; j < nen * 1; j++)
			{
				EMatrixSolve_syn[i][j] = 0.0;
				EMatrixSolve_tub[i][j] = 0.0;
			}
			EVectorSolve_syn[i] = 0.0;
			EVectorSolve_tub[i] = 0.0;

			if (i < 5) 
				eleVal_st[i].resize(nen);
		}
		
		for (int i = 0; i < nen * 1; i++) {
			eleVal_st[0][i] = phi[bzmesh_process[e].IEN[i]] - phi_prev[bzmesh_process[e].IEN[i]];		// elePhiDiff
			eleVal_st[1][i] = syn[bzmesh_process[e].IEN[i]];						// eleSyn
			eleVal_st[2][i] = CellBoundary(phi[bzmesh_process[e].IEN[i]], 0.5);				// elePhi
			eleVal_st[3][i] = CellBoundary(phi_prev[bzmesh_process[e].IEN[i]], 0.5);			// elePhiPrev
			eleVal_st[4][i] = tub[bzmesh_process[e].IEN[i]];						// eleConct
		}

		for (int i = 0; i < Gpt.size(); i++) {
			for (int j = 0; j < Gpt.size(); j++) {
				for (int k = 0; k < Gpt.size(); k++) {

					vector<float> vars_st;
					vars_st.clear();
					vars_st.resize(13);
					//        0      1     2      3	  4      5        6          7          8          9         10          11	     12
					// float elePf, eleS, eleP, dPdx, dPdy, elePprev, eleC, mag_grad_phi0, term_diff, term_alph, term_beta, term_source, dPdz;

					vars_st[7] = pre_mag_grad_phi0[ind];
					vars_st[11] = pre_term_source[ind];

					ElementEvaluationAll_syn_tub(nen, pre_Nx[ind], pre_dNdx[ind], eleVal_st, vars_st);

					for (int m = 0; m < nen; m++) {
						EVectorSolve_syn[m] += (vars_st[1] + kappa * vars_st[0]) * pre_Nx[ind][m] * pre_detJ[ind];
		
						EVectorSolve_tub[m] += (vars_st[6] * vars_st[2] + dt/2 * vars_st[11]) * pre_Nx[ind][m] * pre_detJ[ind];

						for (int n = 0; n < nen; n++) {
							if (judge_syn == 0)
								EMatrixSolve_syn[m][n] += (pre_Nx[ind][m] * pre_Nx[ind][n] + dt/2 * Dc * (pre_dNdx[ind][m][0] * pre_dNdx[ind][n][0] + pre_dNdx[ind][m][1] * pre_dNdx[ind][n][1] + pre_dNdx[ind][m][2] * pre_dNdx[ind][n][2])) * pre_detJ[ind];

							EMatrixSolve_tub[m][n] += (pre_Nx[ind][m] * vars_st[2] * pre_Nx[ind][n]
								- dt/2 * (- Diff * (vars_st[2] * (pre_dNdx[ind][m][0] * pre_dNdx[ind][n][0] + pre_dNdx[ind][m][1] * pre_dNdx[ind][n][1] + pre_dNdx[ind][m][2] * pre_dNdx[ind][n][2]))
									- alphaT * (vars_st[2] * (pre_dNdx[ind][m][0] + pre_dNdx[ind][m][1] + pre_dNdx[ind][m][2]) + (vars_st[3] + vars_st[4] + vars_st[12]) * pre_Nx[ind][m]) * pre_Nx[ind][n]
									- betaT * (vars_st[2] * pre_Nx[ind][m]) * pre_Nx[ind][n] ) 
									) * pre_detJ[ind];
						}
					}
					ind += 1;
				}
			}
		}

		/*Apply Boundary Condition*/
		for (int i = 0; i < nen; i++) {
			int A = bzmesh_process[e].IEN[i];
			if (cpts[A].label == 1) // domain boundary
				ApplyBoundaryCondition(0, i, 0, EMatrixSolve_syn, EVectorSolve_syn);
			if (CellBoundary(phi[A], 0.5) != 1) // outside soma bc for tubulin
				ApplyBoundaryCondition(0, i, 0, EMatrixSolve_tub, EVectorSolve_tub);
		}

		ResidualAssembly(EVectorSolve_syn, bzmesh_process[e].IEN, GR_syn);
		if (judge_syn == 0) // The matrix is the same, so only need to assembly once
			MatrixAssembly(EMatrixSolve_syn, bzmesh_process[e].IEN, GK_syn);

		ResidualAssembly(EVectorSolve_tub, bzmesh_process[e].IEN, GR_tub);
		MatrixAssembly(EMatrixSolve_tub, bzmesh_process[e].IEN, GK_tub);
	}

	VecAssemblyBegin(GR_syn);
	if (judge_syn == 0)
		MatAssemblyBegin(GK_syn, MAT_FINAL_ASSEMBLY);

	VecAssemblyBegin(GR_tub);
	MatAssemblyBegin(GK_tub, MAT_FINAL_ASSEMBLY);
}

int NeuronGrowth::CheckExpansion3D(vector<float> input, const std::vector<Vertex3D>& cpts, int NX, int NY, int NZ, int originX, int originY, int originZ) 
{
	float bc_clearance = 5;
	for (int i = 0; i < cpts.size(); i++) {
		if (CellBoundary(input[i],0) > 0) {
			float currX = cpts[i].coor[0] - originX;
			float currY = cpts[i].coor[1] - originY;
			float currZ = cpts[i].coor[2] - originZ;

			// 0 - left | 1 - top | 2 - right | 3 - bottom | 4 - front | 5 - back | 6 - no action
			if (currX >= (max_x - bc_clearance)) {
				return 3;
			} else if (currX <= (bc_clearance - originX)) {
				return 1;
			}

			if (currY >= (max_y - bc_clearance)) {
				return 0;
			} else if (currY <= (bc_clearance - originY)) {
				return 2;
			}

			if (currZ >= (max_z - bc_clearance)) {
				return 4;
			} else if (currZ <= (bc_clearance - originZ)) {
				return 5;
			}
		}
	}

	return 6;
}

// void NeuronGrowth::ExpandDomain(vector<float> input, vector<float> &expd_var, int edge) 
// {
// 	int length = input.size();
// 	int sz = sqrt(length);

// 	int expd_sz = 10; // directional expanding size
// 	int new_sz = sz + expd_sz;

// 	expd_var.clear(); expd_var.resize(pow(new_sz,2));
// 	for (int i = 0; i < new_sz; i++)
// 		expd_var[i] = 0;
	
// 	int ind;
// 	switch (edge) {
// 		case 0: // left
// 			ind = 0;
// 			for (int i = expd_sz; i < new_sz-1; i++) {
// 				for (int j = (int)(expd_sz/2); j < (new_sz - (int)(expd_sz/2)); j++) {
// 					expd_var[i * new_sz + j] = input[ind];
// 					ind += 1;
// 				}
// 			}
// 			// std::cout << "Expanding left!" << std::endl;
// 			break;
// 		case 1: // top
// 			ind = 0;
// 			for (int i = (int)(expd_sz/2); i < (new_sz - (int)(expd_sz/2)); i++) {
// 				for (int j = 0; j < new_sz-expd_sz; j++) {
// 					expd_var[i * new_sz + j] = input[ind];
// 					ind += 1;
// 				}
// 			}
// 			// std::cout << "Expanding top!" << std::endl;
// 			break;
// 		case 2: // right
// 			ind = 0;
// 			for (int i = 0; i < new_sz-1-expd_sz; i++) {
// 				for (int j = (int)(expd_sz/2); j < (new_sz - (int)(expd_sz/2)); j++) {
// 					expd_var[i * new_sz + j] = input[ind];
// 					ind += 1;
// 				}
// 			}
// 			// std::cout << "Expanding right!" << std::endl;
// 			break;
// 		case 3: // bottom
// 			ind = 0;
// 			for (int i = (int)(expd_sz/2); i < (new_sz - (int)(expd_sz/2)); i++) {
// 				for (int j = expd_sz; j < new_sz; j++) {
// 					expd_var[i * new_sz + j] = input[ind];
// 					ind += 1;
// 				}
// 			}
// 			// std::cout << "Expanding bottom!" << std::endl;
// 			break;
// 		case 4: // all direction
// 			ind = 0;
// 			for (int i = (int)(expd_sz/2); i < (new_sz - (int)(expd_sz/2)); i++) {
// 				for (int j = (int)(expd_sz/2); j < (new_sz - (int)(expd_sz/2)); j++) {
// 					expd_var[i * new_sz + j] = input[ind];
// 					ind += 1;
// 				}
// 			}
// 			// std::cout << "Expanding all direction!" << std::endl;
// 			break;
// 	}
// 	// input.swap(expd_var);
// }

// void NeuronGrowth::ExpandDomain(vector<float> input, vector<float> &expd_var, int edge, int NX, int NY) 
// {
// 	int expd_sz = 10; // directional expanding size
	
// 	expd_var.clear(); expd_var.resize((NX+1) * (NY+1));
// 	for (int i = 0; i < (expd_var.size()); i++)
// 		expd_var[i] = 0;
	
// 	// (0-left|1-top|2-right|3-bottom)
// 	int ind;
// 	switch (edge) {
// 		case 0: // left
// 			ind = 0;
// 			for (int i = expd_sz; i < NX; i++) {
// 				for (int j = 0; j <= NY; j++) {
// 					expd_var[i * (NY+1) + j] = input[ind];
// 					ind += 1;
// 				}
// 			}
// 			// std::cout << "Expanding left!" << std::endl;
// 			break;
// 		case 1: // top
// 			ind = 0;
// 			for (int i = 0; i < NX; i++) {
// 				for (int j = 0; j <= NY-expd_sz; j++) {
// 					expd_var[i * (NY+1) + j] = input[ind];
// 					ind += 1;
// 				}
// 			}
// 			// std::cout << "Expanding top!" << std::endl;
// 			break;
// 		case 2: // right
// 			ind = 0;
// 			for (int i = 0; i < NX-expd_sz; i++) {
// 				for (int j = 0; j <= NY; j++) {
// 					expd_var[i * (NY+1) + j] = input[ind];
// 					ind += 1;
// 				}
// 			}
// 			// std::cout << "Expanding right!" << std::endl;
// 			break;
// 		case 3: // bottom
// 			ind = 0;
// 			for (int i = 0; i < NX; i++) {
// 				for (int j = expd_sz; j <= NY; j++) {
// 					expd_var[i * (NY+1) + j] = input[ind];
// 					ind += 1;
// 				}
// 			}
// 			// std::cout << "Expanding bottom!" << std::endl;
// 			break;
// 	}
// }

// Populates vector elements with random values if they are zero
void NeuronGrowth::PopulateRandom(std::vector<float> &input) {
	std::for_each(input.begin(), input.end(), [](float &value) {
		if (value == 0) {
			value = static_cast<float>(rand() % 100) / 100.0f; // Random float between 0 and 1
		}
	});
}

// Removes outliers from the data vector based on standard deviation and returns the adjusted threshold
float NeuronGrowth::RmOutlier(std::vector<float> &data) {
	float sum = std::accumulate(data.begin(), data.end(), 0.0f);
	float mean = sum / data.size();

	float sq_sum = std::inner_product(data.begin(), data.end(), data.begin(), 0.0,
										[](float acc, float val) { return acc + val; },
										[mean](float a, float b) { return std::pow(a - mean, 2) + b; });
	float standardDeviation = std::sqrt(sq_sum / data.size());

	float threshold = mean + 3 * standardDeviation;
	// Clamp values exceeding threshold to threshold
	std::transform(data.begin(), data.end(), data.begin(), [threshold](float value) {
		return std::min(value, threshold);
	});

	return mean + 2 * standardDeviation;
}

// Determines if a cell's value exceeds a given threshold and returns 1.0f if true, else 0.0f
float NeuronGrowth::CellBoundary(float phi, float threshold) {
	return phi > threshold ? 1.0f : 0.0f;
}

// // Function to apply a simple smoothing operation to a 2D binary variable
// std::vector<float> NeuronGrowth::SmoothBinary2D(const std::vector<float>& binaryData, int rows, int cols) {
// 	// Create a new vector to store the smoothed data
// 	std::vector<float> smoothedData(rows * cols, 0.0f);

// 	// Define a 5x5 smoothing kernel
// 	const std::vector<float> kernel = {
// 		1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
// 		1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
// 		1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
// 		1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
// 	1.0f, 1.0f, 1.0f, 1.0f, 1.0f
// 	};

// 	// Apply the smoothing operation using the defined kernel
// 	for (int i = 2; i < rows - 2; ++i) {
// 		for (int j = 2; j < cols - 2; ++j) {
// 			float sum = 0.0f;
// 			for (int k = -2; k <= 2; ++k) {
// 				for (int l = -2; l <= 2; ++l) {
// 					sum += binaryData[(i + k) * cols + (j + l)] * kernel[(k + 2) * 5 + (l + 2)];
// 				}
// 			}
// 			smoothedData[i * cols + j] = (sum > 12.0f) ? 1.0f : 0.0f; // Apply threshold to convert sum to binary value
// 		}
// 	}

// 	return smoothedData;
// }

// void NeuronGrowth::DetectTipsMulti3D(vector<float> id, int numNeuron, vector<float> &tips, int NX, int NY, int NZ)
// {
// 	float threshold(0.90), maxVal(0);
// 	int ind, length((NX + 1) * (NY + 1) * (NZ + 1));
// 	tips.clear();
// 	tips.resize(length);

// 	for (int i = (5 * (NX + 1) * (NY + 1) + 5); i < (length - 4 * (NX + 1) * (NY + 1) - 4); i++) {
// 		tips[i] = 0;

// 		// Check if the cell is near the boundary
// 		if (CellBoundary(phi[i], 0.25) > 0) {
// 			for (int j = -4; j < 5; j++) {
// 				for (int k = -4; k < 5; k++) {
// 					for (int l = -4; l < 5; l++) {
// 						for (int m = 0; m < numNeuron; m++) {
// 							// Accumulate values based on neighboring cells and neuron IDs
// 							if ((m + 1) == id[i + j * (NY + 1) + k * (NY + 1) * (NZ + 1) + l * (NZ + 1)]) {
// 								tips[i] += CellBoundary(phi[i + j * (NY + 1) + k * (NY + 1) * (NZ + 1) + l * (NZ + 1)], 0.1);
// 							}
// 						}
// 					}
// 				}
// 			}

// 			// Normalize tips value
// 			if (tips[i] > 0)
// 				tips[i] = CellBoundary(phi[i], 0) / tips[i];

// 			// Handle NaN cases
// 			if (std::isnan(tips[i]))
// 				tips[i] = 0;

// 			// Update maxVal
// 			if (tips[i] > maxVal)
// 				maxVal = tips[i];
// 		}
// 	}

// 	// Thresholding and setting tips values
// 	for (int i = 1 + (NX + 1) * (NY + 1); i < length - (NX + 1) * (NY + 1) - 1; i++) {
// 		if (tips[i] < (threshold * maxVal)) {
// 			tips[i] = 0;
// 		} else {
// 			tips[i] = 1;
// 		}
// 	}
// }

// Check if point coordinates are within the bounds defined by center and dx dy dz
bool NeuronGrowth::isInBox(const Vertex3D& point, const Vertex3D& center, float dx, float dy, float dz) {
	return point.coor[0] >= (center.coor[0] - dx/2) && point.coor[0] <= (center.coor[0] + dx/2) &&
		point.coor[1] >= (center.coor[1] - dy/2) && point.coor[1] <= (center.coor[1] + dy/2) &&
		point.coor[2] >= (center.coor[2] - dz/2) && point.coor[2] <= (center.coor[2] + dz/2);
}

/**
 * Calculates the sum of phi values around each vertex in a given set, applying
 * a threshold to identify significant points, potentially indicating neuron growth tips.
 * 
 * This function iterates over a collection of vertices (`cpts`), summing the phi values
 * within a defined vicinity around each vertex. The vicinity is determined by the `dx`, `dy`, and `dz`
 * parameters, which define the dimensions of a box centered on each vertex. Points within
 * this box contribute to the sum. After calculating the sums, the function applies a threshold
 * to these values, normalizing them against the highest value found. Points with summed values
 * above this threshold are marked as potential growth tips by setting their corresponding value
 * in the output vector to 1; all others are set to 0.
 * 
 * @param cpts A vector of Vertex3D objects representing the neuron's vertices.
 * @param dx The delta in the x-direction to define the vicinity around a point.
 * @param dy The delta in the y-direction to define the vicinity around a point.
 * @param dz The delta in the z-direction to define the vicinity around a point.
 * @return std::vector<float> A vector of the same size as `cpts`, where each element is either 0
 *         (indicating the corresponding vertex is not a tip) or 1 (indicating a potential tip),
 *         based on the thresholding of summed phi values.
 * 
 * Note: The function assumes that `phi` is a pre-defined vector accessible within the class
 *       that contains phi values corresponding to each Vertex3D in `cpts`. The function
 *       `isInBox` checks whether a point is within the specified vicinity of another,
 *       and `CellBoundary` computes a boundary-related value for a given phi, with
 *       the second argument presumably allowing for further customization.
 */
void NeuronGrowth::calculatePhiSum(const std::vector<Vertex3D>& cpts, float dx, float dy, float dz) {
	// std::vector<float> tips(cpts.size(), 0); // Initialize the result vector with zeros
	tips.clear(); tips.resize(cpts.size());
	float threshold = 0.95; // Threshold for filtering tips
	float maxVal = 0; // Track the maximum value of tips for normalization

	// Iterate over each center point
	for (int i = 0; i < cpts.size(); ++i) {
		const auto& center = cpts[i];

		tips[i] = 0;

		// Sum phi values for points within the box centered at 'center'
		for (int j = 0; j < phi.size(); ++j) {
			if (isInBox(cpts[j], center, dx, dy, dz)) {
				tips[i] += CellBoundary(phi[j], 0.5) * distI[j];
				// tips[i] += CellBoundary(phi[j], 0.5);
			}
		}

		tips[i] = CellBoundary(phi[i], 0.5) / tips[i] * CellBoundary(phi[i], 0.5);
		// Handle NaN cases, ensuring a valid tips value
		if (std::isnan(tips[i])) {
			tips[i] = 0;
		}
		
		if ((center.coor[0] <= min_x+1) || (center.coor[0] >= max_x-1) 
			|| (center.coor[1] <= min_y+1) || (center.coor[1] >= max_y-1) 
			|| (center.coor[2] <= min_z+1) || (center.coor[2] >= max_z-1)) {
			tips[i] = 0;
		}

		// Update maxVal to the highest tips value found
		if (CellBoundary(phi[i], 0.5) > 0)
		maxVal = max(tips[i], maxVal);
	}

	// for (int i = 0; i < tips.size(); ++i) {
	// 	if (tips[i] > (0.99 * maxVal)) {
	// 		tips[i] = 0; // Set tips below threshold to 0
	// 	}
	// }

	// maxVal = 0;
	// for (int i = 0; i < tips.size(); ++i) {
	// 	maxVal = max(tips[i], maxVal);
	// }

	CheckVar("../io3D/outputs/TIP_", cpts, tips);
	CheckVar("../io3D/outputs/PHI_", cpts, phi);

	// Thresholding and setting tips values based on the maximum value found
	for (int i = 0; i < tips.size(); ++i) {
		if (tips[i] > (threshold * maxVal)) {
		// if (tips[i] > (0.018)) {
			tips[i] = 1; // Set tips below threshold to 0
		} else {
			tips[i] = 0; // Set tips above threshold to 1
		}
		
		// if ((cpts[i].coor[0] >= 9) && (cpts[i].coor[0] <= 11)
		// 	&& (cpts[i].coor[1] >= 9) && (cpts[i].coor[1] <= 11)) {
		// 	tips[i] = 1;
		// } else {
		// 	tips[i] = 0;
		// }

	}

	// // Thresholding and setting tips values based on the maximum value found
	// for (int i = 0; i < tips.size(); ++i) {
	// 	if (tips[i] == 1) {
	// 		for (int j = 0; j < tips.size(); ++i) {
	// 			if (isInBox(cpts[j], center, dx, dy, dz)) {
	// 				tips[i] += CellBoundary(phi[j], 0.5);
	// 			}
	// 		}
	// 	}
	// }

	// return tips; // Return the processed tips vector
}


// Function to interpolate values for a new mesh based on coordinates
vector<float> NeuronGrowth::InterpolateValues3D(const vector<Vertex3D>& cpts_initial, const vector<float>& input,
                                      const vector<Vertex3D>& cpts_new) {
	vector<float> output;
	output.resize(cpts_new.size());

	for (int i = 0; i < cpts_new.size(); i++) {
		// std::cout << i << std::endl;
		float x = cpts_new[i].coor[0];
		float y = cpts_new[i].coor[1];
		float z = cpts_new[i].coor[2];

		int ind;
		if (SearchVertex(cpts_initial, x, y, z, ind)) {
			// Exact match found, no need for interpolation
			output[i] = input[ind];
			// std::cout << ind << " ";
		} 
		else {
			// int numNeighbor = 0;
			// while (numNeighbor < 6) {
			// 	for (int j = 0; j < cpts.size(); ++i) {
			// 		if (isInBox(cpts[j], center, dx, dy, dz)) {
			// 			numNeighbor += 1;
			// 		}
			// 	}
			// }

			// Linear interpolation for the new coordinates
			// Find the vertices around the target coordinates
			float x1 = floor(x), x2 = ceil(x);
			float y1 = floor(y), y2 = ceil(y);
			float z1 = floor(z), z2 = ceil(z);

			// Find the corresponding indices in the initial mesh
			int ind111, ind112, ind121, ind122, ind211, ind212, ind221, ind222;
			SearchVertex(cpts_initial, x1, y1, z1, ind111);
			SearchVertex(cpts_initial, x1, y1, z2, ind112);
			SearchVertex(cpts_initial, x1, y2, z1, ind121);
			SearchVertex(cpts_initial, x1, y2, z2, ind122);
			SearchVertex(cpts_initial, x2, y1, z1, ind211);
			SearchVertex(cpts_initial, x2, y1, z2, ind212);
			SearchVertex(cpts_initial, x2, y2, z1, ind221);
			SearchVertex(cpts_initial, x2, y2, z2, ind222);

			// Interpolate along each dimension separately
			float interpX1 = Lerp(input[ind111], input[ind112], (x - x1));
			float interpX2 = Lerp(input[ind121], input[ind122], (x - x1));
			float interpY1 = Lerp(interpX1, interpX2, (y - y1));

			float interpX3 = Lerp(input[ind211], input[ind212], (x - x2));
			float interpX4 = Lerp(input[ind221], input[ind222], (x - x2));
			float interpY2 = Lerp(interpX3, interpX4, (y - y1));

			// Interpolate along the z dimension
			output[i] = Lerp(interpY1, interpY2, (z - z1));

			// std::cout << ind << " ";
		}
	}
	return output;
}

// Function to calculate the squared distance between two vertices
// Using squared distance to avoid the computational cost of the square root
float SquaredDistance(const Vertex3D& first, const Vertex3D& other) {
	return (first.coor[0] - other.coor[0]) * (first.coor[0] - other.coor[0]) + 
		(first.coor[1] - other.coor[1]) * (first.coor[1] - other.coor[1]) + 
		(first.coor[2] - other.coor[2]) * (first.coor[2] - other.coor[2]);
}

/**
 * Finds and returns the k closest vertices to a specified input vertex from a list of vertices, including their indices.
 * This function is useful for identifying nearby vertices in a mesh or point cloud, which can be essential for 
 * applications requiring spatial analysis or nearest neighbor searches, such as mesh refinement, topology optimization,
 * or even in simulations of physical phenomena like neuron growth.
 *
 * The function uses a priority queue to efficiently manage and sort distances between the input vertex and each vertex
 * in the list. It ensures that only the k closest vertices (based on squared Euclidean distance) are retained and
 * returned in ascending order of distance.
 *
 * @param vertices A vector of Vertex3D objects representing the list of vertices to search through.
 * @param inputVertex The Vertex3D object for which the k closest vertices are to be found.
 * @param k The number of closest vertices to find. If k is larger than the number of available vertices, all vertices
 *          are returned sorted by proximity to the inputVertex.
 *
 * @return A vector of pairs, each containing a Vertex3D object (one of the k closest vertices) and its corresponding
 *         index in the original 'vertices' vector. The pairs are sorted by increasing distance from the inputVertex.
 *
 * @note The function assumes 'Vertex3D' is a class or struct that represents a 3D point with accessible coordinates.
 *       The 'SquaredDistance' function used for comparison must be defined elsewhere, calculating the squared
 *       Euclidean distance between two Vertex3D objects.
 */
std::vector<std::pair<Vertex3D, int>> NeuronGrowth::FindClosestVerticesWithIndices(const std::vector<Vertex3D>& vertices, const Vertex3D& inputVertex, int k) {
	// Custom comparator that prioritizes larger squared distances and considers the vertex index
	auto comp = [&inputVertex](const std::pair<Vertex3D, int>& a, const std::pair<Vertex3D, int>& b) {
		return SquaredDistance(inputVertex, a.first) < SquaredDistance(inputVertex, b.first);
	};

	// Initialize the priority queue with the custom comparator
	std::priority_queue<std::pair<Vertex3D, int>, std::vector<std::pair<Vertex3D, int>>, decltype(comp)> pq(comp);

	for (int i = 0; i < vertices.size(); ++i) {
		pq.emplace(vertices[i], i);
		if (pq.size() > k) pq.pop(); // Remove the vertex with the largest distance
	}

	std::vector<std::pair<Vertex3D, int>> closestVertices;
	while (!pq.empty()) {
		closestVertices.push_back(pq.top()); // Collect the closest vertices and their indices
		pq.pop();
	}
	// Reverse to correct the order from farthest of the closest to nearest
	std::reverse(closestVertices.begin(), closestVertices.end());

	return closestVertices;
}

std::vector<std::tuple<Vertex3D, int, float>> NeuronGrowth::FindClosestVerticesWithIndicesAndDistances(const std::vector<Vertex3D>& vertices, const Vertex3D& inputVertex, int k) {
	// Custom comparator that prioritizes larger squared distances and considers the vertex index
	auto comp = [&inputVertex](const std::tuple<Vertex3D, int, float>& a, const std::tuple<Vertex3D, int, float>& b) {
		return std::get<2>(a) < std::get<2>(b); // Compare based on squared distance
	};

	// Initialize the priority queue with the custom comparator
	std::priority_queue<std::tuple<Vertex3D, int, float>, std::vector<std::tuple<Vertex3D, int, float>>, decltype(comp)> pq(comp);

	for (int i = 0; i < vertices.size(); ++i) {
		float distance = SquaredDistance(inputVertex, vertices[i]);
		pq.emplace(vertices[i], i, distance);
		if (pq.size() > k) pq.pop(); // Remove the vertex with the largest distance
	}

	std::vector<std::tuple<Vertex3D, int, float>> closestVertices;
	while (!pq.empty()) {
		closestVertices.push_back(pq.top()); // Collect the closest vertices, their indices, and squared distances
		pq.pop();
	}
	// Reverse to correct the order from farthest of the closest to nearest
	std::reverse(closestVertices.begin(), closestVertices.end());

	return closestVertices;
}


// Function to perform Breadth-First Search (BFS) for clustering in 3D
void NeuronGrowth::bfs3D(const vector<float>& matrix, int depth, int rows, int cols, int dep, int row, int col,
                          vector<bool>& visited, vector<tuple<int, int, int>>& cluster) {
	static const int dx[] = {-1, 0, 1, 0, 0, 0};  // Directions in x
	static const int dy[] = {0, -1, 0, 1, 0, 0};  // Directions in y
	static const int dz[] = {0, 0, 0, 0, -1, 1};  // Directions in z

	queue<tuple<int, int, int>> q;
	q.push(make_tuple(dep, row, col));
	visited[dep * rows * cols + row * cols + col] = true;
	cluster.push_back(make_tuple(dep, row, col));

	while (!q.empty()) {
		int d = get<0>(q.front());
		int r = get<1>(q.front());
		int c = get<2>(q.front());
		q.pop();

		for (int i = 0; i < 6; ++i) {  // Iterating through all 6 directions
			int newDep = d + dz[i];
			int newRow = r + dx[i];
			int newCol = c + dy[i];

			if (newDep >= 0 && newDep < depth && newRow >= 0 && newRow < rows && newCol >= 0 && newCol < cols &&
				matrix[newDep * rows * cols + newRow * cols + newCol] != 0 && !visited[newDep * rows * cols + newRow * cols + newCol]) {
				visited[newDep * rows * cols + newRow * cols + newCol] = true;
				q.push(make_tuple(newDep, newRow, newCol));
				cluster.push_back(make_tuple(newDep, newRow, newCol));
			}
		}
	}
}

// Function to find connected clusters in the 3D matrix
vector<vector<tuple<int, int, int>>> NeuronGrowth::FindClusters3D(const vector<float>& matrix, int depth, int rows, int cols) {
	vector<bool> visited(depth * rows * cols, false);
	vector<vector<tuple<int, int, int>>> clusters;

	for (int d = 0; d < depth; ++d) {
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				if (matrix[d * rows * cols + i * cols + j] != 0 && !visited[d * rows * cols + i * cols + j]) {
					vector<tuple<int, int, int>> cluster;
					bfs3D(matrix, depth, rows, cols, d, i, j, visited, cluster);

					if (!cluster.empty()) {
						clusters.push_back(cluster);
					}
				}
			}
		}
	}
	return clusters;
}

// Function to find local maxima within connected clusters in the 3D matrix
vector<float> NeuronGrowth::FindLocalMaximaInClusters3D(const vector<float>& matrix, int depth, int rows, int cols) {
    vector<vector<tuple<int, int, int>>> clusters = FindClusters3D(matrix, depth, rows, cols);
    vector<float> localMaxima(depth * rows * cols, 0.0f);

	for (const auto& cluster : clusters) {
		float maxVal = -numeric_limits<float>::max();
		tuple<int, int, int> maxPos = make_tuple(-1, -1, -1);

		for (const auto& pos : cluster) {
			int dep = get<0>(pos);
			int row = get<1>(pos);
			int col = get<2>(pos);

			if (matrix[dep * rows * cols + row * cols + col] > maxVal) {
				maxVal = matrix[dep * rows * cols + row * cols + col];
				maxPos = pos;
			}
		}

		if (get<0>(maxPos) != -1 && get<1>(maxPos) != -1 && get<2>(maxPos) != -1) {
			localMaxima[get<0>(maxPos) * rows * cols + get<1>(maxPos) * cols + get<2>(maxPos)] = 1.0f;

			// Include neighbors
			for (int i = -1; i <= 1; ++i) {
				for (int j = -1; j <= 1; ++j) {
					for (int k = -1; k <= 1; ++k) {
						int newDep = get<0>(maxPos) + i;
						int newRow = get<1>(maxPos) + j;
						int newCol = get<2>(maxPos) + k;

						if (newDep >= 0 && newDep < depth && newRow >= 0 && newRow < rows && newCol >= 0 && newCol < cols) {
							localMaxima[newDep * rows * cols + newRow * cols + newCol] = 1.0f;
						}
					}
				}
			}
		}
	}
	return localMaxima;
}

vector<vector<vector<int>>> NeuronGrowth::ConvertTo3DIntVector(const vector<float> input, int NX, int NY, int NZ) 
{
	vector<vector<vector<int>>> output;

	// std::cout << input.size() << std::endl;
	std::cout << input.size() << " " << NX << " " << NY << " " << NZ << std::endl;

	int k = 0;
	for (int i = 0; i < NX+1; i++) {
		vector<vector<int>> matrix;
		for (int j = 0; j < NY+1; j++) {
			vector<int> row;
			for (int z = 0; z < NZ+1; z++) {
				row.push_back(CellBoundary(input[k], 0));
				// std::cout << k;
				k++;
			}
			matrix.push_back(row);
		}
		output.push_back(matrix);
	}
	return output;
}

vector<vector<vector<float>>> NeuronGrowth::ConvertTo3DFloatVector(const vector<float> input, int NX, int NY, int NZ) 
{
	vector<vector<vector<float>>> output;

	int k = 0;
	for (int i = 0; i < NX; i++) {
		vector<vector<float>> matrix;
		for (int j = 0; j < NY; j++) {
			vector<float> row;
			for (int z = 0; z < NZ; z++) {
				row.push_back(input[k]);
				k++;
			}
			matrix.push_back(row);
		}
		output.push_back(matrix);
	}
	return output;
}

void NeuronGrowth::FloodFill3D(std::vector<std::vector<std::vector<int>>>& image, int x, int y, int z, int newColor, int originalColor) 
{
	// Define the directions: up, down, left, right, front, back
	int dx[] = {0, 0, -1, 1, 0, 0};
	int dy[] = {-1, 1, 0, 0, 0, 0};
	int dz[] = {0, 0, 0, 0, -1, 1};

	if (x < 0 || x >= image.size() || y < 0 || y >= image[0].size() || z < 0 || z >= image[0][0].size() || image[x][y][z] != originalColor || image[x][y][z] == newColor) {
		return;
	}

	image[x][y][z] = newColor;

	// Apply flood fill in all six directions
	for (int i = 0; i < 6; ++i) {
		FloodFill3D(image, x + dx[i], y + dy[i], z + dz[i], newColor, originalColor);
	}
}

void NeuronGrowth::IdentifyNeurons3D(std::vector<std::vector<std::vector<int>>>& neurons, std::vector<std::array<int, 3>> seed, int NX, int NY, int NZ, int originX, int originY, int originZ) 
{
	std::cout << "ck0.1" << std::endl;
	neurons = ConvertTo3DIntVector(phi, NX, NY, NZ);
	std::cout << "ck0.2" << std::endl;

	for (int i = 0; i < seed.size(); i++) {
		std::cout << i << " " << std::endl;	
		int startX = seed[i][0] - originX;
		int startY = seed[i][1] - originY;
		int startZ = seed[i][2] - originZ;

		int newColor = i + 1;
		int originalColor = neurons[startX][startY][startZ];
		FloodFill3D(neurons, startX, startY, startZ, newColor, originalColor);
	}
	std::cout << "ck0.3" << std::endl;

}

bool NeuronGrowth::isValid(int x, int y, int z, int rows, int cols, int depth) 
{
    return (x >= 0 && x < rows && y >= 0 && y < cols && z >= 0 && z < depth);
}

std::vector<std::vector<std::vector<int>>> NeuronGrowth::CalculateGeodesicDistanceFromPoint3D(std::vector<std::vector<std::vector<int>>> neurons, const std::vector<std::array<int, 3>>& seed, int originX, int originY, int originZ) 
{
	int rows = neurons.size();
	int cols = (rows > 0) ? neurons[0].size() : 0;
	int depth = (cols > 0) ? neurons[0][0].size() : 0;

	std::vector<std::vector<std::vector<int>>> distances(rows, std::vector<std::vector<int>>(cols, std::vector<int>(depth, INF)));
	std::vector<std::vector<std::vector<bool>>> visited(rows, std::vector<std::vector<bool>>(cols, std::vector<bool>(depth, false)));

	for (const auto& point : seed) {
		int startX = point[0] - originX;
		int startY = point[1] - originY;
		int startZ = point[2] - originZ;

		distances[startX][startY][startZ] = 0;
		visited[startX][startY][startZ] = true;

		std::queue<std::array<int, 3>> q;
		q.push({startX, startY, startZ});

		int dx[] = {-1, 1, 0, 0, 0, 0};
		int dy[] = {0, 0, -1, 1, 0, 0};
		int dz[] = {0, 0, 0, 0, -1, 1};

		while (!q.empty()) {
			std::array<int, 3> current = q.front();
			q.pop();

			int x = current[0];
			int y = current[1];
			int z = current[2];

			for (int i = 0; i < 6; ++i) {
				int newX = x + dx[i];
				int newY = y + dy[i];
				int newZ = z + dz[i];

				if (isValid(newX, newY, newZ, rows, cols, depth) && neurons[newX][newY][newZ] == 1 && !visited[newX][newY][newZ]) {
					visited[newX][newY][newZ] = true;
					distances[newX][newY][newZ] = distances[x][y][z] + 1;
					q.push({newX, newY, newZ});
				}
			}
		}
	}
	return distances;
}

// // Function to calculate geodesic distance from a point
// vector<vector<array<int, 3>>> NeuronGrowth::NeuriteTracing(vector<vector<float>> distance) 
// {
// 	vector<vector<array<int, 2>>> traces;
// 	return traces;
// }

/**
 * Save neurodevelopmental growth variables (NGvars) to text files.
 *
 * This function saves each variable from a 2D vector of neurodevelopmental growth variables (NGvars)
 * into separate text files for visualization or further analysis. The names of the files are
 * derived from predefined variable names and an index, ensuring a structured and organized output.
 *
 * @param NGvars A 2D vector containing the neurodevelopmental growth variables to be saved. 
 *               Each inner vector represents a different variable (e.g., phi, syn, tub, etc.),
 *               and each element within an inner vector represents the variable's value at a specific point.
 * @param NX The number of points in the x-dimension. This parameter may be used for reshaping or
 *           interpreting the variables in a grid format, but is not directly used in this function's current implementation.
 * @param NY The number of points in the y-dimension. Similar to NX, it may be intended for grid interpretation
 *           but is not utilized within this function.
 * @param fn The base filename or path to which the variable text files will be saved. Each variable file
 *           will be named according to the variable it represents and will be stored in the directory specified by fn.
 *
 * Note:
 * - The function assumes that the directory specified by `fn` exists and is writable.
 * - The variable names are hard-coded within the function and correspond to the indices of the NGvars vector.
 * - Each file will contain data from its corresponding variable in NGvars, formatted for visualization if enabled.
 * - The visualization flag is currently set to true within the function, enabling visualization formatting by default.
 * - The actual saving of data to files is delegated to the `PrintVec2TXT` function, which is assumed to be defined elsewhere.
 * - The function does not currently handle the case where the number of variables in NGvars exceeds the number of predefined variable names.
 */
void NeuronGrowth::SaveNGvars(const vector<vector<float>>& NGvars, int NX, int NY, const string& fn) {
	bool visualization = true; // Flag to enable/disable visualization

	// Array of variable names corresponding to NGvars indices
	const char* varNames[] = {"phi", "syn", "tub", "theta", "phi_0", "tub_0"};

	// Loop through NGvars and save each to a text file with a corresponding name
	for (size_t i = 0; i < NGvars.size() && i < sizeof(varNames)/sizeof(varNames[0]); ++i) {
		string filepath = fn + "/" + varNames[i] + "_" + std::to_string(n) + ".txt";
		PrintVec2TXT(NGvars[i], filepath, visualization);
	}
}

/**
 * Prints a 3D representation of neurons to standard output using PETSc.
 *
 * This function visualizes a 3D grid of neurons, where each neuron's presence
 * is indicated by a hash ('#') symbol and absence by a space (' '). It's designed
 * to work with PETSc for parallel processing environments, outputting the neuron
 * grid to all processes in the PETSc communicator PETSC_COMM_WORLD. The function
 * iterates through a 3D vector of integers, where each integer represents the
 * presence (non-zero) or absence (zero) of a neuron at that grid location.
 *
 * @param neurons A 3D vector of integers representing the neuron grid. A non-zero
 *        value indicates the presence of a neuron, while zero indicates absence.
 *        The outer vector represents the depth (z-axis), each inner vector represents
 *        a row (y-axis), and the innermost vector represents columns (x-axis) in the
 *        3D space.
 *
 * Usage:
 *     vector<vector<vector<int>>> neuronGrid = {{{1, 0}, {0, 1}}, {{0, 1}, {1, 0}}};
 *     NeuronGrowth ng;
 *     ng.PrintOutNeurons3D(neuronGrid);
 *
 * Note:
 *     - The function uses PETSc's PetscPrintf for output, which synchronizes output
 *       across the PETSc communicator PETSC_COMM_WORLD.
 *     - The `dwnRatio` parameter is currently set to 1 and used to control downsampling,
 *       effectively printing every neuron. Adjust `dwnRatio` to skip neurons and reduce
 *       output density if needed.
 *     - Ensure PETSc is initialized before calling this function, as it uses PETSc's
 *       communication routines.
 */
void NeuronGrowth::PrintOutNeurons3D(vector<vector<vector<int>>> neurons) 
{
	ierr = PetscPrintf(PETSC_COMM_WORLD, "-----------------------------------------------------------------------------------------------\n");
	int dwnRatio = 1;

	for (int i = 0; i < neurons.size(); i += 2 * dwnRatio) {
		ierr = PetscPrintf(PETSC_COMM_WORLD, "| ");

		for (int j = 0; j < neurons[i].size(); j += dwnRatio) {
			for (int k = 0; k < neurons[i][j].size(); k += dwnRatio) {
				if (round(neurons[i][j][k]) == 0) {
					ierr = PetscPrintf(PETSC_COMM_WORLD, " ");
				} else {
					ierr = PetscPrintf(PETSC_COMM_WORLD, "#");
				}
			}
		}

		ierr = PetscPrintf(PETSC_COMM_WORLD, "|\n");
	}

	ierr = PetscPrintf(PETSC_COMM_WORLD, "-----------------------------------------------------------------------------------------------\n");
}

/**
 * @brief Computes the residual of the phase field equation and assembles the Jacobian matrix for a given state vector.
 * 
 * This function calculates the residual vector F of the phase field equation for neuron growth modeling, based on the current state vector x. It employs a novel 3D neuron growth model to simulate intricate neurite outgrowth, efficiently handling complex neurite structures through adaptive domain expansion and localized refinement using truncated hierarchical B-splines (THB-splines). The function integrates over each element of a neuron mesh, applying the phase field method and isogeometric analysis (IGA) to evaluate neuron growth dynamics. The resulting residual and Jacobian matrix are crucial for advancing the solution towards neuron morphology understanding and neurodevelopmental disorder treatment development.
 *
 * @param snes The nonlinear solver context.
 * @param x The current solution vector (input).
 * @param F The vector to store the calculated residual (output).
 * @param ctx User-defined context for neuron growth modeling, containing mesh and material properties.
 * @return PetscErrorCode Returns 0 on success, or a nonzero error code on failure.
 *
 * The function performs the following operations:
 * 1. Scatters the global vector x to a sequential vector for localized computation.
 * 2. Initializes and clears necessary data structures for computation.
 * 3. Loops through each element of the neuron mesh to:
 *    a. Prepare element matrices and vectors based on current state values.
 *    b. Evaluate orientation and other properties at Gaussian quadrature points.
 *    c. Assemble local contributions to the global residual vector F and Jacobian matrix.
 * 4. Applies boundary conditions as necessary.
 * 5. Assembles the final global residual vector and Jacobian matrix.
 * 6. Cleans up temporary data structures and restores arrays.
 *
 * @note This function is tailored for use within a PETSc-based finite element method framework for simulating neuron growth.
 */
PetscErrorCode FormFunction_phi(SNES snes, Vec x, Vec F, void *ctx)
{
	PetscErrorCode ierr;
	Vec P_seq;
	PetscScalar *Parray;
	VecScatter scatter_ctx1;
	ierr = VecScatterCreateToAll(x, &scatter_ctx1, &P_seq); CHKERRQ(ierr);
	ierr = VecScatterBegin(scatter_ctx1, x, P_seq, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecScatterEnd(scatter_ctx1, x, P_seq, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecGetArray(P_seq, &Parray); CHKERRQ(ierr);
	ierr = VecSet(F, 0.0); CHKERRQ(ierr);

	NeuronGrowth *user = (NeuronGrowth *)ctx;

	user->pre_mag_grad_phi0.clear();

	/*Build linear system in each process*/
	int ind(0); // pre-calculated variable index
	for (int e = 0; e < user->bzmesh_process.size(); e++) {
		// std::cout << e << std::endl;
		int nen = user->bzmesh_process[e].IEN.size(); // 64 supporting cp for 3D case

		vector<float> eleMphi;
		eleMphi.resize(nen);

		vector<vector<float>> EMatrixSolve;	EMatrixSolve.clear();
		vector<float> EVectorSolve;		EVectorSolve.clear();
		EMatrixSolve.resize(nen);
		EVectorSolve.resize(nen);

		user->eleVal.clear();
		user->eleVal.resize(10);
		for (int i = 0; i < nen * 1; i++) {
			EMatrixSolve[i].resize(nen);
			if (i < 10)
				user->eleVal[i].resize(nen);
		}

		// preparing element stiffness matrix and vector, extract values from control mesh
		// if (user->n >= 50)
		// 	std::cout << user->phi.size() << " " << user->syn.size() << " " << user->tub.size() << " " << user->theta.size() << " " << user->tips.size() << std::endl;

		for (int i = 0; i < nen; i++) {

			EVectorSolve[i] = 0.0;
			for (int j = 0; j < nen; j++) {
				EMatrixSolve[i][j] = 0.0;
			}

			user->eleVal[0][i] = Parray[user->bzmesh_process[e].IEN[i]];		// elePhiGuess
			user->eleVal[1][i] = user->phi[user->bzmesh_process[e].IEN[i]];		// elePhi
			user->eleVal[2][i] = user->syn[user->bzmesh_process[e].IEN[i]];		// eleSyn
			user->eleVal[3][i] = user->tub[user->bzmesh_process[e].IEN[i]];		// eleTubulin
			user->eleVal[4][i] = user->theta[user->bzmesh_process[e].IEN[i]];	// eleTheta
			user->eleVal[5][i] = user->tips[user->bzmesh_process[e].IEN[i]];	// eleTips
			user->eleVal[6][i] = 0;							// eleEpsilon
			user->eleVal[7][i] = 0;							// eleEpsilonP
			// user->eleVal[8][i] = user->polar[user->bzmesh_process[e].IEN[i]];	// elePolar
			// user->eleVal[9][i] = user->azimuth[user->bzmesh_process[e].IEN[i]];	// eleAzimuth
			// eleMphi[i] = user->Mphi[user->bzmesh_process[e].IEN[i]];
		}
		
		// loop through gaussian quadrature points
		for (int i = 0; i < user->Gpt.size(); i++) {
			for (int j = 0; j < user->Gpt.size(); j++) {
				for (int k = 0; k < user->Gpt.size(); k++) {
					
					// user->EvaluateOrientation(nen, user->pre_Nx[ind], user->pre_dNdx[ind], user->eleVal[1], user->eleVal[4], user->eleVal[6], user->eleVal[7]);
					float eleAniso(0), dA_dPdx(0), dA_dPdy(0), dA_dPdz(0);
					if (user->n > 0)
						user->EvaluateOrientation(nen, user->pre_Nx[ind], user->pre_dNdx[ind], user->eleVal[1], user->eleVal[4], eleAniso, dA_dPdx, dA_dPdy, dA_dPdz);

					// float eleEp(0), dEdp(0), dEda(0);
					// user->EvaluateOrientationSpherical(nen, user->pre_Nx[ind], user->pre_dNdx[ind], user->eleVal[1], user->eleVal[8], user->eleVal[9], eleEp, dEdp, dEda);

					user->ElementEvaluationAll_phi(nen, user->pre_Nx[ind], user->pre_dNdx[ind], user->eleVal, user->vars);
					//	 0   1    2     3      4      5     6     7      8     9      10      11      12     13    14    15      16     17   	18     19   20  
					// float C1, C0, elePG, dPGdx, dPGdy, eleP, eleS, eleTb, eleE, eleTp, dThedx, dThedy, eleEP, dAdx, dAdy, eleEEP, dAPdx, dAPdy, dPGdz, dAdz, dAPdz
					
					float eleMp;
					// user->ElementValue(user->pre_Nx[ind], eleMphi, eleMp);
					eleMp = 1;

					// float ck(0);
					// if (user->n < 300) {
					// 	ck = 0.2;
					// } else {
					// 	ck = 0;
					// }

					// adjust rg (assembly rate) and sg (disassembly rate) based on detected tips
					if (user->n < 0) {
						user->vars[8] = user->alphaOverPi*atan(user->gamma * (1 - user->vars[6]));
					} else {
						// if (user->vars[9] > ck) {
						if (user->vars[9] > 0) {
							// user->vars[8] = user->alphaOverPi*atan(user->gamma * user->Regular_Heiviside_fun(50 * user->vars[7] - 0) * (1 - user->vars[6]));
							// user->vars[8] = user->alphaOverPi*atan(user->gamma * 1 * (1 - abs(user->vars[6])));
							user->vars[8] = user->alphaOverPi*atan(user->gamma * 1 * (1 - user->vars[6]));
						} else {
							// user->vars[8] = user->alphaOverPi*atan(user->gamma * user->Regular_Heiviside_fun(user->r * user->vars[7] - user->g) * (1 - user->vars[6]));
							// user->vars[8] = user->alphaOverPi*atan(user->gamma * 0 * (1 - abs(user->vars[6])));
							user->vars[8] = user->alphaOverPi*atan(user->gamma * 0 * (1 - user->vars[6]));
						}
					}
					
					// if (user->n > 50) {
					// 	std::cout << user->alphaOverPi*atan(user->gamma * (1 - user->vars[6])) << " " << std::endl;
					// 	std::cout << user->alphaOverPi*atan(user->gamma * user->Regular_Heiviside_fun(50 * user->vars[7] - 0) * (1 - user->vars[6])) << " " << user->Regular_Heiviside_fun(50 * user->vars[7] - 0) << " " <<  user->vars[7] << std::endl;
					// 	std::cout << user->alphaOverPi*atan(user->gamma * user->Regular_Heiviside_fun(5 * user->vars[7] - user->g) * (1 - user->vars[6])) << " " << user->Regular_Heiviside_fun(50 * user->vars[7] - 0) << " " <<  user->vars[7] << std::endl;
					// }

					// calculate C1 variable for phase field energy term
					user->vars[0] = user->vars[8] - user->pre_C0[ind];

					// if (user->n > 50) {
					// 	std::cout << user->vars[0] << " " << user->vars[8] << " " << user->pre_C0[ind] << std::endl;

					// user->vars[0] = user->vars[8] - user->pre_C0_sp[ind];
					// float tau = sqrt(pow(user->vars[3],2) + pow(user->vars[4],2));
					// float mode_grad_p_2 = pow(sqrt(pow(user->vars[3],2) + pow(user->vars[4],2) + pow(user->vars[18],2)), 2);
					// loop through control points
					for (int m = 0; m < nen; m++) {
						// terma2 = - eleEP * eleEP * (dPGdx * dNdx[m][0] + dPGdy * dNdx[m][1]);
						// termadx = dAdx * eleEEP * dPGdx * dNdx[m][0];
						// termady = dAdy * eleEEP * dPGdy * dNdx[m][1];
						// termadz = dAdz * eleEEP * dPGdz * dNdx[m][2];
						// termdbl = (- elePG * elePG * elePG + (1 - C1) * elePG * elePG + C1 * elePG) * Nx[m];
						// EVectorSolve[m] += (elePG * Nx[m] - user->dt * user->M_phi * (terma2 - termadx + termady + termdbl) - eleP * Nx[m]) * detJ;
						
						EVectorSolve[m] += (user->vars[2] * user->pre_Nx[ind][m] - user->dt * eleMp * (
							(- eleAniso * eleAniso * (user->vars[3] * user->pre_dNdx[ind][m][0] + user->vars[4] * user->pre_dNdx[ind][m][1] + user->vars[18] * user->pre_dNdx[ind][m][2]))
							+ (- user->pre_dNdx[ind][m][0] * eleAniso * dA_dPdx * ( pow(user->vars[3], 2) + pow(user->vars[4], 2) + pow(user->vars[18], 2)))
							+ (- user->pre_dNdx[ind][m][1] * eleAniso * dA_dPdy * ( pow(user->vars[3], 2) + pow(user->vars[4], 2) + pow(user->vars[18], 2)))
							+ (- user->pre_dNdx[ind][m][2] * eleAniso * dA_dPdz * ( pow(user->vars[3], 2) + pow(user->vars[4], 2) + pow(user->vars[18], 2)))
							+ (- user->vars[2] * user->vars[2] * user->vars[2] + (1 - user->vars[0]) * user->vars[2] * user->vars[2] + user->vars[0] * user->vars[2]) * user->pre_Nx[ind][m]
							// - ((6 * user->vars[2] - 6 * (user->vars[2] * user->vars[2])) * user->pre_C0[ind]) * user->pre_Nx[ind][m]
							// - ((pow(user->vars[2], 2) - 2 * pow(user->vars[2], 3) + pow(user->vars[2], 4)) * user->pre_C0[ind]) * user->pre_Nx[ind][m]
							) - user->vars[5] * user->pre_Nx[ind][m]
							) * user->pre_detJ[ind];

						// EVectorSolve[m] += (user->vars[2] * user->pre_Nx[ind][m] - user->dt * eleMp * (
						// 	(- eleEp * eleEp * (user->vars[3] * user->pre_dNdx[ind][m][0] + user->vars[4] * user->pre_dNdx[ind][m][1] + user->vars[18] * user->pre_dNdx[ind][m][2])) -
						// 	(- user->pre_dNdx[ind][m][0] * eleEp / tau * dEdp * user->pre_dNdx[ind][m][0] * user->pre_dNdx[ind][m][2] -
						// 		eleEp / pow(tau, 2) * dEda * mode_grad_p_2 * user->pre_dNdx[ind][m][1]) + 
						// 	(- user->pre_dNdx[ind][m][1] * eleEp / tau * dEdp * user->pre_dNdx[ind][m][1] * user->pre_dNdx[ind][m][2] +
						// 		eleEp / pow(tau, 2) * dEda * mode_grad_p_2 * user->pre_dNdx[ind][m][1]) + 
						// 	(eleEp * dEdp * tau) * user->pre_dNdx[ind][m][2] +
						// 	(- user->vars[2] * user->vars[2] * user->vars[2] + (1 - user->vars[0]) * user->vars[2] * user->vars[2] + user->vars[0] * user->vars[2]) * user->pre_Nx[ind][m])
						// 	- user->vars[5] * user->pre_Nx[ind][m]) * user->pre_detJ[ind];

						// loop through 16 control points
						for (int n = 0; n < nen; n++) {
							// terma2 = - eleEP * eleEP * (dNdx[m][0] * dNdx[n][0] + dNdx[m][1] * dNdx[n][1]);
							// termadx = - dAdx * eleEEP * dNdx[m][1] * dNdx[n][0];
							// termady = - dAdy * eleEEP * dNdx[m][0] * dNdx[n][1];
							// termdbl = (- 3 * elePG * elePG + 2 * (1 - C1) * elePG + C1 * Nx[m]) * Nx[n];
							// EMatrixSolve[m][n] += (Nx[m] * Nx[n] - user->dt * user->M_phi * (terma2 - termadx + termady + termdbl)) * detJ;
							
							EMatrixSolve[m][n] += (user->pre_Nx[ind][m] * user->pre_Nx[ind][n] - user->dt * eleMp * (
								(- eleAniso * eleAniso * (user->pre_dNdx[ind][m][0] * user->pre_dNdx[ind][n][0] + user->pre_dNdx[ind][m][1] * user->pre_dNdx[ind][n][1] + user->pre_dNdx[ind][m][2] * user->pre_dNdx[ind][n][2])) // terma2
								+ (- user->pre_dNdx[ind][m][0] * eleAniso * dA_dPdx * ( 2 * user->vars[3] + 2 * user->vars[4] + 2 * user->vars[18]))
								+ (- user->pre_dNdx[ind][m][1] * eleAniso * dA_dPdy * ( 2 * user->vars[3] + 2 * user->vars[4] + 2 * user->vars[18]))
								+ (- user->pre_dNdx[ind][m][2] * eleAniso * dA_dPdz * ( 2 * user->vars[3] + 2 * user->vars[4] + 2 * user->vars[18]))
								+ (- 3 * user->vars[2] * user->vars[2] + 2 * (1 - user->vars[0]) * user->vars[2] + user->vars[0] * user->pre_Nx[ind][m]) * user->pre_Nx[ind][n] // termdbl
								// - ((6 * user->pre_Nx[ind][m] - 12 * user->vars[2]) * user->pre_C0[ind]) * user->pre_Nx[ind][n]
								// - ((2 * user->vars[2] - 6 * pow(user->vars[2], 2) + 4 * pow(user->vars[2], 3)) * user->pre_C0[ind]) * user->pre_Nx[ind][n])
								)) * user->pre_detJ[ind];

							// EMatrixSolve[m][n] += (user->pre_Nx[ind][m] * user->pre_Nx[ind][n] - user->dt * eleMp * (
							// 	(- eleEp * eleEp * (user->vars[3] * user->pre_dNdx[ind][m][0] + user->vars[4] * user->pre_dNdx[ind][m][1] + user->vars[18] * user->pre_dNdx[ind][m][2])) -
							// 	(- user->pre_dNdx[ind][m][0] * eleEp / tau * dEdp * user->vars[3] * user->vars[18] -
							// 		eleEp / pow(tau, 2) * dEda * (2 * user->vars[3] + 2 * user->vars[4] + 2 * user->vars[18]) * user->vars[4]) + 
							// 	(- user->pre_dNdx[ind][m][1] * eleEp / tau * dEdp * user->vars[4] * user->vars[18] +
							// 		eleEp / pow(tau, 2) * dEda * (2 * user->vars[3] + 2 * user->vars[4] + 2 * user->vars[18]) * user->vars[3]) + 
							// 	+ (- 3 * user->vars[2] * user->vars[2] + 2 * (1 - user->vars[0]) * user->vars[2] + user->vars[0] * user->pre_Nx[ind][m]) * user->pre_Nx[ind][n]) // termdbl
							// 	) * user->pre_detJ[ind];
						}
					}
					ind += 1; // incrementing index for extracting pre-calculated variables
				}
			}
		}

		/*Apply Boundary Condition*/
		for (int i = 0; i < nen; i++) {
			int A = user->bzmesh_process[e].IEN[i];
			if (user->cpts[A].label == 1) {// domain boundary
				user->ApplyBoundaryCondition(0, i, 0, EMatrixSolve, EVectorSolve);
			}
		}

		user->ResidualAssembly(EVectorSolve, user->bzmesh_process[e].IEN, F);
		user->MatrixAssembly(EMatrixSolve, user->bzmesh_process[e].IEN, user->J);
	}

	ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(F); CHKERRQ(ierr);

	ierr = MatAssemblyBegin(user->J, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(user->J, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

	ierr = VecRestoreArray(P_seq, &Parray); CHKERRQ(ierr);
	ierr = VecScatterDestroy(&scatter_ctx1); CHKERRQ(ierr);
	ierr = VecDestroy(&P_seq); CHKERRQ(ierr);
	
	return 0;
}

// Function documentation:
// Forms the Jacobian matrix for the phi component of the NeuronGrowth model.
//
// Parameters:
// - snes: The nonlinear solver context.
// - x: The current solution vector.
// - J: The Jacobian matrix.
// - P: The preconditioner matrix (usually the same as J).
// - ctx: User-defined context, cast to NeuronGrowth structure.
//
// Returns:
// - PetscErrorCode: Error code indicating success or the nature of any failure.
PetscErrorCode FormJacobian_phi(SNES snes, Vec x, Mat J, Mat P, void *ctx)
{
	NeuronGrowth *user = (NeuronGrowth *)ctx; // Cast context to NeuronGrowth structure
	PetscErrorCode ierr;

	// Copy the predefined Jacobian from the user context to the PETSc Jacobian
	ierr = MatCopy(user->J, J, SAME_NONZERO_PATTERN); CHKERRQ(ierr);
	
	return PETSC_SUCCESS; // Indicate successful execution
}

// Function documentation:
// Custom SNES monitor function that displays the default short summary of the SNES progress
// and additionally reports the number of KSP iterations.
//
// Parameters:
// - snes: The nonlinear solver context.
// - its: The number of nonlinear iterations completed.
// - fnorm: The norm of the current function value.
// - vf: Viewer and format context for output.
//
// Returns:
// - PetscErrorCode: Error code indicating success or the nature of any failure.
PetscErrorCode MySNESMonitor(SNES snes, PetscInt its, PetscReal fnorm, PetscViewerAndFormat *vf)
{
	PetscFunctionBeginUser; // PETSc macro for initiating user functions

	// Call the default short summary monitor
	SNESMonitorDefaultShort(snes, its, fnorm, vf);
	
	KSP ksp;
	SNESGetKSP(snes, &ksp); // Retrieve the KSP solver associated with the SNES solver

	// Uncomment the desired monitoring option:
	// PetscOptionsSetValue(NULL, "-ksp_monitor", "");  // Monitor every KSP iteration
	PetscOptionsSetValue(NULL, "-ksp_monitor_singular_value", "");  // Monitor singular values
	// PetscOptionsSetValue(NULL, "-ksp_monitor_solution", "");  // Monitor the solution

	PetscInt numIterations;
	KSPGetTotalIterations(ksp, &numIterations); // Retrieve the total number of iterations from the KSP solver
	PetscPrintf(PETSC_COMM_WORLD, "     - KSP Iterations: %d\n", numIterations); // Print the number of iterations

	PetscFunctionReturn(PETSC_SUCCESS); // Indicate successful execution
}

// Function documentation:
// Cleans up PETSc solvers and related resources for a NeuronGrowth object.
// This includes destroying SNES, KSP, matrices, and vectors associated with
// different aspects of neuron growth modeling (phi, syn, tub).
//
// Parameters:
// - NG: Reference to a NeuronGrowth object containing solver objects and resources.
//
// Returns:
// - PetscErrorCode: Error code indicating if the cleanup was successful or if an error occurred.
PetscErrorCode CleanUpSolvers(NeuronGrowth &NG)
{
	// Destroy SNES solver for phi
	NG.ierr = SNESDestroy(&NG.snes_phi); CHKERRQ(NG.ierr);
	
	// Destroy matrices and vectors associated with phi
	NG.ierr = MatDestroy(&NG.J); CHKERRQ(NG.ierr);
	NG.ierr = VecDestroy(&NG.temp_phi); CHKERRQ(NG.ierr);

	// Destroy KSP solver and resources for synaptogenesis (syn)
	NG.ierr = KSPDestroy(&NG.ksp_syn); CHKERRQ(NG.ierr);
	NG.ierr = MatDestroy(&NG.GK_syn); CHKERRQ(NG.ierr);
	NG.ierr = VecDestroy(&NG.GR_syn); CHKERRQ(NG.ierr);
	NG.ierr = VecDestroy(&NG.temp_syn); CHKERRQ(NG.ierr);

	// Destroy KSP solver and resources for tubules (tub)
	NG.ierr = KSPDestroy(&NG.ksp_tub); CHKERRQ(NG.ierr);	
	NG.ierr = MatDestroy(&NG.GK_tub); CHKERRQ(NG.ierr);
	NG.ierr = VecDestroy(&NG.GR_tub); CHKERRQ(NG.ierr);
	NG.ierr = VecDestroy(&NG.temp_tub); CHKERRQ(NG.ierr);

	// Return the last error code (0 if no errors occurred)
	return NG.ierr;
}

int RunNG(int n_bzmesh, vector<vector<int>> ele_process_in, vector<Vertex3D> cpts_initial, vector<Vertex3D> &cpts, vector<Vertex3D> prev_cpts, string path_in, string path_out, int &iter, int end_iter_in,
	vector<vector<float>> &NGvars, int &NX, int &NY, int &NZ, vector<array<float, 3>> &seed, int &originX, int &originY, int &originZ, bool &localRefine)
{	
	/*========================================================*/
	// Initializations
	NeuronGrowth NG;
	// // NG.SetVariables("simulation_parameters.txt");
	NG.GaussInfo(3);
	NG.n = iter;
	NG.numNeuron = seed.size();
	NG.end_iter = end_iter_in;

	NG.InitializeProblemNG(n_bzmesh, cpts, prev_cpts, NGvars, seed);
	NG.ToPETScVec(NG.phi, NG.temp_phi); // initial guess for SNES (optional)
	PetscPrintf(PETSC_COMM_WORLD, "Set initial guess!-----------------------------------------------------------\n");	
	
	NG.AssignProcessor(ele_process_in);
	// Check MPI element assignments, and print out in orders
	for (int i = 0; i < NG.nProcess; i++){
		if (i == NG.comRank) {
			std::cout << "comRank: " << NG.comRank << "/" << NG.nProcess << " - with element process size: " << NG.ele_process.size() << std::endl;
			NG.ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(NG.ierr);
		} else {
			NG.ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(NG.ierr);
		}
	}				
	// Read bezier mesh and prepare SNES initial guess
	NG.ReadBezierElementProcess(path_in);
	PetscPrintf(PETSC_COMM_WORLD, "Read bzmesh!-----------------------------------------------------------------\n");	

	NG.CheckVar(path_out + "/distI_", cpts, NG.distI);

	/*========================================================*/
	// Initial neuron identifications and geodesic distance calculation
	vector<vector<vector<int>>> neurons, distances;
	vector<float> id, tips;
	// NG.IdentifyNeurons3D(neurons, seed, NX, NY, NZ, originX, originY, originZ);
	// distances = NG.CalculateGeodesicDistanceFromPoint3D(neurons, seed, originX, originY, originZ);
	// PetscPrintf(PETSC_COMM_WORLD, "Calculated geodesic distances!-----------------------------------------------\n");
	// id = Convert3DIntTo1DFloatVector(neurons);
	// NG.DetectTipsMulti3D(id, NG.numNeuron, tips, NX, NY, NZ);
	// // NG.tips = InterpolateVars3D(tips, cpts_initial, cpts, 0);	
	// NG.tips = InterpolateValues3D(cpts_initial, tips, cpts);	
	
	// PetscPrintf(PETSC_COMM_WORLD, "Detected initial tips!-------------------------------------------------------\n");

	/*========================================================*/
	// Write initial variables
	string varName;	
	if (NG.n == 0) {
		NG.VisualizeVTK_PhysicalDomain_All(0, path_out);
		PetscPrintf(PETSC_COMM_WORLD, "Saving all variables!--------------------------------------------------------\n");		
	}

	/*========================================================*/
	// calculate some variable in advance to save computational cost
	NG.prepareBasis();
	PetscPrintf(PETSC_COMM_WORLD, "Prepared basis!--------------------------------------------------------------\n");	
	NG.ierr = MPI_Allreduce(&NG.sum_grad_phi0_local, &NG.sum_grad_phi0_global, 1, MPI_FLOAT, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(NG.ierr);
	PetscPrintf(PETSC_COMM_WORLD, "Calculated sum grad phi0!----------------------------------------------------\n");
	NG.prepareTerm_source();
	// NG.prepareEE();
	PetscPrintf(PETSC_COMM_WORLD, "Prepared variables!----------------------------------------------------------\n");

	while (iter <= NG.end_iter) {
		NG.n = iter;
		tic();		

		/*========================================================*/
		// Neuron identification and tip detection
		// if (NG.n < 100 ) {	

		// 	if ((NG.n % 25 == 0) || (NG.n == 0) || (NG.tips.size() != NG.phi.size())) {	
		// 		PetscPrintf(PETSC_COMM_WORLD, "-----------------------------------------------------------------------------------------\n");	
		// 		PetscPrintf(PETSC_COMM_WORLD, "Identifying neurons, calculating geodesic distances, and detecting tips\n");
		// 		// NG.tips = NG.calculatePhiSum(cpts, 2, 2, 2);
		// 		// NG.tips = NG.calculatePhiSum(cpts, 4, 4, 4);
		// 		NG.calculatePhiSum(cpts, 4, 4, 4);
		// 		NG.CheckVar(path_out + "/tips_", cpts, NG.tips);
		// 		// NG.CheckVar("../io3D/outputs/tips_", cpts, tips);

		// 		toc(t_collect);
		// 		tic();
		// 		t_write += t_collect;
		// 		t_total += t_collect;
		// 	}
		// } else {
		// 	if ((NG.n % 25 == 0) || (NG.n == 0) || (NG.tips.size() != NG.phi.size())) {	
		// 		PetscPrintf(PETSC_COMM_WORLD, "-----------------------------------------------------------------------------------------\n");	
		// 		PetscPrintf(PETSC_COMM_WORLD, "Identifying neurons, calculating geodesic distances, and detecting tips\n");
		// 		// NG.tips = NG.calculatePhiSum(cpts, 2, 2, 2);
		// 		// NG.tips = NG.calculatePhiSum(cpts, 4, 4, 4);
		// 		NG.calculatePhiSum(cpts, 4, 4, 4);
		// 		NG.CheckVar(path_out + "/tips_", cpts, NG.tips);
		// 		// NG.CheckVar("../io3D/outputs/tips_", cpts, tips);

		// 		toc(t_collect);
		// 		tic();
		// 		t_write += t_collect;
		// 		t_total += t_collect;
		// 	}
		// }

		// if ((NG.n % NG.var_save_invl == 0) && (NG.n != 0)) {	
		// if ((NG.n % 25 == 0) && (NG.n != 0)) {
		if ((NG.n % 25 == 0) || (NG.n == 0) || (NG.tips.size() != NG.phi.size())) {	
			PetscPrintf(PETSC_COMM_WORLD, "-----------------------------------------------------------------------------------------\n");	
			
			PetscPrintf(PETSC_COMM_WORLD, "Identifying neurons, calculating geodesic distances, and detecting tips\n");
			// NG.tips = NG.calculatePhiSum(cpts, 2, 2, 2);
			// NG.tips = NG.calculatePhiSum(cpts, 4, 4, 4);
			NG.calculatePhiSum(cpts, 8, 8, 8);
			NG.CheckVar(path_out + "/tips_", cpts, NG.tips);
			// NG.CheckVar("../io3D/outputs/tips_", cpts, tips);

			toc(t_collect);
			tic();
			t_write += t_collect;
			t_total += t_collect;

			/*========================================================*/
			// Writing results to files
			// PetscPrintf(PETSC_COMM_WORLD, "-----------------------------------------------------------------------------------------\n");
			NG.VisualizeVTK_PhysicalDomain_All(NG.n, path_out);
			PetscPrintf(PETSC_COMM_WORLD, "Step: %d/%d | Wrote Physical Domain! | Average time %fs | Total time: %f |\n", 
				NG.n, NG.end_iter, t_write/NG.var_save_invl, t_total); CHKERRQ(NG.ierr);
			PetscPrintf(PETSC_COMM_WORLD, "-----------------------------------------------------------------------------------------\n");
			t_write = 0;
		}

		/*========================================================*/
		/*Implcit Non-liear SNES solver for Phase field equation*/
		NG.phi_prev = NG.phi;		

		if (NG.judge_phi == 0) { 
			NG.ierr = SNESCreate(PETSC_COMM_WORLD, &NG.snes_phi); CHKERRQ(NG.ierr);
			NG.ierr = SNESSetType(NG.snes_phi, SNESNEWTONLS); CHKERRQ(NG.ierr);
			NG.ierr = SNESSetTolerances(NG.snes_phi, 1e-5, 1e-7, 1e-9, 100, 1000); CHKERRQ(NG.ierr);
			// NG.ierr = SNESSetTolerances(NG.snes_phi, 1e-4, 1e-6, 1e-8, 100, 1000); CHKERRQ(NG.ierr);
			PetscPrintf(PETSC_COMM_WORLD, "Set Tolerance!---------------------------------------------------------------\n");
			SNESLineSearch linesearch;
			NG.ierr = SNESGetLineSearch(NG.snes_phi, &linesearch); CHKERRQ(NG.ierr);
			NG.ierr = SNESLineSearchSetType(linesearch, SNESLINESEARCHCP); CHKERRQ(NG.ierr);
			PetscPrintf(PETSC_COMM_WORLD, "Set LineSearch!--------------------------------------------------------------\n");	
			NG.ierr = SNESSetFunction(NG.snes_phi, NULL, FormFunction_phi, &NG); CHKERRQ(NG.ierr);
			NG.ierr = SNESSetJacobian(NG.snes_phi, NULL, NULL, FormJacobian_phi, &NG); CHKERRQ(NG.ierr);
			PetscPrintf(PETSC_COMM_WORLD, "Set FormFunction and FormJacobian!-------------------------------------------\n");

			// PetscViewerAndFormat *vf;
			// PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_DEFAULT, &vf);
			// SNESMonitorSet(NG.snes_phi, (PetscErrorCode(*)(SNES, PetscInt, PetscReal, void *))MySNESMonitor, vf, (PetscErrorCode(*)(void **))PetscViewerAndFormatDestroy);
			// PetscPrintf(PETSC_COMM_WORLD, "Set Monitor!-----------------------------------------------------------------\n");
			
			if (NG.n == 0)
				NG.ierr = SNESView(NG.snes_phi, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(NG.ierr);
			NG.judge_phi = 1;
		}

		NG.ierr = SNESSolve(NG.snes_phi, NULL, NG.temp_phi); CHKERRQ(NG.ierr);
		PetscInt its_phi;
		NG.ierr = SNESGetIterationNumber(NG.snes_phi, &its_phi); CHKERRQ(NG.ierr);
		SNESConvergedReason reason_phi;
		NG.ierr = SNESGetConvergedReason(NG.snes_phi, &reason_phi); CHKERRQ(NG.ierr);
		if (reason_phi < 0) {
			PetscPrintf(PETSC_COMM_WORLD, "SNES phi not converging ........ | SNESreason: %d | SNESiter: %d \n", reason_phi, its_phi); CHKERRQ(NG.ierr);	
			return 3; // divering, ending simulation
		}

		/*========================================================*/
		/*Collecting scattered Phi variable from all processors*/
		Vec temp_phi_seq;
		VecScatter ctx_phi;
		PetscScalar *_p;

		NG.ierr = VecScatterCreateToAll(NG.temp_phi, &ctx_phi, &temp_phi_seq); CHKERRQ(NG.ierr);
		NG.ierr = VecScatterBegin(ctx_phi, NG.temp_phi, temp_phi_seq, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(NG.ierr);
		NG.ierr = VecScatterEnd(ctx_phi, NG.temp_phi, temp_phi_seq, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(NG.ierr);
		NG.ierr = VecGetArray(temp_phi_seq, &_p);

		for (int i = 0; i < NG.phi.size(); i++) {
			NG.phi[i] = PetscRealPart(_p[i]);
		}

		NG.ierr = VecRestoreArray(temp_phi_seq, &_p); CHKERRQ(NG.ierr);
		NG.ierr = VecScatterDestroy(&ctx_phi); CHKERRQ(NG.ierr);
		NG.ierr = VecDestroy(&temp_phi_seq); CHKERRQ(NG.ierr);

		toc(t_phi);
		tic();
		t_write += t_phi;
		t_total += t_phi;
		
		/*========================================================*/
		/*Synaptogenesis and Tubulin equation*/
		NG.ierr = MPI_Allreduce(&NG.sum_grad_phi0_local, &NG.sum_grad_phi0_global, 1, MPI_FLOAT, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(NG.ierr);
		NG.prepareTerm_source();
		// NG.prepareEE();
	
		/*Build Linear System*/
		if (NG.judge_syn == 0) {
			NG.ierr = MatZeroEntries(NG.GK_syn); CHKERRQ(NG.ierr);
		}
		NG.ierr = VecSet(NG.GR_syn, 0); CHKERRQ(NG.ierr);
		NG.ierr = MatZeroEntries(NG.GK_tub); CHKERRQ(NG.ierr);
		NG.ierr = VecSet(NG.GR_tub, 0); CHKERRQ(NG.ierr);
		NG.BuildLinearSystemProcessNG_syn_tub(cpts); // build both Synaptogenesis and tubulin in the same loop
		NG.ierr = VecAssemblyEnd(NG.GR_syn); CHKERRQ(NG.ierr);
		if (NG.judge_syn == 0) {
			NG.ierr = MatAssemblyEnd(NG.GK_syn, MAT_FINAL_ASSEMBLY); CHKERRQ(NG.ierr);
		}
		NG.ierr = VecAssemblyEnd(NG.GR_tub); CHKERRQ(NG.ierr);
		NG.ierr = MatAssemblyEnd(NG.GK_tub, MAT_FINAL_ASSEMBLY); CHKERRQ(NG.ierr);

		/*Petsc solver setting for Synaptogenesis*/
		if (NG.judge_syn == 0) {
			NG.ierr = KSPCreate(PETSC_COMM_WORLD, &NG.ksp_syn); CHKERRQ(NG.ierr);
			NG.ierr = KSPSetOperators(NG.ksp_syn, NG.GK_syn, NG.GK_syn); CHKERRQ(NG.ierr);
			NG.ierr = KSPGetPC(NG.ksp_syn, &NG.pc_syn); CHKERRQ(NG.ierr);
			NG.ierr = PCSetType(NG.pc_syn, PCBJACOBI); CHKERRQ(NG.ierr);
			NG.ierr = KSPSetType(NG.ksp_syn, KSPGMRES); CHKERRQ(NG.ierr);
			// NG.ierr = KSPSetType(NG.ksp_syn, KSPCGS); CHKERRQ(NG.ierr);
			NG.ierr = KSPGMRESSetRestart(NG.ksp_syn, 100); CHKERRQ(NG.ierr);
			NG.ierr = KSPSetInitialGuessNonzero(NG.ksp_syn, PETSC_TRUE); CHKERRQ(NG.ierr);
			NG.ierr = KSPSetTolerances(NG.ksp_syn, 1.e-8, PETSC_DEFAULT, PETSC_DEFAULT, 100000); CHKERRQ(NG.ierr);
			NG.ierr = KSPSetFromOptions(NG.ksp_syn); CHKERRQ(NG.ierr);
			NG.ierr = KSPSetUp(NG.ksp_syn); CHKERRQ(NG.ierr);
			if (NG.n == 0)
				NG.ierr = KSPView(NG.ksp_syn, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(NG.ierr);
			NG.judge_syn = 1;
		}

		NG.ierr = KSPSolve(NG.ksp_syn, NG.GR_syn, NG.temp_syn); CHKERRQ(NG.ierr);
		PetscInt its_syn;
		NG.ierr = KSPGetIterationNumber(NG.ksp_syn, &its_syn); CHKERRQ(NG.ierr);
		KSPConvergedReason reason_syn;
		NG.ierr = KSPGetConvergedReason(NG.ksp_syn, &reason_syn); CHKERRQ(NG.ierr);
		if (reason_syn < 0) {
			PetscPrintf(PETSC_COMM_WORLD, "KSP Synaptogenesis not converging  | KSPreason:  %d | KSPiter: %d \n", reason_syn, its_syn); CHKERRQ(NG.ierr);	
			return 3; // divering, ending simulation
		}

		toc(t_syn);
		tic();
		t_write += t_syn;
		t_total += t_syn;

		/*Petsc solver setting for Tubulin*/
		if (NG.judge_tub == 0) {
			NG.ierr = KSPCreate(PETSC_COMM_WORLD, &NG.ksp_tub); CHKERRQ(NG.ierr);
			NG.ierr = KSPSetOperators(NG.ksp_tub, NG.GK_tub, NG.GK_tub); CHKERRQ(NG.ierr);
			NG.ierr = KSPGetPC(NG.ksp_tub, &NG.pc_tub); CHKERRQ(NG.ierr);
			NG.ierr = PCSetType(NG.pc_tub, PCBJACOBI); CHKERRQ(NG.ierr);
			NG.ierr = KSPSetType(NG.ksp_tub, KSPGMRES); CHKERRQ(NG.ierr);
			// NG.ierr = KSPSetType(NG.ksp_tub, KSPCGS); CHKERRQ(NG.ierr);
			NG.ierr = KSPGMRESSetRestart(NG.ksp_tub, 100); CHKERRQ(NG.ierr);
			NG.ierr = KSPSetInitialGuessNonzero(NG.ksp_tub, PETSC_TRUE); CHKERRQ(NG.ierr);
			NG.ierr = KSPSetTolerances(NG.ksp_tub, 1.e-8, PETSC_DEFAULT, PETSC_DEFAULT, 100000); CHKERRQ(NG.ierr);
			NG.ierr = KSPSetFromOptions(NG.ksp_tub); CHKERRQ(NG.ierr);
			NG.ierr = KSPSetUp(NG.ksp_tub); CHKERRQ(NG.ierr);
			if (NG.n == 0)
				NG.ierr = KSPView(NG.ksp_tub, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(NG.ierr);
			NG.judge_tub = 1;
		}

		NG.ierr = KSPSolve(NG.ksp_tub, NG.GR_tub, NG.temp_tub); CHKERRQ(NG.ierr);
		PetscInt its_tub;
		NG.ierr = KSPGetIterationNumber(NG.ksp_tub, &its_tub); CHKERRQ(NG.ierr);
		KSPConvergedReason reason_tub;
		NG.ierr = KSPGetConvergedReason(NG.ksp_tub, &reason_tub); CHKERRQ(NG.ierr);
		if (reason_tub < 0) {
			PetscPrintf(PETSC_COMM_WORLD, "KSP tubulin not converging ..... | KSPreason:  %d | KSPiter: %d \n", reason_tub, its_tub); CHKERRQ(NG.ierr);		
			return 3; // divering, ending simulation
		}

		toc(t_tub);
		tic();
		t_write += t_tub;
		t_total += t_tub;

		/*========================================================*/
		/*Collecting scattered Synaptogenesis and Tubulin variableS from all processors*/
		Vec temp_syn_seq, temp_tub_seq;
		VecScatter ctx_syn, ctx_tub;
		PetscScalar *_s, *_t;

		NG.ierr = VecScatterCreateToAll(NG.temp_syn, &ctx_syn, &temp_syn_seq); CHKERRQ(NG.ierr);
		NG.ierr = VecScatterBegin(ctx_syn, NG.temp_syn, temp_syn_seq, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(NG.ierr);
		NG.ierr = VecScatterEnd(ctx_syn, NG.temp_syn, temp_syn_seq, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(NG.ierr);
		NG.ierr = VecGetArray(temp_syn_seq, &_s); CHKERRQ(NG.ierr);

		NG.ierr = VecScatterCreateToAll(NG.temp_tub, &ctx_tub, &temp_tub_seq); CHKERRQ(NG.ierr);
		NG.ierr = VecScatterBegin(ctx_tub, NG.temp_tub, temp_tub_seq, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(NG.ierr);
		NG.ierr = VecScatterEnd(ctx_tub, NG.temp_tub, temp_tub_seq, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(NG.ierr);
		NG.ierr = VecGetArray(temp_tub_seq, &_t); CHKERRQ(NG.ierr);

		for (int i = 0; i < NG.syn.size(); i++) {
			NG.syn[i] = PetscRealPart(_s[i]);
			NG.tub[i] = PetscMax(PetscRealPart(_t[i]), 0.0f) * NG.CellBoundary(NG.phi[i], 0.5);
			
			if (isnan(NG.tub[i]) || NG.tub[i] > 1)
				NG.tub[i] = 0;
		}

		NG.ierr = VecRestoreArray(temp_syn_seq, &_s); CHKERRQ(NG.ierr);
		NG.ierr = VecScatterDestroy(&ctx_syn); CHKERRQ(NG.ierr);
		NG.ierr = VecDestroy(&temp_syn_seq); CHKERRQ(NG.ierr);

		NG.ierr = VecRestoreArray(temp_tub_seq, &_t); CHKERRQ(NG.ierr);
		NG.ierr = VecScatterDestroy(&ctx_tub); CHKERRQ(NG.ierr);
		NG.ierr = VecDestroy(&temp_tub_seq); CHKERRQ(NG.ierr);

		PetscPrintf(PETSC_COMM_WORLD, "Step:%d/%d | Phi:%d[%d]%.3fs | Syn:%d[%d]%.3fs | Tub:%d[%d]%.3fs | Mesh:%d |\n",\
			NG.n, NG.end_iter, reason_phi, its_phi, t_phi, reason_syn, its_syn, t_syn, reason_tub, its_tub, t_tub, n_bzmesh); CHKERRQ(NG.ierr);

		/*========================================================*/
		// Obtain initial local refinement information, the very first 25 iterations are purely used 
		// for getting diffused interface for applying local refinements (phi initialization is binary) 
		if ((NG.n == 10) && (localRefine == false)) {
			// NGvars.clear(); NGvars.resize(6);
			// NGvars[0] = NG.phi;	
			// NGvars[1] = NG.syn;
			// NGvars[2] = NG.tub;
			// NGvars[3] = NG.theta;
			// NGvars[4] = NG.phi_0;
			// NGvars[5] = NG.tub_0;

			CleanUpSolvers(NG); // Destroy solvers
			NG.ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(NG.ierr);
			iter = 0; // reset to 0 (beginning of the simulation)
			localRefine = true;

			if (NG.comRank == 0) {
				vector<float> ele_refine = ComputeRefine(NG.phi, NX, NY, NZ);
				writeVectorToFile(ele_refine, path_in + "phi.txt", false);
				NG.CheckVar("../io3D/phi", cpts, NG.phi);
			}	
			NG.ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(NG.ierr);
			// std::cout <<"!!!!!!!" << std::endl;

			return 2;
		}
		
		// /*========================================================*/
		// // Neuron identification and tip detection
		// if ((NG.n % 1 == 0) && (NG.n != 0)) {	
		// 	// PetscPrintf(PETSC_COMM_WORLD, "-----------------------------------------------------------------------------------------\n");	
		// 	// PetscPrintf(PETSC_COMM_WORLD, "Identifying neurons, calculating geodesic distances, and detecting tips\n");
		// 	// std::cout << "ck0" << std::endl;
		// 	// NG.IdentifyNeurons3D(neurons, seed, NX, NY, NZ, originX, originY, originZ);

		// 	// std::cout << "ck1" << std::endl;
		// 	// vector<float> id = Convert3DIntTo1DFloatVector(neurons);
		// 	// std::cout << NX << " " <<  NY << " " << NZ << " " << id.size() << " " << cpts_initial.size() << std::endl;
		// 	// NG.CheckVar("../io3D/id", cpts_initial, id);
		// 	// std::cout << "ck2" << std::endl;
		// 	// NG.DetectTipsMulti3D(id, NG.numNeuron, NG.tips, NX, NY, NZ);
		// 	NG.tips = NG.calculatePhiSum(cpts, 4, 4, 4);
		// 	// std::cout << "ck3" << std::endl;

		// 	// // vector<float> localMaximaMatrix = NG.FindLocalMaximaInClusters(tips, NX+1, NY+1);
		// 	// // NG.tips = InterpolateVars(localMaximaMatrix, cpts_initial, cpts, 0);
		// 	// distances = NG.CalculateGeodesicDistanceFromPoint3D(neurons, seed, originX, originY, originZ);
		// 	// vector<float> geodist = Convert3DIntTo1DFloatVector(distances);

		// 	// NG.CheckVar("../io3D/tips", cpts, NG.tips);
		// 	// NG.CheckVar("../io3D/geodist", cpts, geodist);

		// 	vector<float> Mphi; Mphi.clear(); Mphi.resize(NG.tips.size());
		// 	// int maxGeoInd(0);
		// 	for (int i = 0; i < NG.tips.size(); i++) {
		// 		// if ((geodist[i] != INF) && (geodist[i] > maxGeoInd)) {
		// 		// 	maxGeoInd = i;
		// 		// }
		// 		if (localRefine == true)
		// 			if (NG.tips[i] == 0) {
		// 				Mphi[i] = 0;
		// 			} else {
		// 				Mphi[i] = 1;
		// 			}
		// 		else {
		// 			Mphi[i] = 1;
		// 		}

		// 		// Mphi[i] = 10;

		// 	}
		// 	// for (int i = -1; i < 1; i++) { // increase Mphi at longest tip 
		// 	// 	for (int j = -1; j < 1; j++) {
		// 	// 		Mphi[maxGeoInd+j*(NY+1)-i] = 15;
		// 	// 	}
		// 	// }
		// 	// NG.Mphi = InterpolateVars(Mphi, cpts_initial, cpts, 0);
		// 	// NG.Mphi = InterpolateValues3D(cpts_initial, Mphi, cpts);	
		// 	NG.Mphi = Mphi;

		// 	toc(t_collect);
		// 	tic();
		// 	t_write += t_collect;
		// 	t_total += t_collect;
		// }

		/*========================================================*/
		// Domain expansion and variable passing - back to main.cpp
		if ((NG.n % NG.expandCK_invl == 0) && (NG.n != 0)) {

			localRefine = true;
			// // NG.PrintOutNeurons(neurons); 
			// // neurons = ConvertTo3DIntVector(NG.phi, NX, NY, NZ);
			// int expd_dir_local = NG.CheckExpansion3D(NG.phi, cpts, NX, NY, NZ, originX, originY, originZ);
			// // 0 - left | 1 - top | 2 - right | 3 - bottom | 4 - front | 5 - back | 6 - no action
			// NG.ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(NG.ierr);
			// int expd_dir_global = 6; // 6 - No expansion
			// NG.ierr = MPI_Allreduce(&expd_dir_local, &expd_dir_global, 1, MPI_INT, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(NG.ierr);
			// int expd_dir = expd_dir_global/NG.comSize; // approximate choice, different thread could have different choice 

			// if (expd_dir <= 3) {
			// 	switch (expd_dir) {
			// 		case 0: // left
			// 			NX += 3;	originX -= 1.5/4;	break;
			// 		case 1: // top
			// 			NY += 3;				break;
			// 		case 2: // right
			// 			NX += 3;				break;
			// 		case 3: // bottom
			// 			NY += 3;	originY -= 1.5/4;	break;
			// 		case 4: // back
			// 			NZ += 3;				break;
			// 		case 5: // front
			// 			NZ += 3;	originZ -= 1.5/4;	break;
			// 	}
			// }

			NGvars.clear(); NGvars.resize(6);
			NGvars[0] = NG.phi;	
			NGvars[1] = NG.syn;
			NGvars[2] = NG.tub;
			NGvars[3] = NG.theta;
			NGvars[4] = NG.phi_0;
			NGvars[5] = NG.tub_0;

			if (NG.comRank == 0) {
				vector<float> ele_refine = ComputeRefine(NGvars[0], NX, NY, NZ);
				writeVectorToFile(ele_refine, path_in + "phi.txt", false);
				NG.CheckVar("../io3D/phi", cpts, NG.phi);
			}
			
			CleanUpSolvers(NG); // Destroy solvers
			NG.ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(NG.ierr);
			iter += 1;
			return 2;

		}
		iter += 1;		 
	}

	/*========================================================*/
	CleanUpSolvers(NG); // Destroy solvers
	NG.ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(NG.ierr);

	return 0; // ending simulation
}