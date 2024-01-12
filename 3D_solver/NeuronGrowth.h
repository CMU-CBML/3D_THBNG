#ifndef NeuronGrowth_H
#define NeuronGrowth_H

#include <vector>
#include <array>
#include "BasicDataStructure.h"
#include "utils.h"
#include "time.h"

using namespace std;

// timing function (similar to Matlab tic toc)
void tic();
void toc(float &t);

float MatrixDet(float dxdt[2][2]);
void Matrix2DInverse(float dxdt[2][2], float dtdx[2][2]);

class NeuronGrowth
{
private:

public:
	// MPI parameters
	PetscErrorCode ierr;
	MPI_Comm comm;
	int mpiErr;
	int comRank;
	int comSize;
	int nProcess;

	// Spline parameters
	int n_bzmesh;
	vector<int> ele_process;
	vector<float> Gpt, wght, N_0;
	vector<Vertex3D> cpts;
	vector<Element3D> bzmesh_process;

	// Pre-calculated variables to save computational cost
	vector<vector<float>> pre_Nx;
	vector<vector<array<float, 2>>> pre_dNdx;
	vector<float> pre_detJ, pre_mag_grad_phi0, pre_C0, pre_term_source;
	vector<float> pre_eleEP, pre_eleEEP, pre_dAdx, pre_dAdy, pre_dAPdx, pre_dAPdy;
	vector<float> pre_vars;
	vector<vector<vector<float>>> pre_EMatrixSolve;
	vector<vector<float>> pre_EVectorSolve;

	// element stiffness matrix and load vector
	int nen;
	vector<vector<float>> EMatrixSolve;
	vector<float> EVectorSolve;
	vector<vector<float>> eleVal;
	vector<float> vars;
	vector<float> Nx;
	vector<array<float, 2>> dNdx;

	// Neuron growth variables
	int n; 								// time step
	int judge_phi, judge_syn, judge_tub;				// assembly state	
	vector<float> phi, tub, syn, theta;				// variable to be solved 
	vector<float> phi_prev, phi_0, tub_0, tips, Mphi;			// assisting varibles
	float sum_grad_phi0_local, sum_grad_phi0_global, dP0dx, dP0dy, dP0dz;	
	vector<float> elePhi0, eleTheta;

	// PETSc solvers and variables
	SNES snes_phi;				// PETSc SNES nonlinear solver
	KSP ksp_syn, ksp_tub;			// PETSc KSP linear solver
	PC pc_syn, pc_tub;			// PETSc preconditioner
	Mat GK_syn, GK_tub, J;			// Jacobian matrix
	Vec GR_syn, GR_tub;			// Residual vector
	Vec temp_phi, temp_syn, temp_tub;	// Solution vector

	// Parameters for neuron growth model
	int var_save_invl,expandCK_invl,numNeuron,gc_sz,end_iter,aniso,gamma,seed_radius;
	float kappa,dt,Dc,alpha,alphaOverPi,M_phi,s_coeff,delta,epsilonb,r,g,alphaT,betaT,Diff,source_coeff;

	// Initializations
	NeuronGrowth();
	void AssignProcessor(vector<vector<int>> &ele_proc); // assign elements to different processors
	void SetVariables(string fn_par);
	void InitializeProblemNG(const int n_bz, vector<Vertex3D>& cpts, vector<Vertex3D> prev_cpts, vector<vector<float>> &NGvars, vector<array<int, 2>> &seed);
	void ToPETScVec(vector<float> input, Vec& petscVec); // for SNES phi initial guess

	// Read mesh, calculate basis function value, assemble matrix, etc
	void ReadBezierElementProcess(string fn);
	void GaussInfo(int ng);
	void BasisFunction(float u, float v, const vector<array<float, 2>>& pt, const vector<array<float, 16>> &cmat,
		vector<float> &Nx, vector<array<float, 2>> &dNdx, float dudx[2][2], float& detJ);
	void WeightingFunction(const float velocity[2], const float& s, const float& tau, const vector<float> &Nx, const vector<array<float, 2>> &dNdx, vector<float> &Wx);
	void ApplyBoundaryCondition(const float bc_value, int pt_num, int variable_num, vector<vector<float>>& EMatrixSolve, vector<float>& EVectorSolve);
	void MatrixAssembly(vector<vector<float>>& EMatrixSolve, const vector<int>& IEN, Mat& GK);
	void ResidualAssembly(vector<float>& EVectorSolve, const vector<int>& IEN, Vec& GR);
	void MatrixAssembly_insert(vector<vector<float>>& EMatrixSolve, const vector<int>& IEN, Mat& GK);
	void ResidualAssembly_insert(vector<float>& EVectorSolve, const vector<int>& IEN, Vec& GR);
	void MatrixAssembly_2var(vector<vector<float>>& EMatrixSolve, const vector<int>& IEN, Mat& GK);
	void ResidualAssembly_2var(vector<float>& EVectorSolve, const vector<int>& IEN, Vec& GR);

	// writing files
	void VisualizeVTK_ControlMesh(const vector<Vertex3D> &spt, const vector<Element3D> &mesh, int step, string fn, vector<float> var, string varName);
	void ConcentrationCal_Coupling_Bezier(double u, double v, double w, const Element3D& bzel, double pt[3], double& disp, double dudx[3], double& detJ);
	void VisualizeVTK_PhysicalDomain(int step, string var, string fn);
	void WriteVTK(const vector<array<double, 3>> spt, const vector<double> sdisp, const vector<array<int, 8>> sele, int step, string fn);
	void CalculateVarsForOutput(vector<array<float, 3>> &spt_all, vector<float> &sresult_all, vector<array<int, 8>> &sele_all);
	void VisualizeVTK_PhysicalDomain_All(int step, string fn);
	void WriteVTK_All(const vector<array<float, 3>> spt, const vector<vector<float>> sdisp, const vector<array<int, 8>> sele, int step, string fn);

	// element based operation
	void PointFormValue(vector<float> &Nx, const vector<float> &U, float Value);
	void PointFormGrad(vector<array<float, 2>> &dNdx, const vector<float> &U, float Value[2]);
	void PointFormHess(vector<array<array<float, 2>, 2>>& d2Ndx2, const vector<float> &U, float Value[2][2]);
	void ElementValue(const vector<float> &Nx, const vector<float> value_node, float &value);
	void ElementValueAll(const vector<float> &Nx, const vector<float> elePhiGuess, float &elePG,
		const vector<float> elePhi, float &eleP,
		const vector<float> eleSyn, float &eleS,
		const vector<float> eleTips, float &eleTp,
		const vector<float> eleTubulin, float &eleTb,
		const vector<float> eleEpsilon, float &eleEP,
		const vector<float> eleEpsilonP, float &eleEEP);
	void ElementDeriv(const int nen, vector<array<float, 3>> &dNdx, const vector<float> value_node, float &dVdx, float &dVdy, float &dVdz);
	void ElementDerivAll(const int nen, vector<array<float, 2>> &dNdx,
		const vector<float> elePhiGuess, float &dPGdx, float &dPGdy,
		const vector<float> eleTheta, float &dThedx, float &dThedy,
		const vector<float> eleEpsilon, float &dAdx, float &dAdy,
		const vector<float> eleEpsilonP, float &dAPdx, float &dAPdy);
	void ElementEvaluationAll_phi(const int nen, const vector<float> &Nx, vector<array<float, 3>> &dNdx,
		const vector<float> elePhiGuess, float &elePG,
		const vector<float> elePhi, float &eleP,
		const vector<float> eleSyn, float &eleS,
		const vector<float> eleTips, float &eleTp,
		const vector<float> eleTubulin, float &eleTb,
		float &dPGdx, float &dPGdy, float &dPGdz);
	void ElementEvaluationAll_phi(const int nen, const vector<float> &Nx, vector<array<float, 2>> &dNdx,
		const vector<float> elePhiGuess, float &elePG,
		const vector<float> elePhi, float &eleP,
		const vector<float> eleSyn, float &eleS,
		const vector<float> eleTips, float &eleTp,
		const vector<float> eleTubulin, float &eleTb,
		const vector<float> eleEpsilon, float &eleEP,
		const vector<float> eleEpsilonP, float &eleEEP,
		float &dPGdx, float &dPGdy,
		float &dAdx, float &dAdy,
		float &dAPdx, float &dAPdy);
	void ElementEvaluationAll_phi(const int nen, const vector<float> &Nx, vector<array<float, 3>> &dNdx,
		vector<vector<float>> &eleVal, vector<float> &vars)
	void ElementEvaluationAll_syn_tub(const int nen, const vector<float> &Nx, vector<array<float, 3>> &dNdx,
		const vector<float> elePhiDiff, float &elePf, 
		const vector<float> eleSyn,float &eleS, const vector<float> elePhi, float &eleP,
		const vector<float> elePhiPrev, float &elePprev,
		const vector<float> eleConct, float &eleC,
		float &dPdx, float &dPdy, float &dPdz)
	void ElementEvaluationAll_syn_tub(const int nen, const vector<float> &Nx, vector<array<float, 3>> &dNdx,
		vector<vector<float>> &eleVal, vector<float> &vars);

	// pre-calculate variables to save computational cost
	void prepareBasis();
	void prepareTerm_source();
	void prepareEpsilon();
	void prepareEE();
	// void prepareVars();

	// Phase field equation
	void EvaluateEnergy(const int nen, const vector<float> &Nx, const vector<float> eleS, vector<float>& E);
	float Regular_Heiviside_fun(float x);
	void EvaluateOrientation(const int nen, const vector<float> &Nx, const vector<array<float, 2>> &dNdx, const vector<float> elePhi,
		const vector<float> eleTheta,  vector<float>& eleEpsilon, vector<float>& eleEpsilonP);

	// Synaptogenesis equation
	void BuildLinearSystemProcessNG_syn(const vector<Element3D>& tmesh, const vector<Vertex3D> &cpts);
	void Tangent_syn(const int nen, vector<float>& Nx, vector<array<float, 2>>& dNdx, float detJ, vector<vector<float>>& EMatrixSolve);
	void Residual_syn(const int nen, vector<float>& Nx, vector<array<float, 2>>& dNdx, float detJ, float eleP, float eleS, vector<float>& EVectorSolve);
	// Tubulin equation
	void CalculateSumGradPhi0(const vector<Element3D> &tmesh, const vector<Vertex3D> &cpts);
	void BuildLinearSystemProcessNG_tub(const vector<Element3D>& tmesh, const vector<Vertex3D> &cpts);
	void Tangent_tub(const int nen, vector<float>& Nx, vector<array<float, 2>>& dNdx, float detJ, float eleC, float eleP, float dPdx, float dPdy,
		float elePprev, float mag_grad_phi0, vector<vector<float>>& EMatrixSolve);
	void Residual_tub(const int nen, vector<float>& Nx, vector<array<float, 2>>& dNdx, float detJ, float eleC, float eleP, float elePprev,
		float mag_grad_phi0, vector<float>& EVectorSolve);
	// Build Synaptogenesis and Tubulin together
	void BuildLinearSystemProcessNG_syn_tub(const vector<Element3D> &tmesh, const vector<Vertex3D> &cpts);

	// Domain expansion
	int CheckExpansion(vector<float> input); // for when square domain
	int NeuronGrowth::CheckExpansion3D(vector<float> input, int NX, int NY, int NZ);
	void ExpandDomain(vector<float> input, vector<float> &expd_var, int edge); // for square domain
	void ExpandDomain(vector<float> input, vector<float> &expd_var, int edge, int NX, int NY); // for when NX != NY
	void PopulateRandom(vector<float> &input); // to populate theta with random after expansion

	// Tip detection
	float RmOutlier(vector<float> &data); // standard deviation based outlier remover
	float CellBoundary(float phi, float threshold); // threshould based boundary determination
	vector<float> SmoothBinary2D(const std::vector<float>& binaryData, int rows, int cols);
	void DetectTipsSquare(vector<float> phi, vector<float> &tips); // for square domain
	void DetectTipsRectangle(vector<float> phi, vector<float> &tips, int NX, int NY); // for when NX != NY
	void DetectTipsMulti3D(vector<float> id, int numNeuron, vector<float> &tips, int NX, int NY, int NZ);
	void bfs(const std::vector<float>& matrix, int rows, int cols, int row, int col,
		std::vector<bool>& visited, std::vector<std::pair<int, int>>& cluster);
	std::vector<std::vector<std::pair<int, int>>> FindClusters(const std::vector<float>& matrix, int rows, int cols);
	std::vector<float> FindLocalMaximaInClusters(const std::vector<float>& matrix, int rows, int cols);

	// Neuron detection
	vector<vector<int>> ConvertTo2DIntVector(const vector<float> input, int NX, int NY);
	vector<vector<float>> ConvertTo2DFloatVector(const vector<float> input, int NX, int NY);
	
	void FloodFill3D(std::vector<std::vector<std::vector<int>>>& image, int x, int y, int z, int newColor, int originalColor);
	void IdentifyNeurons3D(std::vector<std::vector<std::vector<int>>>& neurons, std::vector<std::array<int, 3>> seed, int NX, int NY, int NZ, int originX, int originY, int originZ);
	bool isValid(int x, int y, int z, int rows, int cols, int depth);
	std::vector<std::vector<std::vector<int>>> CalculateGeodesicDistanceFromPoint3D(std::vector<std::vector<std::vector<int>>> neurons, const std::vector<std::array<int, 3>>& seed, int originX, int originY, int originZ);
	vector<vector<array<int, 3>>> NeuriteTracing(vector<vector<double>> distance);
	void SaveNGvars(vector<vector<float>> &NGvars, int NX, int NY, string fn);
	void PrintOutNeurons3D(vector<vector<vector<int>>> neurons);
};

// Phase field PETSc Nonlinear SNES solver functions (placing here due to non-static member function error)
PetscErrorCode FormFunction_phi(SNES snes, Vec x, Vec F, void *ctx);
PetscErrorCode FormJacobian_phi(SNES snes, Vec x, Mat J, Mat P, void *ctx);
PetscErrorCode MySNESMonitor(SNES snes, PetscInt its, PetscReal fnorm, PetscViewerAndFormat *vf);
PetscErrorCode CleanUpSolvers(NeuronGrowth &NG);

// Main function
int RunNG(int n_bzmesh, vector<vector<int>> ele_process_in, vector<Vertex3D> cpts_initial, const vector<Element3D> tmesh_initial, vector<Vertex3D> &cpts, vector<Vertex3D> prev_cpts, const vector<Element3D> &tmesh, string path_in, string path_out, int &iter, int end_iter_in,
	vector<vector<float>> &NGvars, int &NX, int &NY, int &NZ, vector<array<int, 3>> &seed, int &originX, int &originY, int &originZ, bool &localRefine);

#endif