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

	float max_x, min_x, max_y, min_y, max_z, min_z;

	// Pre-calculated variables to save computational cost
	vector<vector<float>> pre_Nx;
	vector<vector<array<float, 3>>> pre_dNdx;
	vector<float> pre_detJ, pre_mag_grad_phi0, pre_C0, pre_C0_sp, pre_term_source;
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
	vector<array<float, 3>> dNdx;

	// Neuron growth variables
	int n; 								// time step
	int judge_phi, judge_syn, judge_tub;				// assembly state	
	vector<float> phi, tub, syn, theta;				// variable to be solved 
	//  polar, azimuth;
	vector<float> phi_prev, phi_0, tub_0, tips, Mphi;			// assisting varibles
	float sum_grad_phi0_local, sum_grad_phi0_global, dP0dx, dP0dy, dP0dz;	
	vector<float> elePhi0, eleTheta;

	vector<float> distI;

	// PETSc solvers and variables
	SNES snes_phi;				// PETSc SNES nonlinear solver
	KSP ksp_syn, ksp_tub;			// PETSc KSP linear solver
	PC pc_syn, pc_tub;			// PETSc preconditioner
	Mat GK_syn, GK_tub, J;			// Jacobian matrix
	Vec GR_syn, GR_tub;			// Residual vector
	Vec temp_phi, temp_syn, temp_tub;	// Solution vector

	// Parameters for neuron growth model
	int var_save_invl,expandCK_invl,numNeuron,gc_sz,end_iter,aniso,gamma;
	float seed_radius, kappa,dt,Dc,alpha,alphaOverPi,M_phi,s_coeff,delta,epsilonb,r,g,alphaT,betaT,Diff,source_coeff;

	// Initializations
	NeuronGrowth();
	void AssignProcessor(vector<vector<int>> &ele_proc); // assign elements to different processors
	void SetVariables(string fn_par);
	void InitializeProblemNG(const int n_bz, vector<Vertex3D>& cpts, vector<Vertex3D> prev_cpts, vector<vector<float>> &NGvars, vector<array<float, 3>> &seed);
	void CheckVar(string fn, vector<Vertex3D> cpts, vector<float> input);
	void ToPETScVec(vector<float> input, Vec& petscVec); // for SNES phi initial guess

	// Read mesh, calculate basis function value, assemble matrix, etc
	void ReadBezierElementProcess(string fn);
	void GaussInfo(int ng);
	void BasisFunction(float u, float v, float w, const vector<array<float, 3>>& pt, const vector<array<float, 64>> &cmat,
		vector<float> &Nx, vector<array<float, 3>> &dNdx, float dudx[3][3], float& detJ);
	void ApplyBoundaryCondition(const float bc_value, int pt_num, int variable_num, vector<vector<float>>& EMatrixSolve, vector<float>& EVectorSolve);
	void MatrixAssembly(vector<vector<float>>& EMatrixSolve, const vector<int>& IEN, Mat& GK);
	void ResidualAssembly(vector<float>& EVectorSolve, const vector<int>& IEN, Vec& GR);
	void MatrixAssembly_insert(vector<vector<float>>& EMatrixSolve, const vector<int>& IEN, Mat& GK);
	void ResidualAssembly_insert(vector<float>& EVectorSolve, const vector<int>& IEN, Vec& GR);
	void MatrixAssembly_2var(vector<vector<float>>& EMatrixSolve, const vector<int>& IEN, Mat& GK);
	void ResidualAssembly_2var(vector<float>& EVectorSolve, const vector<int>& IEN, Vec& GR);

	// writing files
	void VisualizeVTK_ControlMesh(const vector<Vertex3D> &spt, const vector<Element3D> &mesh, int step, string fn, vector<float> var, string varName);
	void ConcentrationCal_Coupling_Bezier(float u, float v, float w, const Element3D& bzel, float pt[3], float& disp, float dudx[3], float& detJ);
	void VisualizeVTK_PhysicalDomain(int step, string var, string fn);
	void WriteVTK(const vector<array<float, 3>> spt, const vector<float> sdisp, const vector<array<int, 8>> sele, int step, string fn);
	void CalculateVarsForOutput(vector<array<float, 3>> &spt_all, vector<float> &sresult_all, vector<array<int, 8>> &sele_all);
	void VisualizeVTK_PhysicalDomain_All(int step, string fn);
	void WriteVTK_ALL(const vector<array<float, 3>> spt, const vector<vector<float>> sdisp, const vector<array<int, 8>> sele, int step, string fn);

	// element based operation
	void PointFormValue(vector<float> &Nx, const vector<float> &U, float Value);
	void PointFormGrad(vector<array<float, 3>> &dNdx, const vector<float> &U, float Value[2]);
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
	void ElementDerivAll(const int nen, vector<array<float, 3>> &dNdx,
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
	void ElementEvaluationAll_phi(const int nen, const vector<float> &Nx, vector<array<float, 3>> &dNdx,
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
		vector<vector<float>> &eleVal, vector<float> &vars);
	void ElementEvaluationAll_syn_tub(const int nen, const vector<float> &Nx, vector<array<float, 3>> &dNdx,
		const vector<float> elePhiDiff, float &elePf, 
		const vector<float> eleSyn,float &eleS, const vector<float> elePhi, float &eleP,
		const vector<float> elePhiPrev, float &elePprev,
		const vector<float> eleConct, float &eleC,
		float &dPdx, float &dPdy, float &dPdz);
	void ElementEvaluationAll_syn_tub(const int nen, const vector<float> &Nx, vector<array<float, 3>> &dNdx,
		vector<vector<float>> &eleVal, vector<float> &vars);

	// pre-calculate variables to save computational cost
	void prepareBasis();
	void prepareTerm_source();
	void prepareEpsilon();
	// void prepareEE();

	// Phase field equation
	void EvaluateEnergy(const int nen, const vector<float> &Nx, const vector<float> eleS, vector<float>& E);
	float Regular_Heiviside_fun(float x);
	void EvaluateOrientation(const int nen, const vector<float> &Nx, const vector<array<float, 3>> &dNdx, const vector<float> elePhi,
		const vector<float> eleTheta,  float& eleAniso, float& dA_dPdx, float& dA_dPdy, float& dA_dPdz);
	void EvaluateOrientationSpherical(const int nen, const vector<float> &Nx, const vector<array<float, 3>> &dNdx, const vector<float> elePhi,
		const vector<float> elePolar, const vector<float> eleAzimuth, float& eleEpsilon, float dEdp, float dEda);

	// Build Synaptogenesis and Tubulin together
	void CalculateSumGradPhi0(const vector<Vertex3D> &cpts);
	void BuildLinearSystemProcessNG_syn_tub(const vector<Vertex3D> &cpts);

	// Domain expansion
	int CheckExpansion3D(vector<float> input, const std::vector<Vertex3D>& cpts, int NX, int NY, int NZ, int originX, int originY, int originZ);
	void PopulateRandom(vector<float> &input); // to populate theta with random after expansion

	// Tip detection
	float RmOutlier(vector<float> &data); // standard deviation based outlier remover
	float CellBoundary(float phi, float threshold); // threshould based boundary determination
	void DetectTipsMulti3D(vector<float> id, int numNeuron, vector<float> &tips, int NX, int NY, int NZ);
	
	// Function to check if a point is within the specified box centered at 'center'
	bool isInBox(const Vertex3D& point, const Vertex3D& center, float dx, float dy, float dz);
	// Function to calculate the sum of phi within a specified box for each center point in cpts
	void calculatePhiSum(const std::vector<Vertex3D>& cpts, float dx, float dy, float dz);
	vector<float> InterpolateValues3D(const vector<Vertex3D>& cpts_initial, const vector<float>& input,
                                      const vector<Vertex3D>& cpts_new);
	std::vector<std::pair<Vertex3D, int>> FindClosestVerticesWithIndices(const std::vector<Vertex3D>& vertices, const Vertex3D& inputVertex, int k);
	std::vector<std::tuple<Vertex3D, int, float>> FindClosestVerticesWithIndicesAndDistances(const std::vector<Vertex3D>& vertices, const Vertex3D& inputVertex, int k);

 	void bfs3D(const vector<float>& matrix, int depth, int rows, int cols, int dep, int row, int col,
		vector<bool>& visited, vector<tuple<int, int, int>>& cluster);
	vector<vector<tuple<int, int, int>>> FindClusters3D(const vector<float>& matrix, int depth, int rows, int cols);
	vector<float> FindLocalMaximaInClusters3D(const vector<float>& matrix, int depth, int rows, int cols);

	// Neuron detection
	vector<vector<vector<int>>> ConvertTo3DIntVector(const vector<float> input, int NX, int NY, int NZ);
	vector<vector<vector<float>>> ConvertTo3DFloatVector(const vector<float> input, int NX, int NY, int NZ);

	void FloodFill3D(std::vector<std::vector<std::vector<int>>>& image, int x, int y, int z, int newColor, int originalColor);
	void IdentifyNeurons3D(std::vector<std::vector<std::vector<int>>>& neurons, std::vector<std::array<int, 3>> seed, int NX, int NY, int NZ, int originX, int originY, int originZ);
	bool isValid(int x, int y, int z, int rows, int cols, int depth);
	std::vector<std::vector<std::vector<int>>> CalculateGeodesicDistanceFromPoint3D(std::vector<std::vector<std::vector<int>>> neurons, const std::vector<std::array<int, 3>>& seed, int originX, int originY, int originZ);
	// vector<vector<array<int, 3>>> NeuriteTracing(vector<vector<float>> distance);
	void SaveNGvars(const vector<vector<float>> &NGvars, int NX, int NY, const string& fn);
	void PrintOutNeurons3D(vector<vector<vector<int>>> neurons);

	

};

// Phase field PETSc Nonlinear SNES solver functions (placing here due to non-static member function error)
PetscErrorCode FormFunction_phi(SNES snes, Vec x, Vec F, void *ctx);
PetscErrorCode FormFunction_phi_wip(SNES snes, Vec x, Vec F, void *ctx);
PetscErrorCode FormJacobian_phi(SNES snes, Vec x, Mat J, Mat P, void *ctx);
PetscErrorCode MySNESMonitor(SNES snes, PetscInt its, PetscReal fnorm, PetscViewerAndFormat *vf);
PetscErrorCode CleanUpSolvers(NeuronGrowth &NG);

int RunNG(int n_bzmesh, vector<vector<int>> ele_process_in, vector<Vertex3D> cpts_initial, vector<Vertex3D> &cpts, vector<Vertex3D> prev_cpts, string path_in, string path_out, int &iter, int end_iter_in,
	vector<vector<float>> &NGvars, int &NX, int &NY, int &NZ, vector<array<float, 3>> &seed, int &originX, int &originY, int &originZ, bool &localRefine);

#endif