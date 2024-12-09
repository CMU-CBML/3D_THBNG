#ifndef NEURONGROWTH_H
#define NEURONGROWTH_H

#include <vector>
#include <array>
#include <string>
#include <cmath>
#include <queue>
#include <numeric>
#include "BasicDataStructure.h"
#include "utils.h"
#include "../nanoflann/1.5.5/include/nanoflann.hpp" // KDTree for spatial searches

using namespace std;

// **Utility Functions**: Timing and Matrix Operations
void tic();                              // Start timing
void toc(float &t);                      // End timing and update time

float MatrixDet(float dxdt[2][2]);       // 2x2 matrix determinant
void Matrix2DInverse(float dxdt[2][2], float dtdx[2][2]); // Inverse of a 2x2 matrix

inline float SquaredDistance(const Vertex3D &first, const Vertex3D &other); // Squared Euclidean distance
void CheckAndPrintThresholdExceedance(const vector<float> &input, float threshold); // Print values exceeding a threshold

// **Vertex3DCloud Struct**: KDTree Point Cloud Wrapper
struct Vertex3DCloud {
    const vector<Vertex3D>& pts;

    explicit Vertex3DCloud(const vector<Vertex3D>& pts) : pts(pts) {}

    inline size_t kdtree_get_point_count() const { return pts.size(); }
    inline float kdtree_get_pt(const size_t idx, const size_t dim) const {
        return pts[idx].coor[dim];
    }

    template <class BBOX>
    bool kdtree_get_bbox(BBOX&) const { return false; }
};

using KDTree = nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor<float, Vertex3DCloud>, 
    Vertex3DCloud, 
    3 /* dim */>;

class NeuronGrowth {
public:
    // MPI Parameters
    PetscErrorCode ierr;       // PETSc error code
    MPI_Comm comm;             // MPI communicator
    int mpiErr;                // MPI error code
    int comRank;               // MPI rank
    int comSize;               // Number of processes in communicator
    int nProcess;              // Total processes

    // Spline Parameters
    int n_bzmesh;                      // Number of Bezier mesh elements
    vector<int> ele_process;           // Elements assigned to the process
    vector<float> Gpt, wght, N_0;      // Gauss points, weights, and basis functions
    vector<Vertex3D> cpts;             // Control points
    vector<Element3D> bzmesh_process;  // Processed Bezier elements

    // Spatial Bounds
    float max_x, min_x;  // x-dimension bounds
    float max_y, min_y;  // y-dimension bounds
    float max_z, min_z;  // z-dimension bounds

    // Pre-calculated Variables (for computational efficiency)
    vector<vector<float>> pre_Nx;
    vector<vector<array<float, 3>>> pre_dNdx;
    vector<float> pre_detJ, pre_mag_grad_phi0, pre_C0, pre_C0_sp, pre_term_source;
    vector<float> pre_eleEP, pre_eleEEP, pre_dAdx, pre_dAdy, pre_dAdz, pre_dAPdx, pre_dAPdy, pre_dAPdz;
    vector<float> pre_eleP, pre_eleTh, pre_eleMp, pre_C1;
    vector<float> pre_vars;
    vector<vector<vector<float>>> pre_EMatrixSolve;
    vector<vector<float>> pre_EVectorSolve;

    // Element Stiffness Matrix and Load Vector
    int nen;                         // Number of element nodes
    vector<vector<float>> EMatrixSolve; // Element stiffness matrix
    vector<float> EVectorSolve;        // Element load vector
    vector<vector<float>> eleVal;      // Element values
    vector<float> vars;                // Global variables
    vector<float> Nx;                  // Shape functions
    vector<array<float, 3>> dNdx;      // Shape function gradients

    // Neuron Growth Variables
    int n;                              // Current time step
    int judge_phi, judge_syn, judge_tub; // Assembly state flags
    vector<float> phi, tub, syn, theta; // Variables to be solved
    vector<float> phi_prev, phi_0, tub_0, tips, Mphi; // Supporting variables
    float sum_grad_phi0_local, sum_grad_phi0_global;  // Gradient sums
    float dP0dx, dP0dy, dP0dz;          // Derivatives of pressure
    vector<float> elePhi0, eleTheta;    // Element-specific variables
    vector<float> distI;                // Distances for interpolation

    // PETSc Solvers and Variables
    SNES snes_phi;              // PETSc SNES nonlinear solver
    KSP ksp_syn, ksp_tub;       // PETSc KSP linear solvers
    PC pc_syn, pc_tub;          // PETSc preconditioners
    Mat GK_syn, GK_tub, J;      // PETSc matrices
    Vec GR_syn, GR_tub;         // Residual vectors
    Vec temp_phi, temp_syn, temp_tub; // Temporary solution vectors

    // Parameters for Neuron Growth Model
    int var_save_invl;          // Interval for saving variables
    int expandCK_invl;          // Interval for checking expansion
    int numNeuron;              // Number of neurons
    int gc_sz;                  // Grid cell size
    int end_iter;               // Total number of iterations
    int aniso;                  // Anisotropy parameter
    int gamma;                  // Growth factor
    int seed_radius;            // Radius for neuron seeding

    float kappa;                // Diffusion coefficient
    float dt;                   // Time step
    float Dc;                   // Diffusion constant
    float kp75, k2;             // Material-specific parameters
    float c_opt;                // Optimization parameter
    float alpha;                // Growth rate
    float alphaOverPi;          // Normalized growth rate
    float M_phi, M_axon, M_neurite; // Mobility parameters
    float s_coeff;              // Coefficient for stress
    float delta;                // Growth parameter
    float epsilonb;             // Boundary thickness
    float r;                    // Radius
    float g;                    // Growth anisotropy
    float alphaT, betaT;        // Time-dependent parameters
    float Diff;                 // Diffusion
    float source_coeff;         // Source coefficient

	// Initializations
	NeuronGrowth();
	void AssignProcessor(vector<vector<int>> &ele_proc); // assign elements to different processors
	void SetVariables(string fn_par);
	void InitializeProblemNG(const int n_bz,
							vector<Vertex3D>& cpts,
							const Vertex3DCloud& cloud,
							KDTree& kdTree,
							const vector<Vertex3D>& prev_cpts,
							const Vertex3DCloud& cloud_prev,
							KDTree& kdTree_prev,
							vector<vector<float>>& NGvars,
							vector<array<float, 3>>& seed);
	void CheckVar(const string& fn, const vector<Vertex3D>& cpts, const vector<float>& input);
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
	void PrepareBasis();
	void PreparePhaseField();
	void PrepareTermSource();
	void PrepareEpsilon();
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
	int CheckExpansion3D(vector<float> input, const vector<Vertex3D>& cpts, int NX, int NY, int NZ, int originX, int originY, int originZ);
	void PopulateRandom(vector<float> &input); // to populate theta with random after expansion

	bool KD_SearchPair(const vector<Vertex3D>& cpts, 
					const KDTree& kdTree, 
					float targetX, float targetY, float targetZ, 
					int& ind, float tolerance);
	
	// Tip detection
	float RmOutlier(vector<float> &data); // standard deviation based outlier remover
	float CellBoundary(float phi, float threshold); // threshould based boundary determination
	void DetectTipsMulti3D(vector<float> id, int numNeuron, vector<float> &tips, int NX, int NY, int NZ);
	
	// Function to check if a point is within the specified box centered at 'center'
	bool IsInBox(const Vertex3D& point, const Vertex3D& center, float dx, float dy, float dz);
	// Function to calculate the sum of phi within a specified box for each center point in cpts

	void CalculatePhiSum(const vector<Vertex3D>& cpts, 
						float dx, float dy, float dz, 
						const KDTree& kdTree);
	vector<float> InterpolateValues3D(const vector<Vertex3D>& cpts_initial, const vector<float>& input,
		const vector<Vertex3D>& cpts_new);
	vector<pair<Vertex3D, int>> FindClosestVerticesWithIndices(const vector<Vertex3D>& vertices, const Vertex3D& inputVertex, int k);
	// vector<tuple<Vertex3D, int, float>> FindClosestVerticesWithIndicesAndDistances(const vector<Vertex3D>& vertices, const Vertex3D& inputVertex, int k);
	vector<tuple<Vertex3D, int, float>> FindClosestVerticesWithIndicesAndDistances(const KDTree& kdTree, const Vertex3DCloud& cloud,
																				const Vertex3D& inputVertex, int k);

 	void BFS3D(const vector<float>& matrix, int depth, int rows, int cols, int dep, int row, int col,
		vector<bool>& visited, vector<tuple<int, int, int>>& cluster);
	vector<vector<tuple<int, int, int>>> FindClusters3D(const vector<float>& matrix, int depth, int rows, int cols);
	vector<float> FindLocalMaximaInClusters3D(const vector<float>& matrix, int depth, int rows, int cols);

	// Neuron detection
    vector<vector<vector<int>>> ConvertTo3DIntVector(const vector<float>& input, int NX, int NY, int NZ);
    vector<vector<vector<float>>> ConvertTo3DFloatVector(const vector<float>& input, int NX, int NY, int NZ);

	void FloodFill3DWithKDTree(std::vector<std::vector<std::vector<int>>>& image,
							int x, int y, int z, int newColor, int originalColor,
							const KDTree& kdTree, const Vertex3DCloud& cloud);

	void IdentifyNeurons3DWithKDTree(std::vector<std::vector<std::vector<int>>>& neurons, 
									const std::vector<std::array<int, 3>>& seed,
									int NX, int NY, int NZ, 
									int originX, int originY, int originZ,
									const KDTree& kdTree, const Vertex3DCloud& cloud);							 
	bool IsValid(int x, int y, int z, int rows, int cols, int depth);
	vector<vector<vector<int>>> CalculateGeodesicDistanceFromPoint3D(vector<vector<vector<int>>> neurons, const vector<array<int, 3>>& seed, int originX, int originY, int originZ);
	// vector<vector<array<int, 3>>> NeuriteTracing(vector<vector<float>> distance);
	void SaveNGvars(const vector<vector<float>> &NGvars, int NX, int NY, const string& fn);
	void PrintOutNeurons3D(vector<vector<vector<int>>> neurons);
};

// Phase field PETSc Nonlinear SNES solver functions (placing here due to non-static member function error)
PetscErrorCode SetupSNES(SNES &snes, const char *solverType, void *ctx,
						PetscErrorCode (*formFunction)(SNES, Vec, Vec, void *),
						PetscErrorCode (*formJacobian)(SNES, Vec, Mat, Mat, void *),
						PetscReal rtol, PetscReal atol, PetscReal stol, PetscInt maxIters, PetscInt maxFails);
PetscErrorCode SetupKSP(KSP &ksp, Mat &A, const char *kspType, const char *pcType,
                        PetscReal rtol, PetscReal atol, PetscReal dtol, PetscInt maxIters, PetscInt restart);
PetscErrorCode ScatterVector(Vec src, vector<float>& target, PetscInt size, bool applyBoundary, NeuronGrowth* NG);
PetscErrorCode FormFunction_phi(SNES snes, Vec x, Vec F, void *ctx);
PetscErrorCode FormJacobian_phi(SNES snes, Vec x, Mat J, Mat P, void *ctx);
PetscErrorCode MySNESMonitor(SNES snes, PetscInt its, PetscReal fnorm, PetscViewerAndFormat *vf);
PetscErrorCode CleanUpSolvers(NeuronGrowth &NG);

int RunNG(int n_bzmesh, vector<vector<int>> ele_process_in,
		vector<Vertex3D> cpts_initial, vector<Vertex3D> &cpts, vector<Vertex3D> prev_cpts,
		string path_in, string path_out,
		int &iter, int end_iter_in,
		vector<vector<float>> &NGvars,
		int &NX, int &NY, int &NZ,
		vector<array<float, 3>> &seed, int &originX, int &originY, int &originZ,
		bool &localRefine);
#endif