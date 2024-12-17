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

	// Basis Function Values and Derivatives
	vector<vector<float>> pre_Nx;                   // Shape function values at Gauss points
	vector<vector<array<float, 3>>> pre_dNdx;      // Shape function derivatives at Gauss points

	// Jacobian and Related Properties
	vector<float> pre_detJ;                        // Determinant of Jacobian at Gauss points
	vector<float> pre_mag_grad_phi0;               // Magnitude of gradient of phi at Gauss points

	// Phase Field-Related Parameters
	vector<float> pre_C0;                          // Pre-calculated C0 values
	vector<float> pre_C0_sp;                       // Pre-calculated specific C0 values
	vector<float> pre_term_source;                 // Source term contributions
	vector<float> pre_eleEP;                       // Element epsilon values
	vector<float> pre_eleEEP;                      // Element epsilon derivative values

	// Gradients of Phase Field Variables
	vector<float> pre_dAdx;                        // Derivative of A w.r.t. x
	vector<float> pre_dAdy;                        // Derivative of A w.r.t. y
	vector<float> pre_dAdz;                        // Derivative of A w.r.t. z
	vector<float> pre_dAPdx;                       // Derivative of AP w.r.t. x
	vector<float> pre_dAPdy;                       // Derivative of AP w.r.t. y
	vector<float> pre_dAPdz;                       // Derivative of AP w.r.t. z

	// Element Properties
	vector<float> pre_eleP;                        // Pre-calculated element P values
	vector<float> pre_eleTh;                       // Pre-calculated element theta values
	vector<float> pre_eleMp;                       // Pre-calculated element mass values
	vector<float> pre_C1;                          // Pre-calculated C1 values

	// Miscellaneous Precomputed Variables
	vector<float> pre_vars;                        // General-purpose precomputed variables

	// Matrices and Vectors for Assembly
	vector<vector<vector<float>>> pre_EMatrixSolve; // Pre-calculated element stiffness matrices
	vector<vector<float>> pre_EVectorSolve;         // Pre-calculated element residual vectors

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
    KSP ksp_phi, ksp_syn, ksp_tub;       // PETSc KSP linear solvers
    PC pc_phi, pc_syn, pc_tub;          // PETSc preconditioners
	Mat J;
    Mat GK_phi, GK_syn, GK_tub;      // PETSc matrices
    Vec GR_phi, GR_syn, GR_tub;         // Residual vectors
    Vec temp_phi, temp_syn, temp_tub; // Temporary solution vectors

    // Parameters for Neuron Growth Model
    int var_save_invl;          // Interval for saving variables
    int expandCK_invl;          // Interval for checking expansion
	int refine_invl;			// Interval for local refinement
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

	string phi_solver;
	// Class Methods for NeuronGrowth
	NeuronGrowth(const string& phi_solver_in,
				const int iter_in,
				const int numNeuron_in,
				const int end_iter_in); // Constructor

	// Assign elements to different processors
	void AssignProcessor(vector<vector<int>> &ele_proc); 

	// Load simulation variables from a parameter file
	void SetVariables(const string &fn_par);

	// Initialize problem for Neuron Growth (NG) simulation
	void InitializeProblemNG(
		const int n_bz,                      // Number of Bezier elements
		vector<Vertex3D>& cpts,              // Current control points
		const Vertex3DCloud& cloud,          // Current spatial cloud
		KDTree& kdTree,                      // Current KDTree for spatial search
		const vector<Vertex3D>& prev_cpts,   // Previous control points
		const Vertex3DCloud& cloud_prev,     // Previous spatial cloud
		KDTree& kdTree_prev,                 // Previous KDTree for spatial search
		vector<vector<float>>& NGvars,       // Variables for Neuron Growth
		vector<array<float, 3>>& seed        // Seed data for initialization
	);

	void InterpolateOrFindExact(
		const Vertex3D& cpt, 
		const KDTree& kdTree_prev, 
		const Vertex3DCloud& cloud_prev, 
		const vector<vector<float>>& NGvars, 
		const vector<Vertex3D>& prev_cpts, 
		float& phi, float& syn, float& tub, float& theta, float& phi_0, float& tub_0, 
		float& dist, float& maxDist, 
		bool withinBounds
	);

	// Check and save variables to file
	void CheckVar(
		const string& fn,                    // File name for output
		const vector<Vertex3D>& cpts,        // Control points to save
		const vector<float>& input           // Variable data to save
	);

	// Convert a standard vector to a PETSc vector
	void ToPETScVec(
		vector<float> input,                 // Input standard vector
		Vec& petscVec                        // Output PETSc vector
	); // Typically used for SNES Phi initial guess

	// Methods for mesh reading, basis function calculations, and matrix assembly
	void ReadBezierElementProcess(const string& fn); // Read Bezier mesh data from file
	void GaussInfo(int ng);                          // Initialize Gaussian quadrature points and weights

	// Basis function evaluation and Jacobian computation
	void BasisFunction(
		float u, float v, float w,                   // Parametric coordinates
		const vector<array<float, 3>>& pt,           // Control points
		const vector<array<float, 64>>& cmat,        // Coefficients matrix
		vector<float>& Nx,                           // Basis function values
		vector<array<float, 3>>& dNdx,               // Basis function derivatives
		float dudx[3][3],                            // Jacobian matrix
		float& detJ                                  // Determinant of Jacobian
	);

	// Apply boundary conditions
	void ApplyBoundaryCondition(
		float bc_value,                              // Boundary condition value
		int pt_num,                                  // Point index
		int variable_num,                            // Variable index
		vector<vector<float>>& EMatrixSolve,         // Element stiffness matrix
		vector<float>& EVectorSolve                  // Element residual vector
	);

	// Assembly functions for global stiffness matrix and residual vector
	void MatrixAssembly(
		vector<vector<float>>& EMatrixSolve,         // Element stiffness matrix
		const vector<int>& IEN,                      // Global node indices
		Mat& GK                                      // Global stiffness matrix
	);
	void ResidualAssembly(
		vector<float>& EVectorSolve,                 // Element residual vector
		const vector<int>& IEN,                      // Global node indices
		Vec& GR                                      // Global residual vector
	);

	// Visualize control mesh and save as VTK file
	void VisualizeVTK_ControlMesh(
		const vector<Vertex3D>& spt,          // Control points
		const vector<Element3D>& mesh,        // Mesh elements
		int step,                             // Time step for output
		string fn,                            // File name
		vector<float> var,                    // Variable values to visualize
		string varName                        // Variable name
	);

	// Compute concentration and coupling in Bezier elements
	void ConcentrationCal_Coupling_Bezier(
		float u, float v, float w,            // Parametric coordinates
		const Element3D& bzel,                // Bezier element
		float pt[3],                          // Physical coordinates of the point
		float& disp,                          // Displacement value
		float dudx[3],                        // Derivative of displacement
		float& detJ                           // Determinant of Jacobian
	);

	// Visualize physical domain and save as VTK file
	void VisualizeVTK_PhysicalDomain(
		int step,                             // Time step for output
		string var,                           // Variable name
		string fn                             // File name
	);

	// Write VTK file with specified points, displacements, and elements
	void WriteVTK(
		const vector<array<float, 3>> spt,    // Spatial points
		const vector<float> sdisp,            // Displacement values
		const vector<array<int, 8>> sele,     // Element connectivity
		int step,                             // Time step for output
		string fn                             // File name
	);

	// Calculate variables for output
	void CalculateVarsForOutput(
		vector<array<float, 3>>& spt_all,     // All spatial points
		vector<float>& sresult_all,           // Resultant variable values
		vector<array<int, 8>>& sele_all       // Element connectivity for all elements
	);

	// Visualize the entire physical domain and save as VTK file
	void VisualizeVTK_PhysicalDomain_All(
	int step,							// Time step for output
	string fn  							// File name
	);

	// Write VTK file for all variables with points, displacements, and elements
	void WriteVTK_ALL(
		const vector<array<float, 3>> spt,    // Spatial points
		const vector<vector<float>> sdisp,    // Displacement values for multiple variables
		const vector<array<int, 8>> sele,     // Element connectivity
		int step,                             // Time step for output
		string fn                      // File name
	);

	// Evaluate the value of a field at a point using basis functions
	void PointFormValue(
		vector<float>& Nx,                     // Basis function values at the point
		const vector<float>& U,                // Nodal values of the field
		float& Value                           // Resulting field value at the point
	);

	// Evaluate the gradient of a field at a point using basis function derivatives
	void PointFormGrad(
		vector<array<float, 3>>& dNdx,         // Derivatives of basis functions
		const vector<float>& U,                // Nodal values of the field
		float Value[2]                         // Gradient of the field [dUdx, dUdy]
	);

	// Evaluate the Hessian of a field at a point using second derivatives of basis functions
	void PointFormHess(
		vector<array<array<float, 2>, 2>>& d2Ndx2, // Second derivatives of basis functions
		const vector<float>& U,                    // Nodal values of the field
		float Value[2][2]                          // Hessian of the field
	);

	// Compute a field's value at an element using basis functions
	void ElementValue(
		const vector<float>& Nx,                // Basis function values
		const vector<float>& value_node,        // Nodal values of the field
		float& value                            // Computed field value at the element
	);

	// Compute multiple field values at an element
	void ElementValueAll(
		const vector<float>& Nx,                				// Basis function values
		const vector<float>& elePhiGuess, float& elePG,       	// PhiGuess
		const vector<float>& elePhi, float& eleP,             	// Phi
		const vector<float>& eleSyn, float& eleS,             	// Syn
		const vector<float>& eleTips, float& eleTp,           	// Tips
		const vector<float>& eleTubulin, float& eleTb,        	// Tubulin
		const vector<float>& eleEpsilon, float& eleEP,        	// Epsilon
		const vector<float>& eleEpsilonP, float& eleEEP       	// Epsilon derivative
	);

	// Compute the gradient of a field at an element
	void ElementDeriv(
		const int nen,                        // Number of nodes in the element
		vector<array<float, 3>>& dNdx,        // Derivatives of basis functions
		const vector<float>& value_node,      // Nodal values of the field
		float& dVdx, float& dVdy, float& dVdz // Gradients of the field [dVdx, dVdy, dVdz]
	);

	// Compute gradients of multiple fields at an element
	void ElementDerivAll(
		const int nen,                        // Number of nodes in the element
		vector<array<float, 3>>& dNdx,        // Derivatives of basis functions
		const vector<float>& elePhiGuess, float& dPGdx, float& dPGdy,  // PhiGuess derivatives
		const vector<float>& eleTheta, float& dThedx, float& dThedy,   // Theta derivatives
		const vector<float>& eleEpsilon, float& dAdx, float& dAdy,     // Epsilon derivatives
		const vector<float>& eleEpsilonP, float& dAPdx, float& dAPdy   // EpsilonP derivatives
	);

	// void ElementEvaluationAll_phi(const int nen, const vector<float> &Nx, vector<array<float, 3>> &dNdx,
	// 							const vector<float> elePhiGuess, float &elePG,
	// 							const vector<float> elePhi, float &eleP,
	// 							const vector<float> eleSyn, float &eleS,
	// 							const vector<float> eleTips, float &eleTp,
	// 							const vector<float> eleTubulin, float &eleTb,
	// 							float &dPGdx, float &dPGdy, float &dPGdz);
	// void ElementEvaluationAll_phi(const int nen, const vector<float> &Nx, vector<array<float, 3>> &dNdx,
	// 							const vector<float> elePhiGuess, float &elePG,
	// 							const vector<float> elePhi, float &eleP,
	// 							const vector<float> eleSyn, float &eleS,
	// 							const vector<float> eleTips, float &eleTp,
	// 							const vector<float> eleTubulin, float &eleTb,
	// 							const vector<float> eleEpsilon, float &eleEP,
	// 							const vector<float> eleEpsilonP, float &eleEEP,
	// 							float &dPGdx, float &dPGdy,
	// 							float &dAdx, float &dAdy,
	// 							float &dAPdx, float &dAPdy);
	inline void ElementEvaluationAll_phi(
		int nen,
		const vector<float> &Nx,
		const vector<array<float, 3>> &dNdx,
		const vector<vector<float>> &eleVal,
		vector<float> &vars);

	void ElementEvaluationAll_phi_test(const uint &nen, 
								const vector<float> &Nx, 
								const vector<array<float, 3>> &dNdx, 
								const vector<float> &elePhiGuess, 
								vector<float> &vars);

	void ElementEvaluationAll_phi_opt(const uint &nen, 
								const vector<float> &Nx, 
								const vector<array<float, 3>> &dNdx, 
								const vector<float> &elePhiGuess, 
								vector<float> &vars);

	// void ElementEvaluationAll_phi(const int nen, const vector<float> &Nx, vector<array<float, 3>> &dNdx,
	// 	vector<vector<float>> &eleVal, vector<float> &vars);
	void ElementEvaluationAll_syn_tub(const int nen, const vector<float> &Nx, vector<array<float, 3>> &dNdx,
		const vector<float> elePhiDiff, float &elePf, 
		const vector<float> eleSyn,float &eleS, const vector<float> elePhi, float &eleP,
		const vector<float> elePhiPrev, float &elePprev,
		const vector<float> eleConct, float &eleC,
		float &dPdx, float &dPdy, float &dPdz);
	void ElementEvaluationAll_syn_tub(const int nen, const vector<float> &Nx, vector<array<float, 3>> &dNdx,
		vector<vector<float>> &eleVal, vector<float> &vars);

	// Pre-computation to reduce redundant calculations
	void PrepareBasis();              // Precompute basis functions and derivatives
	void PreparePhaseField();         // Precompute variables specific to the phase field equation
	void PrepareTermSource();         // Precompute source term contributions
	void PrepareEpsilon();            // Precompute epsilon and related values

	// Phase Field Equation Evaluations
	void EvaluateEnergy(
		const int nen,                // Number of nodes in the element
		const vector<float>& Nx,      // Basis function values
		const vector<float>& eleSyn,  // Synaptic field values
		vector<float>& E              // Output energy values
	);
	float Regular_Heiviside_fun(float x); // Smooth Heaviside function implementation

	// Orientation Evaluations
	void EvaluateOrientation_prev(
		const uint& nen,                  		// Number of nodes in the element
		const vector<float>& Nx,          		// Basis function values
		const vector<array<float, 3>>& dNdx, 	// Basis function derivatives
		const vector<float>& elePhi,      		// Phase field values
		const vector<float>& eleTheta,    		// Orientation angle
		vector<float>& eleEpsilon,        		// Output epsilon values
		vector<float>& eleEpsilonP        		// Output epsilon derivative values
	);

	void EvaluateOrientation(
		const int nen,                    		// Number of nodes in the element
		const vector<float>& Nx,          		// Basis function values
		const vector<array<float, 3>>& dNdx, 	// Basis function derivatives
		const vector<float>& elePhi,      		// Phase field values
		const vector<float>& eleTheta,    		// Orientation angle
		float& eleAniso,                  		// Anisotropy factor
		float& dA_dPdx,                   		// Derivative of anisotropy w.r.t x
		float& dA_dPdy,                   		// Derivative of anisotropy w.r.t y
		float& dA_dPdz                   		// Derivative of anisotropy w.r.t z
	);

	void EvaluateOrientationSpherical(
		const int nen,                    		// Number of nodes in the element
		const vector<float>& Nx,          		// Basis function values
		const vector<array<float, 3>>& dNdx, 	// Basis function derivatives
		const vector<float>& elePhi,      		// Phase field values
		const vector<float>& elePolar,    		// Polar angle values
		const vector<float>& eleAzimuth,  		// Azimuthal angle values
		float& eleEpsilon,                		// Output epsilon value
		float dEdp,                       		// Derivative w.r.t polar angle
		float dEda                        		// Derivative w.r.t azimuthal angle
	);

	// Linear System Assembly for the Phase Field Equation
	void BuildLinearSystemProcessNG_phi(
		const vector<Vertex3D>& cpts      // Control points for Bezier mesh
	);

	// Synaptogenesis and Tubulin Operations
	void CalculateSumGradPhi0(const vector<Vertex3D>& cpts); 
		// Calculate the sum of gradient phi0 over the control points

	void BuildLinearSystemProcessNG_syn_tub(const vector<Vertex3D>& cpts); 
		// Assemble the linear system for synaptogenesis and tubulin equations

	// Domain Expansion Operations
	int CheckExpansion3D(
		const vector<float>& input,               							// Input variable for checking expansion
		const vector<Vertex3D>& cpts,      							// Control points for the domain
		const int& NX, const int& NY, const int& NZ,            	// Domain dimensions in X, Y, Z
		const int& originX, const int& originY, const int& originZ 	// Origin coordinates for the domain
	); // Check and handle 3D domain expansion conditions

	void PopulateRandom(vector<float>& input);
		// Populate input vector with random values after domain expansion

	// Interpolate Variables Using Coarse KDTree
	vector<float> InterpolateVars_coarseKDtree(
		const vector<float>& input,           	// Input values from the initial control points
		const vector<Vertex3D>& cpts_initial, 	// Initial control points for interpolation
		const KDTree& kdTree_initial,           // KDTree built from the initial control points
		const vector<Vertex3D>& cpts,         	// Current control points for interpolation
		int type,                               // Interpolation type: 0 (max), 1 (average), 2 (zero)
		int isTheta                             // Special flag for theta handling: 1 for special out-of-bound handling
	); // Interpolates values from the coarse KDTree to the current control points.

	// KDTree Search for Pair Matching
	bool KD_SearchPair(
		const vector<Vertex3D>& cpts,      				// Control points for searching
		const KDTree& kdTree,              				// KDTree for spatial indexing
		float targetX, float targetY, float targetZ, 	// Target coordinates to search
		int& ind,                          				// Output: index of the closest match
		float tolerance = 1.0f             				// Tolerance for pair matching
	); // Search for a pair in KDTree within a tolerance
	
	// Tip Detection Functions
	float RmOutlier(vector<float>& data);
		// Remove outliers from the dataset using a standard deviation-based method

	float CellBoundary(float phi, float threshold);
		// Determine cell boundary based on a given threshold value

	void DetectTipsMulti3D(
		vector<float> id,                // Input neuron IDs
		int numNeuron,                   // Total number of neurons
		vector<float>& tips,             // Output vector to store detected tips
		int NX, int NY, int NZ           // Dimensions of the 3D domain
	); // Detect tips in a 3D multi-neuron setup

	// Spatial Operations
	bool IsInBox(
		const Vertex3D& point,           // Point to check
		const Vertex3D& center,          // Center of the box
		float dx, float dy, float dz     // Half-dimensions of the box
	); // Check if a given point is within a specified 3D box centered at 'center'

	// Sum Calculation for Phi within a Specified Box
	void CalculatePhiSum(
		const vector<Vertex3D>& cpts,    // Control points representing center points
		float dx, float dy, float dz,    // Half-dimensions of the box
		const KDTree& kdTree             // KDTree for spatial indexing
	); // Calculate the sum of phi values within a 3D box for each center point in `cpts`

	vector<pair<Vertex3D, int>> FindClosestVerticesWithIndices(
		const vector<Vertex3D>& vertices,     // List of vertices
		const Vertex3D& inputVertex,          // Target vertex
		int k                                 // Number of closest vertices to find
	); // Find the k closest vertices to a given vertex along with their indices
		
	// vector<tuple<Vertex3D, int, float>> FindClosestVerticesWithIndicesAndDistances(const vector<Vertex3D>& vertices, const Vertex3D& inputVertex, int k);
	
	// 3D Vertex and Cluster Operations
	vector<tuple<Vertex3D, int, float>> FindClosestVerticesWithIndicesAndDistances(
		const KDTree& kdTree, const Vertex3DCloud& cloud, const Vertex3D& inputVertex, int k
	); // Find the k closest vertices along with indices and distances using KDTree.

	// Refinement and Clustering
	vector<float> ComputeRefine(const vector<float>& phi_in,
		int NX, int NY, int NZ,
		int &originX, int &originY, int &originZ,
		const KDTree& kdTree, const Vertex3DCloud& cloud);

	void BFS3D(
		const vector<float>& matrix, int depth, int rows, int cols, 
		int dep, int row, int col, vector<bool>& visited, 
		vector<tuple<int, int, int>>& cluster
	); // Perform Breadth-First Search in a 3D matrix to find connected clusters.

	vector<vector<tuple<int, int, int>>> FindClusters3D(
		const vector<float>& matrix, int depth, int rows, int cols
	); // Identify clusters of connected points in a 3D matrix.

	vector<float> FindLocalMaximaInClusters3D(
		const vector<float>& matrix, int depth, int rows, int cols
	); // Locate local maxima in identified 3D clusters.

	// Neuron Detection and Processing
	vector<vector<vector<int>>> ConvertTo3DIntVector(
		const vector<float>& input, int NX, int NY, int NZ
	); // Convert a 1D float vector to a 3D integer vector.

	vector<vector<vector<float>>> ConvertTo3DFloatVector(
		const vector<float>& input, int NX, int NY, int NZ
	); // Convert a 1D float vector to a 3D float vector.

	void FloodFill3DWithKDTree(
		vector<vector<vector<int>>>& image, int x, int y, int z, 
		int newColor, int originalColor, const KDTree& kdTree, 
		const Vertex3DCloud& cloud
	); // Perform 3D flood fill with KDTree for spatial connectivity.

	void IdentifyNeurons3DWithKDTree(
		vector<vector<vector<int>>>& neurons, const vector<array<int, 3>>& seed,
		int NX, int NY, int NZ, int originX, int originY, int originZ,
		const KDTree& kdTree, const Vertex3DCloud& cloud
	); // Identify neurons in a 3D grid using KDTree and seed points.

	bool IsValid(const int& x, const int& y, const int& z, 
				const int& rows, const int& cols, 
				const int& depth
	); // Check if a 3D point is valid within specified bounds.

	vector<vector<vector<int>>> CalculateGeodesicDistanceFromPoint3D(
		vector<vector<vector<int>>> neurons, const vector<array<int, 3>>& seed, 
		int originX, int originY, int originZ
	); // Calculate geodesic distances from a given seed point in 3D.

	// Save and Output Operations
	void SaveNGvars(
		const vector<vector<float>>& NGvars, int NX, int NY, const string& fn
	); // Save NGvars to a file with specified dimensions.

	void PrintOutNeurons3D(
		const vector<vector<vector<int>>>& neurons
	); // Print the neuron 3D structure.

	void PrintStatus(
		int n, int end_iter, int reason_phi, int its_phi, double t_phi,
		int reason_syn, int its_syn, double t_syn, 
		int reason_tub, int its_tub, double t_tub, int n_bzmesh
	); // Print the current status of the simulation with aligned output.
};

// Phase Field PETSc Nonlinear SNES Solver Functions
PetscErrorCode SetupSNES(
    SNES &snes, const char *solverType, void *ctx,
    PetscErrorCode (*formFunction)(SNES, Vec, Vec, void *),
    PetscErrorCode (*formJacobian)(SNES, Vec, Mat, Mat, void *),
    PetscReal rtol, PetscReal atol, PetscReal stol, PetscInt maxIters, PetscInt maxFails
); // Configures and initializes a PETSc SNES solver for nonlinear systems.
// PetscErrorCode SetupSNES(
//     SNES &snes, 
//     const char *solverType, 
//     void *ctx,
//     PetscErrorCode (*formFunction)(SNES, Vec, Vec, void *),
//     PetscErrorCode (*formJacobian)(SNES, Vec, Mat, Mat, void *),
//     PetscReal rtol = 1e-5, 
//     PetscReal atol = 1e-7, 
//     PetscReal dtol = 1e-9,
//     PetscInt maxIters = 100, 
//     PetscInt maxFails = 1000, 
//     PetscReal lineSearchDamping = 0.8, 
//     SNESLineSearchType lineSearchType = SNESLINESEARCHCP
// );

PetscErrorCode SetupKSP(
    KSP &ksp, Mat &A, const char *kspType, const char *pcType,
    PetscReal rtol, PetscReal atol, PetscReal dtol, PetscInt maxIters, PetscInt restart
); // Configures and initializes a PETSc KSP solver for linear systems.

PetscErrorCode ScatterVector(
    Vec src, vector<float>& target, PetscInt size, bool applyBoundary, NeuronGrowth* NG
); // Scatters a PETSc vector into a local float vector, optionally applying boundary conditions.

PetscErrorCode FormFunction_phi(
    SNES snes, Vec x, Vec F, void *ctx
); // Defines the nonlinear residual function for the phase field equation.

PetscErrorCode FormJacobian_phi(
    SNES snes, Vec x, Mat J, Mat P, void *ctx
); // Defines the Jacobian matrix for the phase field equation.

PetscErrorCode MySNESMonitor(
    SNES snes, PetscInt its, PetscReal fnorm, PetscViewerAndFormat *vf
); // Custom monitor for SNES solver to track progress during iterations.

PetscErrorCode CleanUpSolvers(NeuronGrowth &NG);
	// Cleans up solver-related memory allocations for NeuronGrowth object.

// Main Simulation Driver
int RunNG(
    const int n_bzmesh, vector<vector<int>> ele_process_in,
    vector<Vertex3D> &cpts_initial, vector<Vertex3D> &cpts, vector<Vertex3D> prev_cpts,
    string path_in, string path_out,
    int &iter, int end_iter_in,
    vector<vector<float>> &NGvars,
    int &NX, int &NY, int &NZ,
    vector<array<float, 3>> &seed,
	int &originX, int &originY, int &originZ,
    bool &localRefine,
	const string& phi_solver
); // Runs the Neuron Growth simulation for the specified input parameters.

#endif