#include "NeuronGrowth.h"               // Main header for NeuronGrowth class
#include "utils.h"                      // Utility functions
#include <cmath>                        // For mathematical operations
#include <fstream>                      // For file operations
#include <iostream>                     // For console I/O
#include <sstream>
#include <iomanip>                      // For formatted I/O
#include <algorithm>                    // For STL algorithms
#include <numeric>                      // For numeric operations
#include <stack>                        // For tic-toc timing
#include <ctime>                        // For srand
#include <queue>                        // For geodesic distance
#include <cfloat>                       // For FLT_MAX
#include <limits>                       // For numeric limits
#include <tuple>                        // For handling tuples
#include <chrono>                       // For modern timing
#include "../nanoflann/1.5.5/include/nanoflann.hpp" // For KD-tree implementation
#include <petsctime.h> // PETSc timing functions

using namespace std;

// Type Definitions
typedef unsigned int uint;

// Constants
constexpr float PI = 3.1415926f;        // Pi constant
constexpr int INF = 1e9;               // Infinity constant for large integers

stack<double> tictoc_stack; // Use double for higher precision
double t_phi(0), t_syn_tub(0), t_collect(0), t_write(0);

// Start timing
void tic() {
    PetscLogDouble current_time;
    PetscTime(&current_time); // Get the current time
    tictoc_stack.push(current_time);
}

// End timing and calculate elapsed time
void toc(double &elapsed_time) {
    if (!tictoc_stack.empty()) {
        PetscLogDouble current_time;
        PetscTime(&current_time); // Get the current time
        elapsed_time = current_time - tictoc_stack.top(); // Calculate elapsed time
        tictoc_stack.pop();
    } else {
        elapsed_time = 0.0; // Handle case where tic was not called
    }
}

void UpdateSimulationTimers(double& t_check, double& t_write, double& t_global) {
    toc(t_check);       // Measure total time for the current solver/process
    tic();            	// Start timing for the next process
    t_write += t_check; // Accumulate write time
    t_global += t_check; // Accumulate total simulation time
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

inline float SquaredDistance(const Vertex3D& first, const Vertex3D& other) {
    const float dx = first.coor[0] - other.coor[0];
    const float dy = first.coor[1] - other.coor[1];
    const float dz = first.coor[2] - other.coor[2];
    return dx * dx + dy * dy + dz * dz;
}

void CheckAndPrintThresholdExceedance(const vector<float>& input, float threshold) {
    bool found = false;  // Flag to track if any value exceeds the threshold

    for (size_t i = 0; i < input.size(); ++i) {
        if (input[i] > threshold) {
            cout << "Value " << input[i] << " at index " << i << " exceeds the threshold " << threshold << ".\n";
            found = true;
        }
    }

    if (!found) {
        cout << "No values exceed the threshold " << threshold << ".\n";
    }
}

NeuronGrowth::NeuronGrowth(const string& phi_solver,
						const int iter,
						const int numNeuron_in,
						const int end_iter) 
{
    // MPI setup
    comm = PETSC_COMM_WORLD; // Initialize MPI communicator
    mpiErr = MPI_Comm_rank(comm, &comRank); // Get rank of the process
    mpiErr = MPI_Comm_size(comm, &comSize); // Get size of the communicator
    nProcess = comSize; // Total number of processes

    // Flags and initialization
    judge_phi = 0; // Phase field judgment flag
    judge_syn = 0; // Synapse judgment flag
    judge_tub = 0; // Tubule judgment flag
    sum_grad_phi0_local = 0; // Local gradient sum of phi
    sum_grad_phi0_global = 0; // Global gradient sum of phi

	this->phi_solver = phi_solver;
	this->n = iter;
	this->numNeuron = numNeuron_in;
	this->end_iter = end_iter;

	if (phi_solver == "ksp") {
		// Simulation parameters
		var_save_invl   = 100;       // Interval for saving variables
		expandCK_invl   = 5000;     // Interval for expanding control knots
		expand_sz 		= 3.0f;		// Domain expansion size
		tip_detect_invl = 250;		// interval for tip detection

		// Neuron-specific parameters
		aniso           = 6;         // Anisotropy constant
		numNeuron       = 1;         // Number of neurons
		gc_sz           = 2;         // Growth cone size

		// Diffusion and mobility
		Diff            = 4;         // Diffusion coefficient
		M_axon          = 100;       // Mobility for axon
		M_neurite       = 50;        // Mobility for neurite
		// M_phi           = 60;        // Mobility for phase field

		// Growth-related parameters
		alpha           = 0.9;       // Growth rate scaling factor
		alphaT          = 0.001;     // Tubulin production scaling factor
		betaT           = 0.001;     // Tubulin degradation scaling factor
		c_opt           = 1;         // Optimization constant for growth
		delta           = 0.20;      // Growth anisotropy coefficient

		// Time-stepping
		dt              = 0.0005;    // Time step size for no-local refinement
		// dt              = 0.0001;    // Time step size for local refinement

		// Phase field model parameters
		epsilonb        = 0.04;      // Baseline epsilon for growth anisotropy
		g               = 0.1;       // Threshold for growth regulation
		gamma           = 10;        // Growth feedback factor
		k2              = 0;         // Secondary feedback term
		kappa           = 2;         // Stiffness constant for bending
		Dc              = 6;         // Diffusion coefficient for concentration

		// Boundary and seed-related parameters
		kp75            = 0;         // Placeholder constant for tuning
		r               = 5;         // Seed radius
		s_coeff         = 0.007;     // Source coefficient for growth
		
		seed_radius     = 5;         // Initial seed radius
		source_coeff    = 15;        // Source term coefficient

	} else if (phi_solver == "snes") {
		// Simulation settings
		expandCK_invl   = 5000;    	// Interval for expanding control knots
		var_save_invl   = 100;      // Interval for saving variables
		refine_invl		= 5000;    	// Interval for local refinement
		tip_detect_invl = 250;		// interval for tip detection
		tip_threshold	= 0.0;		// tip thresold in formfunction
		expand_sz 		= 3.0f;		// Domain expansion size

		numNeuron       = 1;        // Number of neurons
		seed_radius     = 4;        // Initial seed radius for neuron growth
		dt              = 5e-3;     // Time step size

		// Phase field growth-related parameters
		aniso           = 6;        // Anisotropy constant
		gamma           = 10;       // Growth feedback factor
		kappa           = 1.8;      // Stiffness constant for bending
		Dc              = 3;        // Diffusion coefficient for synaptic concentration
		alpha           = 0.9;      // Growth rate scaling factor
		alphaOverPi     = alpha / PI; // Alpha normalized over π
		M_phi           = 1;       // Mobility for phase field
		M_axon          = 10;       // Mobility for phase field
		M_neurite       = 5;       // Mobility for phase field
		s_coeff         = 0.007;    // Source coefficient for growth
		delta           = 0.50;     // Growth anisotropy coefficient
		epsilonb        = 0.01;     // Baseline epsilon for anisotropy

		// Tubulin & synaptogensis parameters
		r               = 5;        // Radial growth parameter
		g               = 0.1;      // Growth sensitivity factor
		alphaT          = 0.001;    // Tubulin production scaling factor
		betaT           = 0.001;    // Tubulin degradation scaling factor
		Diff            = 4;        // Diffusion coefficient
		source_coeff    = 15;       // Source term coefficient
	}
}

// Assign processor-specific elements to the local process
void NeuronGrowth::AssignProcessor(vector<vector<int>> &ele_proc) {
    ele_process = ele_proc[comRank]; // Direct assignment is more efficient
}

// Load variables from a parameter file
void NeuronGrowth::SetVariables(const string &fn_par) {
    ifstream inputFile(fn_par);
    if (!inputFile.is_open()) {
        PetscPrintf(PETSC_COMM_WORLD, "Error: Unable to open file %s\n", fn_par.c_str());
        return;
    }

    string line;
    while (getline(inputFile, line)) {
        istringstream iss(line);
        string variableName;
        char equalsSign;
        float value;

        // Parse each line for a variable and its value
        if (!(iss >> variableName >> equalsSign >> value) || equalsSign != '=') {
            PetscPrintf(PETSC_COMM_WORLD, "Warning: Invalid format in line: %s\n", line.c_str());
            continue;
        }

        // Map variable names to class members
        if (variableName == "expandCK_invl") expandCK_invl = value;
        else if (variableName == "refine_invl") refine_invl = value;
        else if (variableName == "var_save_invl") var_save_invl = value;
        else if (variableName == "tip_detect_invl") tip_detect_invl = value;
        else if (variableName == "tip_threshold") tip_threshold = value;
        else if (variableName == "expand_sz") expand_sz = value;
        else if (variableName == "numNeuron") numNeuron = value;
        else if (variableName == "gc_sz") gc_sz = value;
        else if (variableName == "aniso") aniso = value;
        else if (variableName == "gamma") gamma = value;
        else if (variableName == "seed_radius") seed_radius = value;
        else if (variableName == "kappa") kappa = value;
        else if (variableName == "dt") dt = value;
        else if (variableName == "Dc") Dc = value;
        else if (variableName == "kp75") kp75 = value;
        else if (variableName == "k2") k2 = value;
        else if (variableName == "c_opt") c_opt = value;
        else if (variableName == "alpha") alpha = value;
        else if (variableName == "M_phi") M_phi = value;
        else if (variableName == "M_axon") M_axon = value;
        else if (variableName == "M_neurite") M_neurite = value;
        else if (variableName == "s_coeff") s_coeff = value;
        else if (variableName == "delta") delta = value;
        else if (variableName == "epsilonb") epsilonb = value;
        else if (variableName == "r") r = value;
        else if (variableName == "g") g = value;
        else if (variableName == "alphaT") alphaT = value;
        else if (variableName == "betaT") betaT = value;
        else if (variableName == "Diff") Diff = value;
        else if (variableName == "source_coeff") source_coeff = value;
        else {
            PetscPrintf(PETSC_COMM_WORLD, "Warning: Unknown variable %s in file.\n", variableName.c_str());
        }
    }

    inputFile.close();
    PetscPrintf(PETSC_COMM_WORLD, "Parameters loaded successfully.\n");

    // Recompute dependent parameters
    alphaOverPi = alpha / PI;
}

void NeuronGrowth::InitializeProblemNG(const int n_bz,
									vector<Vertex3D>& cpts,
									const Vertex3DCloud& cloud,
									KDTree& kdTree,
									const vector<Vertex3D>& prev_cpts,
									const Vertex3DCloud& cloud_prev,
									KDTree& kdTree_prev,
									vector<vector<float>>& NGvars,
									vector<array<float, 3>>& seed)
{
	// Synchronize all processes and start initialization
	MPI_Barrier(PETSC_COMM_WORLD);
	PetscPrintf(PETSC_COMM_WORLD, "Initializing Neuron Growth Problem--------------------------------------------\n");

	// Initialize parameters
	int cpt_sz = cpts.size(); // Number of control points
	GaussInfo(3);
	n_bzmesh = n_bz;          // Number of Bezier elements
	N_0.resize(cpt_sz);       // Resize N_0 to match control points

	// Determine bounds for current control points
	// float max_x, min_x, max_y, min_y, max_z, min_z;
	for (int i = 0; i < cpts.size(); ++i) {
		const auto& coord = cpts[i].coor;
		if (i == 0) {
			// Initialize bounds with the first control point
			max_x = min_x = coord[0];
			max_y = min_y = coord[1];
			max_z = min_z = coord[2];
		} else {
			// Update bounds for subsequent control points
			max_x = max(coord[0], max_x); min_x = min(coord[0], min_x);
			max_y = max(coord[1], max_y); min_y = min(coord[1], min_y);
			max_z = max(coord[2], max_z); min_z = min(coord[2], min_z);
		}
	}

	// Determine bounds for previous control points
	float max_x_prev, min_x_prev;
	float max_y_prev, min_y_prev;
	float max_z_prev, min_z_prev;
	for (int i = 0; i < prev_cpts.size(); ++i) {
		const auto& coord = prev_cpts[i].coor;
		if (i == 0) {
			// Initialize bounds with the first previous control point
			max_x_prev = min_x_prev = coord[0];
			max_y_prev = min_y_prev = coord[1];
			max_z_prev = min_z_prev = coord[2];
		} else {
			// Update bounds for subsequent previous control points
			max_x_prev = max(coord[0], max_x_prev); min_x_prev = min(coord[0], min_x_prev);
			max_y_prev = max(coord[1], max_y_prev); min_y_prev = min(coord[1], min_y_prev);
			max_z_prev = max(coord[2], max_z_prev); min_z_prev = min(coord[2], min_z_prev);
		}
	}

	// Log the calculated bounds for debugging
	PetscPrintf(PETSC_COMM_WORLD, "Bounds for current control points:\n");
	PetscPrintf(PETSC_COMM_WORLD, "max x: %.2f, min x: %.2f | max y: %.2f, min y: %.2f | max z: %.2f, min z: %.2f\n", 
				max_x, min_x, max_y, min_y, max_z, min_z);
	PetscPrintf(PETSC_COMM_WORLD, "Bounds for previous control points:\n");
	PetscPrintf(PETSC_COMM_WORLD, "max px: %.2f, min px: %.2f | max py: %.2f, min py: %.2f | max pz: %.2f, min pz: %.2f\n", 
				max_x_prev, min_x_prev, max_y_prev, min_y_prev, max_z_prev, min_z_prev);

	// Initialize parameters for the simulation
	float r, x, y, z;  // Radius and spatial coordinates
	srand(static_cast<unsigned int>(time(NULL)));  // Seed for random number generator

	if (n == 0) {  // Initialize simulation variables for the first run
		// Initialize control point variables
		for (int i = 0; i < cpts.size(); ++i) {
			phi.push_back(0.0f);
			tub.push_back(0.0f);
			theta.push_back(static_cast<float>(rand() % 100) / 100.0f);  // Random orientation (theta)
			syn.push_back(0.0f);  // Synaptogenesis variable
		}

		// Assign values and boundary labels for each control point
		for (int i = 0; i < cpts.size(); ++i) {
			x = cpts[i].coor[0];
			y = cpts[i].coor[1];
			z = cpts[i].coor[2];

			// Assign boundary labels based on spatial bounds
			cpts[i].label = (x == min_x || x == max_x || y == min_y || y == max_y || z == min_z || z == max_z) ? 1 : 0;

			// Check for initial soma placement around seeds
			for (const auto& seed_point : seed) {
				r = sqrt(pow(x - seed_point[0], 2) + pow(y - seed_point[1], 2) + pow(z - seed_point[2], 2));
				if (r <= seed_radius) {  // Inside seed radius
					phi[i] = 1.0f;  // Initial phi
					tub[i] = 0.5f + 0.5f * tanh((sqrt(seed_radius) - r) / 2.0f);  // Based on literature equation
					// cpts[i].label = 1;
				}
			}

			// vector<tuple<Vertex3D, int, float>> closestVertices = FindClosestVerticesWithIndicesAndDistances(kdTree, cloud, cpts[i], 4);
		}

		// Save initial phi and tub states
		phi_0 = phi;
		tub_0 = tub;

	} else {  // Continue simulation using NGvars
		// Initialize variables
		phi.assign(cpt_sz, 0.0f);
		syn.assign(cpt_sz, 0.0f);
		tub.assign(cpt_sz, 0.0f);
		theta.assign(cpt_sz, 0.0f);
		phi_0.assign(cpt_sz, 0.0f);
		tub_0.assign(cpt_sz, 0.0f);

// #pragma omp parallel for
		PetscPrintf(PETSC_COMM_WORLD, "Vector Sizes: phi: %zu, syn: %zu, tub: %zu, theta: %zu, phi_0: %zu, tub_0: %zu\n", phi.size(), syn.size(), tub.size(), theta.size(), phi_0.size(), tub_0.size());
		PetscPrintf(PETSC_COMM_WORLD, "cpt Sizes: phi: %d\n", cpt_sz);
		for (size_t i = 0; i < cpt_sz; ++i) {
			const auto& cpt = cpts[i];
			const auto& [x, y, z] = cpt.coor;

			// Assign boundary labels based on spatial bounds
			cpts[i].label = (x == min_x || x == max_x ||
								y == min_y || y == max_y ||
								z == min_z || z == max_z) ? 1 : 0;

			// // Check for initial soma placement around seeds
			// for (const auto& seed_point : seed) {
			// 	r = sqrt(pow(x - seed_point[0], 2) + pow(y - seed_point[1], 2) + pow(z - seed_point[2], 2));
			// 	if (r <= seed_radius) {  // Inside seed radius
			// 		cpts[i].label = 1;
			// 	}
			// }

			// Check if the point is within bounds
			bool withinBounds = (x > min_x_prev && x < max_x_prev &&
									y > min_y_prev && y < max_y_prev &&
									z > min_z_prev && z < max_z_prev);

			if (withinBounds) {
				
				// Interpolate or find exact match for current point
				InterpolateOrFindExact(
					cpt, kdTree_prev, cloud_prev, NGvars, prev_cpts, 
					phi[i], syn[i], tub[i], theta[i], phi_0[i], tub_0[i]);	
			} else {
				// Out-of-bounds handling (defaults)
				phi[i] = syn[i] = tub[i] = phi_0[i] = tub_0[i] = 0.0f;
				theta[i] = static_cast<float>(rand() % 100) / 100.0f;
			}
		}
	}
	// Assign control points for further processing
	this->cpts = cpts;

	/* Initialize PETSc vectors and matrices */
	PetscInt mat_dim = cpt_sz; // Dimension of PETSc matrices and vectors

	// Temporary vector for phi
	ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, mat_dim, &temp_phi);

	// Initialize matrix J (Jacobian matrix or general-purpose matrix)
	// 250 = 2 * 5^3 (3D: 2 variables, 5 basis functions) || 25 = 1 * 5^2 (2D: 1 variable, 5 basis functions)
	ierr = MatCreate(PETSC_COMM_WORLD, &J);
	ierr = MatSetSizes(J, PETSC_DECIDE, PETSC_DECIDE, mat_dim, mat_dim);
	ierr = MatSetType(J, MATMPIAIJ); // Matrix type: MPIAIJ (sparse matrix)
	ierr = MatMPIAIJSetPreallocation(J, 250, NULL, 250, NULL); // Preallocate non-zero entries
	ierr = MatSetOption(J, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE); // Allow dynamic allocation
	ierr = MatSetUp(J);

	// Initialize matrices and vectors for synaptogenesis (GK_phi and GR_phi)
	ierr = MatCreate(PETSC_COMM_WORLD, &GK_phi);
	ierr = MatSetSizes(GK_phi, PETSC_DECIDE, PETSC_DECIDE, mat_dim, mat_dim);
	ierr = MatSetType(GK_phi, MATMPIAIJ);
	ierr = MatMPIAIJSetPreallocation(GK_phi, 250, NULL, 250, NULL);
	ierr = MatSetOption(GK_phi, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
	ierr = MatSetUp(GK_phi);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, mat_dim, &GR_phi);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, mat_dim, &temp_phi);

	// Initialize matrices and vectors for synaptogenesis (GK_syn and GR_syn)
	ierr = MatCreate(PETSC_COMM_WORLD, &GK_syn);
	ierr = MatSetSizes(GK_syn, PETSC_DECIDE, PETSC_DECIDE, mat_dim, mat_dim);
	ierr = MatSetType(GK_syn, MATMPIAIJ);
	ierr = MatMPIAIJSetPreallocation(GK_syn, 250, NULL, 250, NULL);
	ierr = MatSetOption(GK_syn, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
	ierr = MatSetUp(GK_syn);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, mat_dim, &GR_syn);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, mat_dim, &temp_syn);

	// Initialize matrices and vectors for tubule dynamics (GK_tub and GR_tub)
	ierr = MatCreate(PETSC_COMM_WORLD, &GK_tub);
	ierr = MatSetSizes(GK_tub, PETSC_DECIDE, PETSC_DECIDE, mat_dim, mat_dim);
	ierr = MatSetType(GK_tub, MATMPIAIJ);
	ierr = MatMPIAIJSetPreallocation(GK_tub, 250, NULL, 250, NULL);
	ierr = MatSetOption(GK_tub, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
	ierr = MatSetUp(GK_tub);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, mat_dim, &GR_tub);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, mat_dim, &temp_tub);

	// Reserve memory for additional variables and clear pre-solve containers
	vars.reserve(21 * sizeof(float)); // Preallocate space for 21 variables
	vars.resize(21);                  // Resize to hold 21 entries

	pre_EMatrixSolve.clear(); // Clear pre-solve matrix container
	pre_EVectorSolve.clear(); // Clear pre-solve vector container
    pre_eleVal.clear();      // Element values

	int maxNen = 200; // estimated maximum nen in locally refined mesh - obtained from testing
    EMatrixSolve.resize(maxNen, vector<float>(maxNen, 0.0f));
    EVectorSolve.resize(maxNen, 0.0f);
    eleVal.resize(8, vector<float>(maxNen, 0.0f)); // 8 fields, preallocated
}

// Function to interpolate or find exact phi values for a given point
void NeuronGrowth::InterpolateOrFindExact(
    const Vertex3D& cpt, 
    const KDTree& kdTree_prev, 
    const Vertex3DCloud& cloud_prev, 
    const vector<vector<float>>& NGvars, 
    const vector<Vertex3D>& prev_cpts, 
    float& phi, float& syn, float& tub, float& theta, float& phi_0, float& tub_0
) {
    phi = syn = tub = theta = phi_0 = tub_0 = 0.0f;

    int exactIndex;
    if (KD_SearchPair(prev_cpts, kdTree_prev, cpt.coor[0], cpt.coor[1], cpt.coor[2], exactIndex)) {
        // Exact match: copy values directly
        phi   = NGvars[0][exactIndex];
        syn   = NGvars[1][exactIndex];
        tub   = NGvars[2][exactIndex];
        theta = NGvars[3][exactIndex];
        phi_0 = NGvars[4][exactIndex];
        tub_0 = NGvars[5][exactIndex];
		return;
    } else {
        // Interpolate using nearest neighbors
        auto closestVertices = FindClosestVerticesWithIndicesAndDistances(kdTree_prev, cloud_prev, cpt, 4);

        float totalWeight = 0.0f;
        for (const auto& [_, idx, distance] : closestVertices) {
            float weight = 1.0f / (distance + 1e-6f); // Inverse-distance weight
            phi   += NGvars[0][idx] * weight;
            syn   += NGvars[1][idx] * weight;
            tub   += NGvars[2][idx] * weight;
            // theta += NGvars[3][idx] * weight;
            phi_0 += NGvars[4][idx] * weight;
            tub_0 += NGvars[5][idx] * weight;
            totalWeight += weight;
        }

        // Normalize interpolated values
        phi   /= totalWeight;
        syn   /= totalWeight;
        tub   /= totalWeight;
        // theta /= totalWeight;
		theta = static_cast<float>(rand() % 100) / 100.0f;
        phi_0 /= totalWeight;
        tub_0 /= totalWeight;
    }
}

void NeuronGrowth::InterpolateOrFindExact_singleVar(
    const Vertex3D& cpt, 
    const KDTree& kdTree_prev, 
    const Vertex3DCloud& cloud_prev, 
    const vector<float>& original_var, 
    const vector<Vertex3D>& prev_cpts,
    float& output_var,
	bool weighted
) {
    int exactIndex(0);

    if (KD_SearchPair(prev_cpts, kdTree_prev, cpt.coor[0], cpt.coor[1], cpt.coor[2], exactIndex)) {
        // Exact match: copy values directly
        output_var = original_var[exactIndex];
    } else {
        // Interpolate using nearest neighbors
        auto closestVertices = FindClosestVerticesWithIndicesAndDistances(kdTree_prev, cloud_prev, cpt, 6);

		if (closestVertices.empty()) {
			cerr << "KD-Tree validation failed: No neighbors found for sample point." << endl;
		} 
		// else {
			// cout << "KD-Tree validation successful: Found " << closestVertices.size() << " neighbors for sample point." << endl;
		// }

		float totalWeight = 0.0f;
		output_var = 0.0f; // Initialize before accumulation

		for (const auto& [_, idx, distance] : closestVertices) {
			float weight;
			if (weight) {
				weight = 1.0f / (distance + 1e-6f); // Inverse-distance weight
			} else {
				weight = 1.0f;
			}
			output_var += original_var[idx] * weight;
			totalWeight += weight;
		}

		// Normalize interpolated values
		if (totalWeight > 0.0f) {
			output_var /= totalWeight;
		} else {
			output_var = 0.0f; // Handle zero totalWeight gracefully
		}
    }
}

void NeuronGrowth::CheckVar(const string& fn, const vector<Vertex3D>& cpts, const vector<float>& input) 
{
	// Construct the output file name
	// string fname = fn + "_" + to_string(n) + ".vtk";
	string fname = path_out + "/" + fn + "_" + to_string(n) + ".vtk";
	
	// Open the output file
	ofstream fout(fname.c_str());
	if (fout.is_open()) {
		// Write VTK file header
		fout << "# vtk DataFile Version 2.0\n";
		fout << "Square plate test\n";
		fout << "ASCII\n";
		fout << "DATASET UNSTRUCTURED_GRID\n";
		
		// Write control point coordinates
		fout << "POINTS " << cpts.size() << " float\n";
		for (const auto& point : cpts) {
			fout << point.coor[0] << " " << point.coor[1] << " " << point.coor[2] << "\n";
		}
		
		// Write scalar data associated with points
		fout << "POINT_DATA " << input.size() << "\n";
		fout << "SCALARS pact float 1\n";
		fout << "LOOKUP_TABLE default\n";
		for (const auto& value : input) {
			fout << value << "\n";
		}
		
		fout.close(); // Close the file
	} else {
		// Print an error message if the file cannot be opened
		cerr << "Error: Cannot open file " << fname << " for writing.\n";
	}
}

void NeuronGrowth::ToPETScVec(const vector<float>& input, Vec& petscVec)
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

void NeuronGrowth::ReadBezierElementProcess(const string& fn)
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

void NeuronGrowth::GaussInfo(int ng) {
    // Clear existing data
    Gpt.clear();
    wght.clear();

    // Resize vectors based on the number of Gauss points
    Gpt.resize(ng);
    wght.resize(ng);

    switch (ng) {
        case 2:
            Gpt = {0.2113248654051871, 0.7886751345948129};
            wght = {1.0, 1.0};
            break;

        case 3:
            Gpt = {0.1127016653792583, 0.5, 0.8872983346207417};
            wght = {0.5555555555555556, 0.8888888888888889, 0.5555555555555556};
            break;

        case 4:
            Gpt = {0.06943184420297371, 0.33000947820757187, 0.6699905217924281, 0.9305681557970262};
            wght = {0.3478548451374539, 0.6521451548625461, 0.6521451548625461, 0.3478548451374539};
            break;

        case 5:
            Gpt = {0.046910077030668, 0.2307653449471585, 0.5, 0.7692346550528415, 0.953089922969332};
            wght = {0.2369268850561891, 0.4786286704993665, 0.5688888888888889, 0.4786286704993665, 0.2369268850561891};
            break;

        default: // Fallback to 2 Gauss points
            Gpt = {0.2113248654051871, 0.7886751345948129};
            wght = {1.0, 1.0};
            break;
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

void NeuronGrowth::ApplyBoundaryCondition(
    const float bc_value, int pt_num, int variable_num,
    vector<vector<float>> &EMatrixSolve, vector<float> &EVectorSolve)
{
    int nen = EVectorSolve.size();
    int boundary_index = pt_num + variable_num * nen;

    // Update the RHS vector with the boundary condition value
    for (int j = 0; j < nen; ++j) {
        EVectorSolve[j] -= bc_value * EMatrixSolve[j][boundary_index];
        EMatrixSolve[j][boundary_index] = 0.0f; // Clear column
    }

    // Clear row and set diagonal
    fill(EMatrixSolve[boundary_index].begin(), EMatrixSolve[boundary_index].end(), 0.0f);
    EMatrixSolve[boundary_index][boundary_index] = 1.0f;

    // Set RHS boundary value
    EVectorSolve[boundary_index] = bc_value;
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

void NeuronGrowth::VisualizeVTK_ControlMesh(const vector<Vertex3D>& spt, int step, string fn)
{
	string fname;
	stringstream ss;
	ss << setw(6) << setfill('0') << step;
	ofstream fout;
	unsigned int i;
	fname = fn + "/controlmesh_" + ss.str() + ".vtk";
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << spt.size() << " float\n";
		for (i = 0; i < spt.size(); i++)
		{
			fout << spt[i].coor[0] << " " << spt[i].coor[1] << " " << spt[i].coor[2] << "\n";
		}
		// fout << "\nCELLS " << mesh.size() << " " << 9 * mesh.size() << '\n';
		// for (i = 0; i < mesh.size(); i++)
		// {
		// 	fout << "8 " << mesh[i].IEN[0] << " " << mesh[i].IEN[1] << " " << mesh[i].IEN[2] << " " << mesh[i].IEN[3]
		// 	     << " " << mesh[i].IEN[4] << " " << mesh[i].IEN[5] << " " << mesh[i].IEN[6] << " " << mesh[i].IEN[7] << '\n';
		// }
		// fout << "\nCELL_TYPES " << mesh.size() << '\n';
		// for (i = 0; i < mesh.size(); i++)
		// {
		// 	fout << "12\n";
		// }
		fout << "\nPOINT_DATA " << phi.size() << "\nSCALARS phi float 1\nLOOKUP_TABLE default\n";
		for (uint i = 0; i < phi.size(); i++)
		{
			fout << phi[i] /* +N_plus[i]+N_minus[i] */ << "\n";
		}
		fout << "\nSCALARS synaptogenesis float 1\nLOOKUP_TABLE default\n";
		for (uint i = 0; i < syn.size(); i++)
		{
			fout << syn[i] /* +N_plus[i]+N_minus[i] */ << "\n";
		}
		fout << "\nSCALARS tips float 1\nLOOKUP_TABLE default\n";
		for (uint i = 0; i < tips.size(); i++)
		{
			fout << tips[i] /* +N_plus[i]+N_minus[i] */ << "\n";
		}
		fout << "\nSCALARS tubulin float 1\nLOOKUP_TABLE default\n";
		for (uint i = 0; i < tub.size(); i++)
		{
			fout << tub[i] /* +N_plus[i]+N_minus[i] */ << "\n";
		}
		fout << "\nSCALARS theta float 1\nLOOKUP_TABLE default\n";
		for (uint i = 0; i < theta.size(); i++)
		{
			fout << theta[i] /* +N_plus[i]+N_minus[i] */ << "\n";
		}
		fout << "\nSCALARS tub_0 float 1\nLOOKUP_TABLE default\n";
		for (uint i = 0; i < tub_0.size(); i++)
		{
			fout << tub_0[i] /* +N_plus[i]+N_minus[i] */ << "\n";
		}
		fout << "\nSCALARS phi_0 float 1\nLOOKUP_TABLE default\n";
		for (uint i = 0; i < phi_0.size(); i++)
		{
			fout << phi_0[i] /* +N_plus[i]+N_minus[i] */ << "\n";
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

void NeuronGrowth::VisualizeVTK_PhysicalDomain_All(int step, string fn) {
	// Initialize containers for sample points, simulation results, and element connectivity for all variables
	vector<vector<array<float, 3>>> spt_all_5var(5); // Sample points for all 5 variables
	vector<vector<float>> sresult_all_5var(5); // Simulation results for all 5 variables
	vector<vector<array<int, 8>>> sele_all_5var(5); // Element connectivity for all 4 variables

	// Set current variable to phi and calculate its output variables
	N_0 = phi;
	CalculateVarsForOutput(spt_all_5var[0], sresult_all_5var[0], sele_all_5var[0]);

	// Set current variable to syn and calculate its output variables
	N_0 = syn;
	CalculateVarsForOutput(spt_all_5var[1], sresult_all_5var[1], sele_all_5var[1]);

	// Set current variable to tub and calculate its output variables
	N_0 = tub;
	CalculateVarsForOutput(spt_all_5var[2], sresult_all_5var[2], sele_all_5var[2]);

	// Set current variable to tips and calculate its output variables
	N_0 = tips;
	CalculateVarsForOutput(spt_all_5var[3], sresult_all_5var[3], sele_all_5var[3]);

	N_0 = theta;
	CalculateVarsForOutput(spt_all_5var[4], sresult_all_5var[4], sele_all_5var[4]);

	// If running on the master process (comRank == 0), write the VTK file for visualization
	if (comRank == 0) {
		WriteVTK_ALL(spt_all_5var[0], sresult_all_5var, sele_all_5var[0], step, fn);
	}
}

/**
 * @brief Writes an unstructured-grid VTK file with hexahedral cells and scalar fields.
 *
 * @param[in]  spt   A vector of 3D points (size = number of points).
 * @param[in]  sdisp A vector of scalar fields, each of size = number of points.
 *                   E.g., sdisp[0] is phi, sdisp[1] is synaptogenesis, etc.
 * @param[in]  sele  A vector of elements (each is an array<int,8> for a HEX),
 *                   specifying the 8 node indices of each cell.
 * @param[in]  step  Current timestep or iteration number, used in the filename.
 * @param[in]  fn    Base directory or filename prefix to which we append the VTK file name.
 *
 * This function produces an ASCII VTK file containing:
 *   - Unstructured grid with "POINTS", "CELLS", "CELL_TYPES"
 *   - Point-based scalar data ("POINT_DATA")
 *
 * The file is named:  "<fn>/physics_allparticle_<step>.vtk"
 *
 * @note If you need additional cell data or vector fields, you can extend this function accordingly.
 */
void NeuronGrowth::WriteVTK_ALL(const vector<array<float, 3>> &spt,
                                const vector<vector<float>> &sdisp,
                                const vector<array<int, 8>> &sele,
                                int step,
                                const string &fn)
{
    //--------------------------------------------------------------------------
    // 1) Construct the output file name
    //--------------------------------------------------------------------------
    stringstream ss;
    ss << step;  // Convert the step to string
    string fname = fn + "/physics_allparticle_" + ss.str() + ".vtk";

    //--------------------------------------------------------------------------
    // 2) Open the file for writing
    //--------------------------------------------------------------------------
    ofstream fout(fname);
    if (!fout.is_open()) {
        cerr << "Cannot open " << fname << " for writing!\n";
        return;
    }

    //--------------------------------------------------------------------------
    // 3) Preliminary checks to avoid mismatched data
    //--------------------------------------------------------------------------
    // Ensure we have at least one scalar field
    if (sdisp.empty()) {
        cerr << "Warning: sdisp is empty. No scalar fields will be written.\n";
    }

    // Check that all scalar fields have the same size as the number of points
    for (size_t varIndex = 0; varIndex < sdisp.size(); ++varIndex) {
        if (sdisp[varIndex].size() != spt.size()) {
            cerr << "Warning: sdisp[" << varIndex << "] has size "
                      << sdisp[varIndex].size() << " != " << spt.size()
                      << " (number of points). The output may be inconsistent.\n";
        }
    }

    //--------------------------------------------------------------------------
    // 4) Write the VTK file header
    //--------------------------------------------------------------------------
    fout << "# vtk DataFile Version 2.0\n";
    fout << "Hex test\n";
    fout << "ASCII\n";
    fout << "DATASET UNSTRUCTURED_GRID\n";

    //--------------------------------------------------------------------------
    // 5) Write the points
    //--------------------------------------------------------------------------
    fout << "POINTS " << spt.size() << " float\n";
    for (size_t i = 0; i < spt.size(); i++) {
        fout << spt[i][0] << " " << spt[i][1] << " " << spt[i][2] << "\n";
    }

    //--------------------------------------------------------------------------
    // 6) Write the cells (VTK requires "numCells" followed by "9*numCells" for HEX)
    //--------------------------------------------------------------------------
    fout << "\nCELLS " << sele.size() << " " << sele.size() * 9 << "\n";
    for (size_t i = 0; i < sele.size(); i++) {
        fout << "8 "     // each cell has 8 vertices
             << sele[i][0] << " " << sele[i][1] << " " << sele[i][2] << " " << sele[i][3] << " "
             << sele[i][4] << " " << sele[i][5] << " " << sele[i][6] << " " << sele[i][7] << "\n";
    }

    //--------------------------------------------------------------------------
    // 7) Write the cell types (12 corresponds to VTK_HEXAHEDRON)
    //--------------------------------------------------------------------------
    fout << "\nCELL_TYPES " << sele.size() << "\n";
    for (size_t i = 0; i < sele.size(); i++) {
        fout << "12\n";
    }

    //--------------------------------------------------------------------------
    // 8) Write scalar fields: "POINT_DATA" indicates data attached to each point
    //--------------------------------------------------------------------------
    // We'll map each sub-vector in sdisp to a named scalar field.
    // Adjust or extend if you have more or fewer than 5 scalars.
    static const char *scalarNames[] = {"phi", "synaptogenesis", "tubulin", "tips", "theta"};

    // Only write POINT_DATA if there's at least one scalar
    if (!sdisp.empty()) {
        fout << "\nPOINT_DATA " << spt.size() << "\n";

        // Write each field. If sdisp.size() > 5, you'll need more names or a dynamic approach.
        size_t numScalars = min(sdisp.size(), static_cast<size_t>(5));
        for (size_t varIndex = 0; varIndex < numScalars; ++varIndex) {
            fout << "\nSCALARS " << scalarNames[varIndex] << " float 1\n";
            fout << "LOOKUP_TABLE default\n";
            for (size_t i = 0; i < sdisp[varIndex].size(); i++) {
                fout << sdisp[varIndex][i] << "\n";
            }
        }
    }

    //--------------------------------------------------------------------------
    // 9) Clean up and close
    //--------------------------------------------------------------------------
    fout.close();
}

/**
 * @brief Reads a VTK file in ASCII format, extracting 3D points (cpts)
 *        and up to five scalar fields (phi, syn, tub, tips, theta).
 */
bool NeuronGrowth::ReadVTK(const string &filename)
{
    ifstream fin(filename);
    if (!fin.is_open()) {
        cerr << "Cannot open VTK file: " << filename << endl;
        return false;
    }

    cpts.clear();
    phi.clear();
    syn.clear();
    tub.clear();
    tips.clear();
    theta.clear();
	phi_0.clear();
	tub_0.clear();

    string line;
    bool foundPoints = false;
    size_t numPoints = 0;

    // We'll store the names in the order we expect them
    const vector<string> scalarNames = {
        "phi", "synaptogenesis", "tubulin", "tips", "theta", "phi_0", "tub_0"
    };
    // We'll keep references in an array so we can fill them in by index
    vector<vector<float>*> scalarArrays = {
        &phi, &syn, &tub, &tips, &theta, &phi_0, &tub_0
    };

    //--------------------------------------------------------------------------
    // 1) Parse lines until we find "POINTS <N> float"
    //--------------------------------------------------------------------------
    while (getline(fin, line)) {
        // e.g. "POINTS 123 float"
        if (line.rfind("POINTS ", 0) == 0) {
            // Parse the number of points
            stringstream ss(line);
            string dummy;
            ss >> dummy;            // "POINTS"
            ss >> numPoints;        // e.g. 123
            // skip "float"
            foundPoints = true;

            // read the points
            cpts.resize(numPoints);
            for (size_t i = 0; i < numPoints; i++) {
                if (!getline(fin, line)) {
                    cerr << "Error reading points (EOF encountered)\n";
                    return false;
                }
                stringstream ssp(line);
                float x, y, z;
                ssp >> x >> y >> z;
				cpts[i].coor[0] = x;
				cpts[i].coor[1] = y;
				cpts[i].coor[2] = z;
            }
            break; // done reading points
        }
    }

    if (!foundPoints || numPoints == 0) {
        cerr << "No POINTS section found or zero points.\n";
        return false;
    }

    //--------------------------------------------------------------------------
    // 2) Now look for "POINT_DATA <numPoints>"
    //--------------------------------------------------------------------------
    bool foundPointData = false;
    while (getline(fin, line)) {
        if (line.rfind("POINT_DATA ", 0) == 0) {
            // parse how many
            stringstream ss(line);
            string dummy;
            size_t checkN;
            ss >> dummy;      // "POINT_DATA"
            ss >> checkN;     // e.g. 123

            if (checkN != numPoints) {
                cerr << "Warning: POINT_DATA " << checkN
                          << " != #points (" << numPoints << ")\n";
            }
            foundPointData = true;
            break;
        }
    }
    if (!foundPointData) {
        cerr << "No POINT_DATA section found.\n";
        return false;
    }

    //--------------------------------------------------------------------------
    // 3) Read the scalars in the known order (phi, syn, tub, tips, theta, phi_0, tub_0)
    //    For each scalar: 
    //         SCALARS <name> float 1
    //         LOOKUP_TABLE default
    //         <numPoints> lines
    //--------------------------------------------------------------------------
    size_t foundScalarsCount = 0;
    while (foundScalarsCount < scalarNames.size()) {
        if (!getline(fin, line)) break; // no more lines

        // We expect: "SCALARS phi float 1"
        if (line.rfind("SCALARS ", 0) == 0) {
            // parse the name
            stringstream ss(line);
            string dummy, scalarName, type;
            int components = 0;
            ss >> dummy;       // "SCALARS"
            ss >> scalarName;  // e.g. "phi"
            ss >> type;        // e.g. "float"
            ss >> components;  // e.g. 1

            // next line should be "LOOKUP_TABLE default"
            if (!getline(fin, line)) {
                cerr << "EOF reading SCALARS " << scalarName << endl;
                return false;
            }
            // ignore the "LOOKUP_TABLE" line

            // Find which scalar index we have
            auto it = find(scalarNames.begin(), scalarNames.end(), scalarName);
            if (it == scalarNames.end()) {
                // Unrecognized scalar; skip lines anyway
                for (size_t i = 0; i < numPoints; i++) {
                    if (!getline(fin, line)) {
                        cerr << "EOF skipping unknown scalar data.\n";
                        return false;
                    }
                }
            } else {
                size_t index = distance(scalarNames.begin(), it);
                scalarArrays[index]->resize(numPoints);

                // read numPoints lines
                for (size_t i = 0; i < numPoints; i++) {
                    if (!getline(fin, line)) {
                        cerr << "EOF reading scalar " << scalarName << "\n";
                        return false;
                    }
                    float val = stof(line);
                    (*scalarArrays[index])[i] = val;
                }
                foundScalarsCount++;
            }
        }
        // else skip lines that are not "SCALARS ..."
    }

    fin.close();
    return true;
}

void NeuronGrowth::PointFormValue(vector<float>& Nx,
								const vector<float>& U,
								float& Value)
{
	Value = 0;

	for (int j = 0; j < Nx.size(); j++)
		Value += U[j] * Nx[j];
}

void NeuronGrowth::PointFormGrad(vector<array<float, 3>>& dNdx,
								const vector<float>& U,
								float Value[2])
{
	for (int j = 0; j < 2; j++)
		Value[j] = 0.;

	for (int j = 0; j < 2; j++)
		for (int k = 0; k < dNdx.size(); k++)
			Value[j] += U[k] * dNdx[k][j];
}

void NeuronGrowth::PointFormHess(vector<array<array<float, 2>, 2>>& d2Ndx2,
								const vector<float>& U,
								float Value[2][2])
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

void NeuronGrowth::ElementValue(const vector<float>& Nx,
							const vector<float>& value_node,
							float& value)
{
	value = 0.;
	for (int i = 0; i < Nx.size(); i++)
		value += value_node[i] * Nx[i];
}

void NeuronGrowth::ElementValueAll(const vector<float>& Nx,
								const vector<float>& elePhiGuess, float& elePG,
								const vector<float>& elePhi, float& eleP,
								const vector<float>& eleSyn, float& eleS,
								const vector<float>& eleTips, float& eleTp,
								const vector<float>& eleTubulin, float& eleTb,
								const vector<float>& eleEpsilon, float& eleEP,
								const vector<float>& eleEpsilonP, float& eleEEP)
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

void NeuronGrowth::ElementDeriv(const int nen,
							vector<array<float, 3>>& dNdx,
							const vector<float>& value_node,
							float& dVdx, float& dVdy, float& dVdz)
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

void NeuronGrowth::ElementDerivAll(const int nen,
								vector<array<float, 3>>& dNdx,
								const vector<float>& elePhiGuess, float& dPGdx, float& dPGdy,
								const vector<float>& eleTheta, float& dThedx, float& dThedy,
								const vector<float>& eleEpsilon, float& dAdx, float& dAdy,
								const vector<float>& eleEpsilonP, float& dAPdx, float& dAPdy)
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

inline void NeuronGrowth::ElementEvaluationAll_phi(
    int nen,
    const vector<float> &Nx,
    const vector<array<float, 3>> &dNdx,
    const vector<vector<float>> &eleVal,
    vector<float> &vars)
{
    // Reset vars
    vars[2] = vars[5] = vars[6] = vars[9] = vars[7] = vars[12] = vars[15] = 0.0f;
    vars[3] = vars[4] = vars[18] = vars[13] = vars[14] = vars[19] = vars[16] = vars[17] = vars[20] = 0.0f;

    // Local references to reduce indexing overhead
    float &phi_g     = vars[2];  // elePG
    float &phi       = vars[5];  // eleP
    float &syn       = vars[6];  // eleS
    float &tips      = vars[9];  // eleTp
    float &tubulin   = vars[7];  // eleTb
    // float &epsilon  = vars[12]; // If needed
    // float &epsilonP = vars[15]; // If needed

    float &dPGdx = vars[3];
    float &dPGdy = vars[4];
    float &dPGdz = vars[18];
    // Similarly for dAdx, dAdy, dAdz if you need them

    // Unroll loops if nen is fixed and known (optional)
    for (int i = 0; i < nen; i++) {
        float Nx_i = Nx[i];
        phi_g   += eleVal[0][i] * Nx_i;
        phi     += eleVal[1][i] * Nx_i;
        syn     += eleVal[2][i] * Nx_i;
        tips    += eleVal[5][i] * Nx_i;
        tubulin += eleVal[3][i] * Nx_i;
        // epsilon += eleVal[6][i] * Nx_i; 
        // epsilonP+= eleVal[7][i] * Nx_i;

        const float dNdx_i0 = dNdx[i][0];
        const float dNdx_i1 = dNdx[i][1];
        const float dNdx_i2 = dNdx[i][2];
        float val0 = eleVal[0][i];
        dPGdx += val0 * dNdx_i0;
        dPGdy += val0 * dNdx_i1;
        dPGdz += val0 * dNdx_i2;
        // If needed for epsilon and epsilonP derivatives, do similarly here.
    }
}

void NeuronGrowth::ElementEvaluationAll_phi_test(const uint &nen, 
	const vector<float> &Nx, 
	const vector<array<float, 3>> &dNdx, 
	const vector<float> &elePhiGuess, 
	vector<float> &vars) 
{
    // Reset variables
    vars[0] = 0.;  // elePG (element value for PhiGuess)
    vars[1] = 0.;  // dPGdx (derivative of PhiGuess w.r.t. x)
    vars[2] = 0.;  // dPGdy (derivative of PhiGuess w.r.t. y)
    vars[3] = 0.;  // dPGdz (derivative of PhiGuess w.r.t. z)

    for (size_t i = 0; i < nen; i++) {
        // Accumulate element values and derivatives
        vars[0] += elePhiGuess[i] * Nx[i];           // Element value for PhiGuess
        vars[1] += elePhiGuess[i] * dNdx[i][0];     // Derivative w.r.t. x
        vars[2] += elePhiGuess[i] * dNdx[i][1];     // Derivative w.r.t. y
        vars[3] += elePhiGuess[i] * dNdx[i][2];     // Derivative w.r.t. z
    }
}

void NeuronGrowth::ElementEvaluationAll_phi_opt(
    const uint &nen,
    const vector<float> &Nx,
    const vector<array<float, 3>> &dNdx,
    const vector<float> &elePhiGuess,
    vector<float> &vars
) {
    // Initialize variables in a single step
    vars[0] = vars[1] = vars[2] = vars[3] = 0.0f;

    // Use loop unrolling for improved performance
    size_t i = 0;
    for (; i + 3 < nen; i += 4) {
        vars[0] += elePhiGuess[i] * Nx[i] +
                   elePhiGuess[i + 1] * Nx[i + 1] +
                   elePhiGuess[i + 2] * Nx[i + 2] +
                   elePhiGuess[i + 3] * Nx[i + 3];

        vars[1] += elePhiGuess[i] * dNdx[i][0] +
                   elePhiGuess[i + 1] * dNdx[i + 1][0] +
                   elePhiGuess[i + 2] * dNdx[i + 2][0] +
                   elePhiGuess[i + 3] * dNdx[i + 3][0];

        vars[2] += elePhiGuess[i] * dNdx[i][1] +
                   elePhiGuess[i + 1] * dNdx[i + 1][1] +
                   elePhiGuess[i + 2] * dNdx[i + 2][1] +
                   elePhiGuess[i + 3] * dNdx[i + 3][1];

        vars[3] += elePhiGuess[i] * dNdx[i][2] +
                   elePhiGuess[i + 1] * dNdx[i + 1][2] +
                   elePhiGuess[i + 2] * dNdx[i + 2][2] +
                   elePhiGuess[i + 3] * dNdx[i + 3][2];
    }

    // Handle remaining iterations
    for (; i < nen; ++i) {
        vars[0] += elePhiGuess[i] * Nx[i];
        vars[1] += elePhiGuess[i] * dNdx[i][0];
        vars[2] += elePhiGuess[i] * dNdx[i][1];
        vars[3] += elePhiGuess[i] * dNdx[i][2];
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

void NeuronGrowth::PrepareBasis() {

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

void NeuronGrowth::PreparePhaseField() 
{
    // Clear previous data
    pre_eleEP.clear();
    pre_eleEEP.clear();
    pre_dAdx.clear();
    pre_dAdy.clear();
    pre_dAdz.clear();  // Added for 3D
    pre_eleP.clear();
    pre_eleTh.clear();
    pre_eleMp.clear();
    pre_C1.clear();

    // Initialize variables
    float eleEP(0), eleEEP(0), dAdx(0), dAdy(0), dAdz(0), // Added dAdz
          eleP(0), eleTh(0), eleS(0), eleTb(0), eleTp(0), eleE(0);
    uint ind(0);

    // Loop over elements
    for (size_t e = 0; e < bzmesh_process.size(); e++) {
        size_t nen = bzmesh_process[e].IEN.size();

        // Initialize element-wise variables
        vector<float> elePhi(nen, 0), eleTheta(nen, 0), eleEpsilon(nen, 0), eleEpsilonP(nen, 0);
        vector<float> eleSyn(nen, 0), eleTubulin(nen, 0), eleTips(nen, 0);

        // Extract nodal values for the current element
        for (size_t i = 0; i < nen; i++) {
            size_t nodeIndex = bzmesh_process[e].IEN[i];
            elePhi[i]    = phi[nodeIndex];       // elePhi
            eleTheta[i]  = theta[nodeIndex];     // eleTheta
            eleSyn[i]    = syn[nodeIndex];       // eleSyn
            eleTubulin[i]= tub[nodeIndex];       // eleTubulin
            eleTips[i]   = tips[nodeIndex];      // eleTips
        }

        // Loop over Gaussian quadrature points in 3D
        for (size_t i = 0; i < Gpt.size(); i++) {
            for (size_t j = 0; j < Gpt.size(); j++) {
                for (size_t k = 0; k < Gpt.size(); k++) { // Added k-loop for 3D

                    // Evaluate orientation and compute element variables
                    // EvaluateOrientation(nen, pre_Nx[ind], pre_dNdx[ind], elePhi, eleTheta, eleEpsilon, eleEpsilonP);
					// EvaluateOrientation(nen, pre_Nx[ind], pre_dNdx[ind], elePhi, eleTheta, eleAniso, dA_dPdx, dA_dPdy, dA_dPdz)
					EvaluateOrientation_prev(nen, Nx, dNdx, elePhi, eleTheta, eleEpsilon, eleEpsilonP);
                    // Compute element values and derivatives
                    ElementValue(pre_Nx[ind], eleEpsilon, eleEP);
                    pre_eleEP.push_back(eleEP);

                    ElementValue(pre_Nx[ind], eleEpsilonP, eleEEP);
                    pre_eleEEP.push_back(eleEEP);

                    ElementDeriv(nen, pre_dNdx[ind], eleEpsilonP, dAdx, dAdy, dAdz); // Modified to include dAdz
                    pre_dAdx.push_back(dAdx);
                    pre_dAdy.push_back(dAdy);
                    pre_dAdz.push_back(dAdz); // Added for 3D

                    ElementValue(pre_Nx[ind], elePhi, eleP);
                    pre_eleP.push_back(eleP);

                    ElementValue(pre_Nx[ind], eleTheta, eleTh);
                    pre_eleTh.push_back(eleTh);

                    ElementValue(pre_Nx[ind], eleSyn, eleS);
                    ElementValue(pre_Nx[ind], eleTubulin, eleTb);
                    ElementValue(pre_Nx[ind], eleTips, eleTp);

                    // Adjust assembly and disassembly rates based on detected tips
                    if (n < 50) { // Assuming 'n' is a time step or iteration count; ensure it's properly defined
                        eleE = alphaOverPi * atan(gamma * (c_opt - eleS));
                        pre_eleMp.push_back(M_neurite);
                    } else {
                        if (eleTp != 0) {
                            eleE = alphaOverPi * atan(gamma * Regular_Heiviside_fun(50 * eleTb) * (c_opt - eleS));
                            if (eleTp < 0) {
                                pre_eleMp.push_back(M_axon);
                            } else {
                                pre_eleMp.push_back(M_neurite);
                            }
                        } else {
                            eleE = alphaOverPi * atan(gamma * Regular_Heiviside_fun(r * eleTb - g) * (c_opt - eleS));
                            pre_eleMp.push_back(5);
                        }
                    }

                    // Calculate C1 variable for the phase-field energy term
                    float C1 = eleE - pre_C0[ind]; // Update pre_C0[ind] if necessary to include 3D effects
                    pre_C1.push_back(C1);

                    // Increment the index for precomputed variables
                    ind += 1;
                }
            }
        }
    }
}

void NeuronGrowth::PreparePhaseField_SNES() 
{
    uint ind(0);

    // Resize `pre_eleVal` to match the structure: [elements][fields][nodes]
    size_t numElements = bzmesh_process.size();
    size_t numFields = 8; // Number of fields as indicated by your code

    // Determine maximum `nen` (number of nodes per element) dynamically
    size_t maxNen = 0;
    for (const auto& elem : bzmesh_process) {
        maxNen = max(maxNen, elem.IEN.size());
    }

    pre_eleVal.resize(numElements, vector<vector<float>>(numFields, vector<float>(maxNen, 0.0f)));

    // Loop over elements
    for (size_t e = 0; e < numElements; e++) {
        const auto& IEN = bzmesh_process[e].IEN;
        size_t nen = IEN.size(); // Get number of nodes for the current element

        for (size_t ii = 0; ii < nen; ii++) {
            int A = IEN[ii];

            pre_eleVal[e][1][ii] = phi[A];    // phi
            pre_eleVal[e][2][ii] = syn[A];    // syn
            pre_eleVal[e][3][ii] = tub[A];    // tub
            pre_eleVal[e][4][ii] = theta[A];  // theta
            pre_eleVal[e][5][ii] = tips[A];   // tips
            pre_eleVal[e][6][ii] = 0.0f;      // epsilon
            pre_eleVal[e][7][ii] = 0.0f;      // epsilonP
        }
    }
}

/**
 * @brief Precomputes element-level data that remain constant
 *        during each solve iteration, except for the final assembled phi guess.
 *
 * This includes:
 *  - Element shape function data for old fields (phi, syn, tub, theta, tips).
 *  - Possible anisotropy/orientation factors (if n > 0).
 *  - Computed partial terms like pre_eleP (old-phase accumulation).
 *  - Precomputed assembly/disassembly rates (pre_eleMp) and phase coefficients (pre_C1).
 *
 * The only thing that changes during the solve (in FormFunction_phi) is (float)Parray[A],
 * so we cache everything else here to reduce computation cost inside the solver loop.
 */
void NeuronGrowth::PreparePhaseField_SNES_preComputed() 
{
    //--------------------------------------------------------------------------
    // 1. Clear old data in precomputed arrays to avoid stale or leftover values
    //--------------------------------------------------------------------------
    pre_eleP.clear();
    pre_eleAniso.clear();
    pre_dA_dPdx.clear();
    pre_dA_dPdy.clear();
    pre_dA_dPdz.clear();
    pre_eleMp.clear();
    pre_C1.clear();

    //--------------------------------------------------------------------------
    // 2. Prepare the size for pre_eleVal and other containers
    //--------------------------------------------------------------------------
    size_t numElements = bzmesh_process.size();

    // Determine maximum number of nodes per element (nen)
    size_t maxNen = 0;
    for (const auto &elem : bzmesh_process) {
        maxNen = max(maxNen, elem.IEN.size());
    }

    // pre_eleVal structure: [element_index][field_index][node_index]
    // 6 fields: 0->(phi guess), 1->phi, 2->syn, 3->tub, 4->theta, 5->tips
    pre_eleVal.resize(numElements,
                      vector<vector<float>>(8, vector<float>(maxNen, 0.0f)));

    // Resize precomputed arrays to match the total number of Gauss points
    // (which is pre_Nx.size())
    pre_eleAniso.resize(pre_Nx.size(), 0.0f);
    pre_dA_dPdx.resize(pre_Nx.size(), 0.0f);
    pre_dA_dPdy.resize(pre_Nx.size(), 0.0f);
    pre_dA_dPdz.resize(pre_Nx.size(), 0.0f);
    pre_eleP.resize(pre_Nx.size(), 0.0f);
    pre_eleMp.resize(pre_Nx.size(), 0.0f);
    pre_C1.resize(pre_Nx.size(), 0.0f);

    const size_t gptSize = Gpt.size();

    //--------------------------------------------------------------------------
    // 3. Fill element fields with OLD data for phi, syn, tub, tips, etc.
    //    Then, for each Gauss point, compute orientation & partial PDE terms
    //--------------------------------------------------------------------------
    size_t ind = 0;  // This will index into pre_Nx, pre_dNdx, etc. for each Gauss point

    // Loop over elements
    for (size_t e = 0; e < numElements; e++) {
        // (a) Retrieve number of nodes for this element
        auto &IEN = bzmesh_process[e].IEN;
        const size_t nen = IEN.size();

        // (b) For each node, fill in old fields in pre_eleVal (phi, syn, tub, etc.)
        for (size_t ii = 0; ii < nen; ii++) {
            const size_t A = IEN[ii];
            pre_eleVal[e][1][ii] = phi[A];		// Field 1 -> phi
            pre_eleVal[e][2][ii] = syn[A];		// Field 2 -> syn
            pre_eleVal[e][3][ii] = tub[A];		// Field 3 -> tub
            pre_eleVal[e][4][ii] = theta[A];	// Field 4 -> theta
            pre_eleVal[e][5][ii] = tips[A];		// Field 5 -> tips
        }

        // (c) For each Gauss point (i,j,k), compute orientation (if needed) and
        //     evaluate old-phase data in ElementEvaluationAll_phi.
        for (size_t i = 0; i < gptSize; i++) {
            for (size_t j = 0; j < gptSize; j++) {
                for (size_t k = 0; k < gptSize; k++) {

                    // i. Evaluate anisotropy if n > 0
                    if (n > 0) {
                        EvaluateOrientation(static_cast<int>(nen),
                                            pre_Nx[ind],
                                            pre_dNdx[ind],
                                            pre_eleVal[e][1],    // old phi
                                            pre_eleVal[e][4],    // old theta
                                            pre_eleAniso[ind],
                                            pre_dA_dPdx[ind],
                                            pre_dA_dPdy[ind],
                                            pre_dA_dPdz[ind]);
                    }

                    // ii. Evaluate old-phase fields in vars[5], etc.
                    //     This also accumulates phi_guess if needed (in vars[2]),
                    //     but here we typically only use the old-phase portion.
                    ElementEvaluationAll_phi(static_cast<int>(nen),
                                             pre_Nx[ind],
                                             pre_dNdx[ind],
                                             pre_eleVal[e],
                                             vars);

                    // store old-phase value in pre_eleP
                    pre_eleP[ind] = vars[5]; // vars[5] -> old phi

                    // iii. Decide assembly rate based on tips (vars[9]) & syn (vars[6])
                    float eleE = 0.0f;
                    if (n < 0) {
                        // negative n logic
                        eleE = alphaOverPi * atan(gamma * (1 - vars[6]));
						pre_eleMp[ind] = M_neurite;
                    } else {
                        // positive n logic: check tips
                        if (vars[9] > tip_threshold) {
                            eleE = alphaOverPi * atan(gamma * 1.0f * (1 - vars[6]));
                            // pre_eleMp[ind] = M_axon;
                            pre_eleMp[ind] = M_neurite;
                        } else {
                            eleE = alphaOverPi * atan(gamma * 0.05f * (1 - vars[6]));
                            pre_eleMp[ind] = M_phi;
                        }
                    }

                    // iv. Compute pre_C1 as (eleE - pre_C0[ind]) if pre_C0 is an array of size == pre_Nx.size().
                    pre_C1[ind] = eleE - pre_C0[ind];

                    // (d) Increment index to move to the next Gauss point
                    ind++;
                }
            }
        }
    }
}

void NeuronGrowth::PrepareTermSource() {
    // Precompute values outside of loops for efficiency
    const size_t pointsPerElement = Gpt.size() * Gpt.size() * Gpt.size();
    const size_t totalPoints = bzmesh_process.size() * pointsPerElement;

    // Reserve space for efficiency
    pre_term_source.clear();
    pre_term_source.reserve(totalPoints);

    // Avoid redundant computation of `source_coeff / sum_grad_phi0_global`
    const double normalizedSourceCoeff = source_coeff / sum_grad_phi0_global;

    // Single loop to flatten index computations
    for (size_t ind = 0; ind < totalPoints; ++ind) {
        // Compute and append term source
        pre_term_source.emplace_back(pre_mag_grad_phi0[ind] * normalizedSourceCoeff);
    }
}

// Evaluate energy based on synaptogenesis
void NeuronGrowth::EvaluateEnergy(const int nen, const vector<float>& Nx, const vector<float>& eleSyn, vector<float>& E) {
    // Precompute constants for efficiency
    const float scale = alphaOverPi * gamma;

    for (int i = 0; i < nen; i++) {
        // Compute energy contribution for each element
        E[i] = alphaOverPi * atan(scale * (1 - eleSyn[i]));
    }
}

// Smooth approximation of the Heaviside function
float NeuronGrowth::Regular_Heiviside_fun(float x) {
    const float epsilon = 1e-5; // Epsilon to ensure stability and smoothness
    return 0.5 * (1 + (2 / PI) * atan(x / epsilon));
}

void NeuronGrowth::EvaluateOrientation_prev(const uint& nen, const vector<float>& Nx, const vector<array<float, 3>>& dNdx,
	const vector<float>& elePhi, const vector<float>& eleTheta, vector<float>& eleEpsilon, vector<float>& eleEpsilonP)
{
	for (size_t i = 0; i < nen; i++) {
		// Compute the gradient magnitude and the orientation (phi gradient direction)
		float gradX = elePhi[i] * dNdx[i][0];
		float gradY = elePhi[i] * dNdx[i][1];
		float gradZ = elePhi[i] * dNdx[i][2];
		float magnitude = sqrt(gradX * gradX + gradY * gradY + gradZ * gradZ + 1e-8f); // Adding a small value to avoid division by zero
		
		// Compute theta using spherical coordinates
		float theta = atan2(sqrt(gradX * gradX + gradY * gradY), gradZ);  // Polar angle
		float phi = atan2(gradY, gradX);                                  // Azimuthal angle

		// Orientation terms
		float angleDifference = phi - eleTheta[i] * Nx[i];

		// Epsilon and its derivative
		eleEpsilon[i] += epsilonb * (1.0f + delta * cos(aniso * angleDifference));
		eleEpsilonP[i] += -epsilonb * (aniso * delta * sin(aniso * angleDifference));
	}
}

void NeuronGrowth::EvaluateOrientation(
    const int nen, 
    const vector<float>& Nx, 
    const vector<array<float, 3>>& dNdx, 
    const vector<float>& elePhi, 
    const vector<float>& eleTheta, 
    float& eleAniso, 
    float& dA_dPdx, 
    float& dA_dPdy, 
    float& dA_dPdz
) {
    // Initialization
    dA_dPdx = dA_dPdy = dA_dPdz = eleAniso = 0.0f;
    const float stable = 1e-2;

    auto sanitize_value = [](float value, const string& name) -> float {
        if (isnan(value) || abs(value) > 10.0f) {
            cout << "Invalid " << name << ": " << value << endl;
            return 0.0f;
        }
        return value;
    };

    for (int i = 0; i < nen; ++i) {
        // Compute first-order derivatives
        float dPdx1 = elePhi[i] * dNdx[i][0];
        float dPdy1 = elePhi[i] * dNdx[i][1];
        float dPdz1 = elePhi[i] * dNdx[i][2];

        // Compute higher-order derivatives
        float dPdx2 = dPdx1 * dPdx1, dPdy2 = dPdy1 * dPdy1, dPdz2 = dPdz1 * dPdz1;
        float dPdx3 = dPdx2 * dPdx1, dPdy3 = dPdy2 * dPdy1, dPdz3 = dPdz2 * dPdz1;
        float dPdx4 = dPdx3 * dPdx1, dPdy4 = dPdy3 * dPdy1, dPdz4 = dPdz3 * dPdz1;

        // Aggregate contributions for fourth-order and magnitude terms
        float dP4 = dPdx4 + dPdy4 + dPdz4;
        float md4 = (dPdx2 + dPdy2 + dPdz2) * (dPdx2 + dPdy2 + dPdz2);

        // Compute eleAniso contribution
        float tmp_aniso = epsilonb * (1 - 3 * delta) + epsilonb * 4 * delta * dP4 / (md4 + stable);
        eleAniso += sanitize_value(tmp_aniso, "eleAniso");

        // Compute contributions to dA_dPdx
        float C1x = dPdy4 + dPdz4, C2x = dPdy2 + dPdz2;
        float tmp_dA_dPdx = epsilonb * 4 * delta * (
            (4 * dPdx3) / ((C2x + dPdx2) * (C2x + dPdx2) + stable) -
            (4 * dPdx1 * (C1x + dPdx4)) / ((C2x + dPdx2) * (C2x + dPdx2) * (C2x + dPdx2) + stable)
        );
        dA_dPdx += sanitize_value(tmp_dA_dPdx, "dA_dPdx");

        // Compute contributions to dA_dPdy
        float C1y = dPdx4 + dPdz4, C2y = dPdx2 + dPdz2;
        float tmp_dA_dPdy = epsilonb * 4 * delta * (
            (4 * dPdy3) / ((C2y + dPdy2) * (C2y + dPdy2) + stable) -
            (4 * dPdy1 * (C1y + dPdy4)) / ((C2y + dPdy2) * (C2y + dPdy2) * (C2y + dPdy2) + stable)
        );
        dA_dPdy += sanitize_value(tmp_dA_dPdy, "dA_dPdy");

        // Compute contributions to dA_dPdz
        float C1z = dPdx4 + dPdy4, C2z = dPdx2 + dPdy2;
        float tmp_dA_dPdz = epsilonb * 4 * delta * (
            (4 * dPdz3) / ((C2z + dPdz2) * (C2z + dPdz2) + stable) -
            (4 * dPdz1 * (C1z + dPdz4)) / ((C2z + dPdz2) * (C2z + dPdz2) * (C2z + dPdz2) + stable)
        );
        dA_dPdz += sanitize_value(tmp_dA_dPdz, "dA_dPdz");
    }
}

void NeuronGrowth::EvaluateOrientation_old(const int nen, const vector<float> &Nx, const vector<array<float, 3>> &dNdx, const vector<float> elePhi, const vector<float> eleTheta,  float& eleAniso, float& dA_dPdx, float& dA_dPdy, float& dA_dPdz)
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
				cout << "nan 1: " << tmp << endl;
				tmp = 0;
			}
			eleAniso += tmp;

			float C1x = dPdy4 + dPdz4;
			float C2x = dPdy2 + dPdz2;
			tmp = epsilonb * 4 * delta * ( (4 * dPdx3) / (pow(C2x + dPdx2, 2) + stable)
				- (4 * dPdx1 * (C1x + dPdx4)) / (pow((C2x + dPdx2), 3) + stable) );
			if ((isnan(tmp) == 1) || (abs(tmp) > 10)) {
				cout << "nan 2: " << tmp << endl;
				tmp = 0;
			}
			dA_dPdx += tmp;

			float C1y = dPdx4 + dPdz4;
			float C2y = dPdx2 + dPdz2;
			tmp =  epsilonb * 4 * delta * ( (4 * dPdy3) / (pow(C2y + dPdy2, 2) + stable)
				- (4 * dPdy1 * (C1y + dPdy4)) / (pow((C2y + dPdy2), 3) + stable) );
			if ((isnan(tmp) == 1) || (abs(tmp) > 10)) {
				cout << "nan 3: " << tmp << endl;
				tmp = 0;
			}
			dA_dPdy += tmp;
		
			float C1z = dPdx4 + dPdy4;
			float C2z = dPdx2 + dPdy2;
			tmp =  epsilonb * 4 * delta * ( (4 * dPdz3) / (pow(C2z + dPdz2, 2) + stable)
				- (4 * dPdz1 * (C1z + dPdz4)) / (pow((C2z + dPdz2), 3) + stable) );
			if ((isnan(tmp) == 1) || (abs(tmp) > 10)) {
				cout << "nan 4: " << tmp << endl;
				tmp = 0;
			}
			dA_dPdz += tmp;
		// }
	}
}

void NeuronGrowth::EvaluateOrientationSpherical(
		const int nen,                    // Number of nodes in the element
		const vector<float>& Nx,          // Basis function values
		const vector<array<float, 3>>& dNdx, // Basis function derivatives
		const vector<float>& elePhi,      // Phase field values
		const vector<float>& elePolar,    // Polar angle values
		const vector<float>& eleAzimuth,  // Azimuthal angle values
		float& eleEpsilon,                // Output epsilon value
		float dEdp,                       // Derivative w.r.t polar angle
		float dEda                        // Derivative w.r.t azimuthal angle
)
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

void NeuronGrowth::BuildLinearSystemProcessNG_phi(const vector<Vertex3D> &cpts) {
    /* Build linear system in each process */
    int ind = 0; // Pre-calculated variable index

    vector<vector<float>> EMatrixSolve;
    vector<float> EVectorSolve;
    vector<float> elePhiGuess;

    eleVal.resize(4); // Assuming 4 fields are used in eleVal

    for (auto &val : eleVal) {
        val.resize(64, 0.0f); // Pre-allocate for 64 nodes (adjust as needed)
    }

    for (size_t e = 0; e < bzmesh_process.size(); ++e) {
        const int nen = bzmesh_process[e].IEN.size(); // Number of nodes in the element

        // Resize vectors for current element
        // eleMphi.assign(nen, 0.0f);
        // elePhiGuess.assign(nen, 0.0f);
        EVectorSolve.assign(nen, 0.0f);
        EMatrixSolve.assign(nen, vector<float>(nen, 0.0f));

        // Resize and initialize vectors
        eleVal[0].resize(nen, 0.0f);
        eleVal[1].assign(nen, 0.0f); // Initialize once if always 0
        eleVal[2].assign(nen, 0.0f);
        eleVal[3].assign(nen, 0.0f);

        // Initialize element values for PhiGuess
        for (int i = 0; i < nen; ++i) {
            // elePhiGuess[i] = phi[bzmesh_process[e].IEN[i]];
			eleVal[0][i] = phi[bzmesh_process[e].IEN[i]]; // Set PhiGuess
        }

		// loop through gaussian quadrature points
		const size_t GptSize = Gpt.size();
		for (int i = 0; i < GptSize; i++) {
			for (int j = 0; j < GptSize; j++) {
				for (int k = 0; k < GptSize; k++) {
					// ElementEvaluationAll_phi_test(nen, pre_Nx[ind], pre_dNdx[ind], eleVal[0], vars);
					ElementEvaluationAll_phi_opt(nen, pre_Nx[ind], pre_dNdx[ind], eleVal[0], vars);
					// vars[0] = vars[8] - pre_C0_sp[ind];
					// float tau = sqrt(pow(vars[3],2) + pow(vars[4],2));
					// float mode_grad_p_2 = pow(sqrt(pow(vars[3],2) + pow(vars[4],2) + pow(vars[18],2)), 2);
					// loop through control points
					for (int m = 0; m < nen; m++) {
						// float Nx_m = pre_Nx[ind][m];
						EVectorSolve[m] += (vars[0] * pre_Nx[ind][m] - dt * pre_eleMp[ind] *\
							(
							// (- pre_eleEP[ind] * pre_eleEP[ind] * (vars[1] * pre_dNdx[ind][m][0] + vars[2] * pre_dNdx[ind][m][1])) -\
							// (- pre_dAdx[ind] * pre_eleEEP[ind] * vars[2] * pre_dNdx[ind][m][0]) +
							// (- pre_dAdy[ind] * pre_eleEEP[ind] * vars[1] * pre_dNdx[ind][m][1]) +
							(- vars[0] * vars[0] * vars[0] + (1 - pre_C1[ind]) * vars[0] * vars[0] + pre_C1[ind] * vars[0]) * pre_Nx[ind][m])
							- pre_eleP[ind] * pre_Nx[ind][m]) * pre_detJ[ind];

						// loop through 16 control points
						for (int n = 0; n < nen; n++) {
								// float Nx_n = pre_Nx[ind][n];
							EMatrixSolve[m][n] += (pre_Nx[ind][m] * pre_Nx[ind][n] - dt * pre_eleMp[ind] *
								(
								// (- pre_eleEP[ind] * pre_eleEP[ind] * (pre_dNdx[ind][m][0] * pre_dNdx[ind][n][0] + pre_dNdx[ind][m][1] * pre_dNdx[ind][n][1])) // terma2
								// - (- pre_dAdx[ind] * pre_eleEEP[ind] * pre_dNdx[ind][m][1] * pre_dNdx[ind][n][0]) // termadx
								// + (- pre_dAdy[ind] * pre_eleEEP[ind] * pre_dNdx[ind][m][0] * pre_dNdx[ind][n][1]) // termady
								+ (- 3 * vars[0] * vars[0] + 2 * (1 - pre_C1[ind]) * vars[0] + pre_C1[ind] * pre_Nx[ind][m]) * pre_Nx[ind][n]) // termdbl
								) * pre_detJ[ind];
						}
					}
					ind += 1; // incrementing index for extracting pre-calculated variables
				}
			}
		}

		/*Apply Boundary Condition*/
		for (int i = 0; i < nen; i++) {
			int A = bzmesh_process[e].IEN[i];
			if (cpts[A].label == 1) {// domain boundary
				ApplyBoundaryCondition(0, i, 0, EMatrixSolve, EVectorSolve);
			}
		}

		ResidualAssembly(EVectorSolve, bzmesh_process[e].IEN, GR_phi);
		MatrixAssembly(EMatrixSolve, bzmesh_process[e].IEN, GK_phi);
	}

	VecAssemblyBegin(GR_phi);
	MatAssemblyBegin(GK_phi, MAT_FINAL_ASSEMBLY);
}

void NeuronGrowth::BuildLinearSystemProcessNG_syn_tub(const vector<Vertex3D> &cpts) {
    int ind = 0; // Index for pre-calculated variables

    // Loop over Bezier mesh elements assigned to the current process
    for (int e = 0; e < bzmesh_process.size(); ++e) {
        int nen = bzmesh_process[e].IEN.size(); // Number of element nodes

        // Initialize element matrices and vectors for synaptogenesis (syn) and tubulin (tub)
        vector<vector<float>> EMatrixSolve_syn(nen, vector<float>(nen, 0.0f)); // Element stiffness matrix for syn
        vector<float> EVectorSolve_syn(nen, 0.0f);                             // Element load vector for syn

        vector<vector<float>> EMatrixSolve_tub(nen, vector<float>(nen, 0.0f)); // Element stiffness matrix for tub
        vector<float> EVectorSolve_tub(nen, 0.0f);                             // Element load vector for tub

        // Initialize storage for element values
        vector<vector<float>> eleVal_st(5, vector<float>(nen, 0.0f)); // Stores element-specific quantities

        // Pre-compute element values for syn and tub
        for (int i = 0; i < nen; ++i) {
            int idx = bzmesh_process[e].IEN[i];
            eleVal_st[0][i] = phi[idx] - phi_prev[idx]; // Difference in phi (elePhiDiff)
            eleVal_st[1][i] = syn[idx];                // Synaptic concentration (eleSyn)
            eleVal_st[2][i] = CellBoundary(phi[idx], 0.5); // Cell boundary condition for current phi (elePhi)
            eleVal_st[3][i] = CellBoundary(phi_prev[idx], 0.5); // Cell boundary condition for previous phi (elePhiPrev)
            eleVal_st[4][i] = tub[idx];               // Tubulin concentration (eleConct)
        }

        // Gaussian quadrature: Loop over Gauss points for numerical integration
        for (int i = 0; i < Gpt.size(); ++i) {
            for (int j = 0; j < Gpt.size(); ++j) {
                for (int k = 0; k < Gpt.size(); ++k) {

                    // Initialize storage for calculated quantities
                    vector<float> vars_st(13, 0.0f); // Stores variables needed for calculations
                    vars_st[7] = pre_mag_grad_phi0[ind]; // Precomputed magnitude of gradient of phi
                    vars_st[11] = pre_term_source[ind];  // Precomputed source term
					// Description of `vars_st` indices:
					//        [0] elePf       - Element-level Pf
					//        [1] eleS        - Element-level Synaptogenesis
					//        [2] eleP        - Element-level Pressure
					//        [3] dPdx        - Derivative of Pressure w.r.t. x
					//        [4] dPdy        - Derivative of Pressure w.r.t. y
					//        [5] elePprev    - Previous step's Pressure
					//        [6] eleC        - Concentration variable
					//        [7] mag_grad_phi0 - Magnitude of gradient of Phi
					//        [8] term_diff   - Diffusion term
					//        [9] term_alph   - Alpha term
					//        [10] term_beta  - Beta term
					//        [11] term_source - Source term
					//        [12] dPdz       - Derivative of Pressure w.r.t. z

                    // Evaluate element-level quantities for syn and tub
                    ElementEvaluationAll_syn_tub(nen, pre_Nx[ind], pre_dNdx[ind], eleVal_st, vars_st);

                    // Assemble element load vectors and stiffness matrices
                    for (int m = 0; m < nen; ++m) {
                        // Load vector for syn
                        EVectorSolve_syn[m] += (vars_st[1] + kappa * vars_st[0]) * pre_Nx[ind][m] * pre_detJ[ind];

                        // Load vector for tub
                        EVectorSolve_tub[m] += (vars_st[6] * vars_st[2] + dt / 2 * vars_st[11]) * pre_Nx[ind][m] * pre_detJ[ind];

                        // Stiffness matrices
                        for (int n = 0; n < nen; ++n) {
                            // Stiffness matrix for syn (assembled only once if judge_syn == 0)
                            if (judge_syn == 0) {
                                EMatrixSolve_syn[m][n] += (pre_Nx[ind][m] * pre_Nx[ind][n] +
                                                           dt / 2 * Dc * (pre_dNdx[ind][m][0] * pre_dNdx[ind][n][0] +
                                                                          pre_dNdx[ind][m][1] * pre_dNdx[ind][n][1] +
                                                                          pre_dNdx[ind][m][2] * pre_dNdx[ind][n][2])) *
                                                          pre_detJ[ind];
                            }

                            // Stiffness matrix for tub
                            EMatrixSolve_tub[m][n] += (pre_Nx[ind][m] * vars_st[2] * pre_Nx[ind][n] -
                                                       dt / 2 * (-Diff * (vars_st[2] * (pre_dNdx[ind][m][0] * pre_dNdx[ind][n][0] +
                                                                                       pre_dNdx[ind][m][1] * pre_dNdx[ind][n][1] +
                                                                                       pre_dNdx[ind][m][2] * pre_dNdx[ind][n][2])) -
                                                                 alphaT * (vars_st[2] * (pre_dNdx[ind][m][0] + pre_dNdx[ind][m][1] + pre_dNdx[ind][m][2]) +
                                                                           (vars_st[3] + vars_st[4] + vars_st[12]) * pre_Nx[ind][m]) *
                                                                     pre_Nx[ind][n] -
                                                                 betaT * (vars_st[2] * pre_Nx[ind][m]) * pre_Nx[ind][n])) *
                                                      pre_detJ[ind];
                        }
                    }
                    ind++; // Move to the next Gauss point
                }
            }
        }

        // Apply boundary conditions
        for (int i = 0; i < nen; ++i) {
            int A = bzmesh_process[e].IEN[i];
            if (cpts[A].label == 1) { // Boundary condition for domain boundary
                ApplyBoundaryCondition(0, i, 0, EMatrixSolve_syn, EVectorSolve_syn);
            }
            if (CellBoundary(phi[A], 0.5) != 1) { // Boundary condition for tubulin outside soma
                ApplyBoundaryCondition(0, i, 0, EMatrixSolve_tub, EVectorSolve_tub);
            }
        }

        // Assemble global residuals and stiffness matrices
        ResidualAssembly(EVectorSolve_syn, bzmesh_process[e].IEN, GR_syn);
        if (judge_syn == 0) { // Matrix remains unchanged, assemble only once
            MatrixAssembly(EMatrixSolve_syn, bzmesh_process[e].IEN, GK_syn);
        }

        ResidualAssembly(EVectorSolve_tub, bzmesh_process[e].IEN, GR_tub);
        MatrixAssembly(EMatrixSolve_tub, bzmesh_process[e].IEN, GK_tub);
    }

    // Finalize assembly
    VecAssemblyBegin(GR_syn);
    if (judge_syn == 0) MatAssemblyBegin(GK_syn, MAT_FINAL_ASSEMBLY);

    VecAssemblyBegin(GR_tub);
    MatAssemblyBegin(GK_tub, MAT_FINAL_ASSEMBLY);
}

void NeuronGrowth::HandleExpansion(const vector<float>& phi_in,
    int& NX, int& NY, int& NZ,
    int& originX, int& originY, int& originZ) 
{
    int expd_dir_global = 6; // Default: No expansion
    // Determine expansion direction on rank 0
    if (comRank == 0) {
        expd_dir_global = CheckExpansion3D(phi_in, cpts, originX, originY, originZ);
	}

    // Broadcast expansion direction to all ranks
    MPI_Bcast(&expd_dir_global, 1, MPI_INT, 0, PETSC_COMM_WORLD);

    // Print detailed expansion information
    PetscPrintf(PETSC_COMM_WORLD, "\n");
    PetscPrintf(PETSC_COMM_WORLD, "============================================================\n");
    PetscPrintf(PETSC_COMM_WORLD, "Expanding in direction: %d\n", expd_dir_global);

    // Define the offset for expansion
    float offset = expand_sz * 4.0; // Half of expansion * delta element (original coarse level)

    // Apply expansion based on global direction
    switch (expd_dir_global) {
        case 0: 
            NX += expand_sz;
            PetscPrintf(PETSC_COMM_WORLD, "Expanding in +X direction. Updated NX: %d\n", NX);
            break; // Positive x boundary
        case 1: 
            NX += expand_sz;
            originX -= offset;
            PetscPrintf(PETSC_COMM_WORLD, "Expanding in -X direction. Updated NX: %d, OriginX: %.2f\n", NX, originX);
            break; // Negative x boundary
        case 2: 
            NY += expand_sz;
            PetscPrintf(PETSC_COMM_WORLD, "Expanding in +Y direction. Updated NY: %d\n", NY);
            break; // Positive y boundary
        case 3: 
            NY += expand_sz;
            originY -= offset;
            PetscPrintf(PETSC_COMM_WORLD, "Expanding in -Y direction. Updated NY: %d, OriginY: %.2f\n", NY, originY);
            break; // Negative y boundary
        case 4: 
            NZ += expand_sz;
            PetscPrintf(PETSC_COMM_WORLD, "Expanding in +Z direction. Updated NZ: %d\n", NZ);
            break; // Positive z boundary
        case 5: 
            NZ += expand_sz;
            originZ -= offset;
            PetscPrintf(PETSC_COMM_WORLD, "Expanding in -Z direction. Updated NZ: %d, OriginZ: %.2f\n", NZ, originZ);
            break; // Negative z boundary
        case 6: 
            PetscPrintf(PETSC_COMM_WORLD, "No expansion required. Global direction: %d\n", expd_dir_global);
            break; // No expansion
        case 7: 
            // Expand all directions
            NX += expand_sz;
            NY += expand_sz;
            NZ += expand_sz;
            originX -= expand_sz;
            originY -= expand_sz;
            originZ -= expand_sz;
            PetscPrintf(PETSC_COMM_WORLD, 
                        "Expanding all directions. Updated NX: %d, NY: %d, NZ: %d, "
                        "OriginX: %.2f, OriginY: %.2f, OriginZ: %.2f\n", 
                        NX, NY, NZ, originX, originY, originZ);
            break; // Expand all directions
        default: 
            PetscPrintf(PETSC_COMM_WORLD, "Invalid expansion direction: %d\n", expd_dir_global);
            break;
    }
}

int NeuronGrowth::CheckExpansion3D(const vector<float>& input,
                                   const vector<Vertex3D>& cpts, 
                                   const int& originX, const int& originY, const int& originZ) 
{
    constexpr float bc_clearance = 2.0f; // Distance from boundary to trigger expansion
    constexpr float epsilon = 1e-5f;    // Small value to account for floating-point precision

    // Debugging: Print grid bounds
    // CheckVar("CheckExp", cpts, input);
    PetscPrintf(PETSC_COMM_WORLD, "Grid bounds: max_x: %.2f, max_y: %.2f, max_z: %.2f\n", 
                max_x, max_y, max_z);

    // Iterate over all control points to detect boundary conditions
    for (size_t i = 0; i < cpts.size(); ++i) {
        // Only consider points where phi exceeds the threshold
        if (input[i] > 0.01f) {
            float currX = cpts[i].coor[0];
            float currY = cpts[i].coor[1];
            float currZ = cpts[i].coor[2];

            // // Debugging: Print each control point and phi value
            // PetscPrintf(PETSC_COMM_WORLD, "Point %zu: (%.2f, %.2f, %.2f), Phi: %.5f\n", 
            //             i, currX, currY, currZ, input[i]);

            // Check each boundary direction in axis order (x, y, z)
            if (currX >= (max_x - bc_clearance - epsilon)) return 0; // Positive x boundary
            if (currX <= (originX + bc_clearance + epsilon)) return 1; // Negative x boundary

            if (currY >= (max_y - bc_clearance - epsilon)) return 2; // Positive y boundary
            if (currY <= (originY + bc_clearance + epsilon)) return 3; // Negative y boundary

            if (currZ >= (max_z - bc_clearance - epsilon)) return 4; // Positive z boundary
            if (currZ <= (originZ + bc_clearance + epsilon)) return 5; // Negative z boundary
        }
    }

    // No boundary condition met
    PetscPrintf(PETSC_COMM_WORLD, "No expansion required.\n");
    return 6; // No expansion
}

void NeuronGrowth::PopulateRandom(vector<float> &input) {
    for (float &value : input) {
        if (value == 0.0f) {
            // Generate a random float between 0 and 1
            value = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
        }
    }
}

bool NeuronGrowth::KD_SearchPair(const vector<Vertex3D>& cpts, 
								const KDTree& kdTree, 
								float targetX, float targetY, float targetZ, 
								int& ind, 
								float tolerance) 
{
    // Ensure the input cloud is not empty
    if (cpts.empty()) {
        PetscPrintf(PETSC_COMM_WORLD, "Error: cpts is empty.\n");
        return false;
    }

    // Query point
    const float query_pt[3] = {targetX, targetY, targetZ};
    size_t closestIndex;
    float out_dist_sqr;

    // Set up KNN result set for 1 neighbor
    nanoflann::KNNResultSet<float> resultSet(1);
    resultSet.init(&closestIndex, &out_dist_sqr);

    // Perform the nearest-neighbor search
    if (!kdTree.findNeighbors(resultSet, query_pt, nanoflann::SearchParameters(0))) {
        PetscPrintf(PETSC_COMM_WORLD, "Error: KDTree search failed.\n");
        return false;
    }

    // Retrieve the coordinates of the closest point
    const float x = cpts[closestIndex].coor[0];
    const float y = cpts[closestIndex].coor[1];
    const float z = cpts[closestIndex].coor[2];

    // Check if the closest point is within the acceptable distance
    if (max({abs(x - targetX), abs(y - targetY), abs(z - targetZ)}) <= tolerance) {
        ind = static_cast<int>(closestIndex);
        return true; // Valid neighbor found
    }

    // No valid neighbor within the tolerance
    return false;
}

float NeuronGrowth::RmOutlier(vector<float> &data) {
	if (data.empty()) {
        return 0.0f; // Handle empty input
    }

    // 1. Filter out zeros and create a new vector
    std::vector<float> non_zero_data;
    for (float value : data) {
        if (value != 0.0f) {
            non_zero_data.push_back(value);
        }
    }

    if (non_zero_data.empty()) {
        return 0.0f; // Or handle the case where all values were zero.
    }


    // 2. Calculate mean and standard deviation using non-zero data
    float sum = std::accumulate(non_zero_data.begin(), non_zero_data.end(), 0.0f);
    float mean = sum / non_zero_data.size();

    float sq_sum = std::inner_product(non_zero_data.begin(), non_zero_data.end(), non_zero_data.begin(), 0.0f,
                               [](float acc, float val) { return acc + val; },
                               [mean](float a, float b) { return std::pow(a - mean, 2) + b; });
    float standardDeviation = std::sqrt(sq_sum / non_zero_data.size());

    // 3. Define the threshold
    float threshold = mean + 3 * standardDeviation;

    // 4. Clamp values in the ORIGINAL data vector (not the filtered one)
    std::transform(data.begin(), data.end(), data.begin(), [threshold](float value) {
        if (value != 0.0f) { // Only clamp non-zero values
          return std::min(value, threshold);
        }
        return value; // Leave zeros unchanged
    });

    // 5. Return the adjusted threshold (based on non-zero data)
    return mean + 2 * standardDeviation;
}

float NeuronGrowth::CellBoundary(float phi, float threshold) {
    return (phi > threshold) ? 1.0f : 0.0f;
}

bool NeuronGrowth::IsInBox(const Vertex3D& point, const Vertex3D& center, float dx, float dy, float dz) {
    if (point.coor[0] < (center.coor[0] - dx/2) || point.coor[0] > (center.coor[0] + dx/2)) return false;
    if (point.coor[1] < (center.coor[1] - dy/2) || point.coor[1] > (center.coor[1] + dy/2)) return false;
    if (point.coor[2] < (center.coor[2] - dz/2) || point.coor[2] > (center.coor[2] + dz/2)) return false;
    return true;
}

bool NeuronGrowth::IsWithinRadius(const Vertex3D& point, const Vertex3D& center, float radius) {
    float dx = point.coor[0] - center.coor[0];
    float dy = point.coor[1] - center.coor[1];
    float dz = point.coor[2] - center.coor[2];
    float distanceSquared = dx * dx + dy * dy + dz * dz;

    return distanceSquared <= (radius * radius); // Compare squared distances to avoid sqrt
}

vector<int> NeuronGrowth::GetBoxNeighbors(
    int idx,
    const vector<Vertex3D> &cpts_fine,
    float dx, float dy, float dz)
{
    vector<int> neighbors;
    const Vertex3D &center = cpts_fine[idx];

    // Naive approach: check every point to see if it's "in box"
    for (int j = 0; j < (int)cpts_fine.size(); j++) {
        if (j == idx) continue; // skip self
        if (IsInBox(cpts_fine[j], center, dx, dy, dz)) {
            neighbors.push_back(j);
        }
    }
    return neighbors;
}

void NeuronGrowth::FindLocalMaximaClusters_box(
    vector<float> &tips_fine,
    const vector<Vertex3D> &cpts_fine,
    float dx, float dy, float dz)
{
    const size_t nPoints = cpts_fine.size();
    if (nPoints == 0) return;

    // Keep track of which points have been visited
    vector<bool> visited(nPoints, false);

    // We'll define a comparison lambda for picking
    // the local maximum based on tips_fine values:
    auto compareTips = [&](int a, int b) {
        return (tips_fine[a] < tips_fine[b]);
    };

	int numClusters(0);
	
    // Loop over each point to find clusters
    for (size_t startIdx = 0; startIdx < nPoints; ++startIdx) {
        // Only proceed if tips_fine[startIdx] != 0 (part of a cluster)
        // and we haven't already visited it.
        if (tips_fine[startIdx] != 0.0f && !visited[startIdx]) {
			numClusters++;
            // We'll gather all points in this cluster
            queue<int> Q;
            vector<int> clusterIndices;

            // Start BFS
            Q.push((int)startIdx);
            visited[startIdx] = true;
            clusterIndices.push_back((int)startIdx);

            while (!Q.empty()) {
                int current = Q.front();
                Q.pop();

                // 1) Get neighbors within the bounding box around cpts_fine[current]
                auto neighbors = GetBoxNeighbors(current, cpts_fine, dx, dy, dz);

                // 2) For each neighbor, if it's part of tips_fine (non-zero) and unvisited,
                //    add it to BFS.
                for (int nb : neighbors) {
                    if (tips_fine[nb] != 0.0f && !visited[nb]) {
                        visited[nb] = true;
                        Q.push(nb);
                        clusterIndices.push_back(nb);
                    }
                }
            }

            // Now clusterIndices contains all points of this connected component
            // => find the local maximum by tips_fine value
            auto maxIt = max_element(clusterIndices.begin(),
                                          clusterIndices.end(),
                                          compareTips);
            int localMaxIdx = *maxIt;  // index with highest tips_fine value

            // => label only that local max as 1, set all other indices to 0
            for (int idx : clusterIndices) {
                tips_fine[idx] = 0.0f; 
            }
            tips_fine[localMaxIdx] = 1.0f; 
        }
    }

    // Print the number of clusters using PetscPrintf
    PetscPrintf(PETSC_COMM_WORLD, "Number of clusters identified: %d\n", numClusters);
}

void NeuronGrowth::DetectTips(const vector<Vertex3D>& cpts_fine, 
							const Vertex3DCloud& cloud_fine,
							const KDTree& kdTree_fine,
							const float& tip_I_sz,
							const vector<Vertex3D>& cpts,
							const Vertex3DCloud& cloud,
							const KDTree& kdTree,
							vector<array<float, 3>>& seed,
							const int NX, const int NY, const int NZ,
							const int originX, const int originY, const int originZ
)
{
	vector<float> phi_fine(cpts_fine.size(), 0.0f);

	for (size_t i = 0; i < cpts_fine.size(); ++i) {
		InterpolateOrFindExact_singleVar(
			cpts_fine[i], kdTree, cloud, phi, cpts, phi_fine[i], true);
	}

	// if (n % var_save_invl == 0) CheckVar("PHI_FINE", cpts_fine, phi_fine);
	// CheckVar("PHI_FINE", cpts_fine, phi_fine);

    const float threshold = 0.95f;    // Threshold for tip detection
    float maxTipValue = 0.0f;        // Tracks maximum tip value for normalization

    // Precompute transformed phi values
    vector<float> phiTransformed(phi_fine.size());
    for (size_t j = 0; j < phi_fine.size(); ++j) {
        phiTransformed[j] = CellBoundary(phi_fine[j], 0.50f);
    }

    // Clear and resize tips to match the number of control points
    vector<float> tips_fine(cpts_fine.size(), 0.0f);
	
    // compute tip intensity
    for (size_t id = 0; id < seed.size(); ++id) {
		for (size_t i = 0; i < cpts_fine.size(); ++i) {
			const auto& center = cpts_fine[i];
			float localSum = 0.0f; // Sum of phi values within the box

			// Compute the sum of phi values for points within the vicinity
			for (size_t j = 0; j < phi_fine.size(); ++j) {
				// if (IsInBox(cpts_fine[j], center, tip_I_sz, tip_I_sz, tip_I_sz) && neurons_flattern[i] == id) {
				if (IsInBox(cpts_fine[j], center, tip_I_sz, tip_I_sz, tip_I_sz)) {
				// if (IsWithinRadius(cpts_fine[j], center, tip_I_sz)) {
					localSum += phiTransformed[j];
				}
			}

			// Compute the tip score for the current vertex
			// if (localSum > 0.0f && tmp[i] < INF) {
			if (localSum > 0.0f) {
				tips_fine[i] = (phiTransformed[i] / localSum) * phiTransformed[i];
				//  * tmp[i];
			} else {
				tips_fine[i] = 0.0f; // Avoid division by zero
			}
			maxTipValue = max(maxTipValue, tips_fine[i]);
		}
	}

    // CheckVar("TIP_FINE_", cpts_fine, tips_fine);

	maxTipValue = min(maxTipValue, 0.00130f);
	// maxTipValue = 0.00133;

    // Thresholding and normalization
    for (float& tip : tips_fine) {
        tip = (tip > threshold * maxTipValue) ? 1.0f : 0.0f;
    }
	// CheckVar("TIP_FINE_cutoff_", cpts_fine, tips_fine);

	FindLocalMaximaClusters_box(tips_fine, cpts_fine, 4, 4, 4);

    // Debugging and visualization
    // if (n % var_save_invl == 0) CheckVar("TIP_FINE_", cpts_fine, tips_fine);
    // CheckVar("TIP_local_FINE_", cpts_fine, tips_fine);

    // Clear and resize tips to match the number of control points
    tips.clear();
    tips.resize(cpts.size(), 0.0f);

	for (size_t i = 0; i < cpts.size(); ++i) {
		InterpolateOrFindExact_singleVar(
			cpts[i], kdTree_fine, cloud_fine, tips_fine, cpts_fine, tips[i], false);
	}
	// CheckVar("TIP_FINAL", cpts, tips);

}

void NeuronGrowth::DetectTips_multi(const vector<Vertex3D>& cpts_fine, 
	const Vertex3DCloud& cloud_fine,
	const KDTree& kdTree_fine,
	const float& tip_I_sz,
	const vector<Vertex3D>& cpts,
	const Vertex3DCloud& cloud,
	const KDTree& kdTree,
	vector<array<float, 3>>& seed,
	const int NX, const int NY, const int NZ,
	const int originX, const int originY, const int originZ
) {
	vector<float> phi_fine(cpts_fine.size(), 0.0f);

	for (size_t i = 0; i < cpts_fine.size(); ++i) {
		InterpolateOrFindExact_singleVar(
		cpts_fine[i], kdTree, cloud, phi, cpts, phi_fine[i], true);
	}

	if (n % var_save_invl == 0) CheckVar("PHI_FINE", cpts_fine, phi_fine);

	const float threshold = 0.95f;    // Threshold for tip detection
	float maxTipValue = 0.0f;        // Tracks maximum tip value for normalization

	// Precompute transformed phi values
	vector<float> phiTransformed(phi_fine.size());
	for (size_t j = 0; j < phi_fine.size(); ++j) {
		phiTransformed[j] = CellBoundary(phi_fine[j], 0.50f);
	}

	// Clear and resize tips to match the number of control points
	vector<float> tips_fine(cpts_fine.size(), 0.0f);

	// Initialize neurons

	vector<vector<vector<int>>> neurons;
	IdentifyNeurons3DWithKDTree(phi_fine, neurons, seed, NX, NY, NZ, originX, originY, originZ, kdTree_fine, cloud_fine);
	auto neurons_flattern = Convert3DIntTo1DFloatVector(neurons);
	CheckVar("NEURONS_", cpts_fine, neurons_flattern);
	// vector<vector<vector<int>>> geoDist = CalculateGeodesicDistanceFromPoint3D(neurons, seed, originX, originY, originZ);
	// auto geoDist_flattern = Convert3DIntTo1DFloatVector(geoDist);
	// CheckVar("DIST_", cpts_fine, geoDist_flattern);
	
	// Clear and resize tips to match the number of control points
	tips.clear();
	tips.resize(cpts.size(), 0.0f);	
	// compute tip intensity
	for (size_t id = 0; id < seed.size(); ++id) {
		vector<float> tips_fine(cpts_fine.size(), 0.0f);
		for (size_t i = 0; i < cpts_fine.size(); ++i) {
			const auto& center = cpts_fine[i];
			float localSum = 0.0f; // Sum of phi values within the box

			// Compute the sum of phi values for points within the vicinity
			for (size_t j = 0; j < phi_fine.size(); ++j) {
				if (IsInBox(cpts_fine[j], center, tip_I_sz, tip_I_sz, tip_I_sz) && neurons_flattern[j] == neurons_flattern[i]) {
				// if (IsInBox(cpts_fine[j], center, tip_I_sz, tip_I_sz, tip_I_sz)) {
					localSum += phiTransformed[j];
				}
			}

			// Compute the tip score for the current vertex
			if (localSum > 0.0f) {
				tips_fine[i] = (phiTransformed[i] / localSum) * phiTransformed[i];
				//  * tmp[i];
			} else {
				tips_fine[i] = 0.0f; // Avoid division by zero
			}
		}

		vector<float> tips_tmp(cpts.size(), 0.0f);
		for (size_t i = 0; i < cpts.size(); ++i) {
			InterpolateOrFindExact_singleVar(
			cpts[i], kdTree_fine, cloud_fine, tips_fine, cpts_fine, tips_tmp[i], false);
			maxTipValue = max(maxTipValue, tips[i]);
		}
		// cout << maxTipValue << endl;
		maxTipValue = 0.0065;
		// if (n % var_save_invl == 0) CheckVar("TIP_", cpts, tips);

		// Thresholding and normalization
		for (float& tip : tips_tmp) {
			tip = (tip > threshold * maxTipValue) ? 1.0f : 0.0f;
		}
		for (int i = 0; i < tips_tmp.size(); i++) {
			tips_tmp[i] = (tips_tmp[i] > threshold * maxTipValue) ? 1.0f : 0.0f;
			tips[i] += tips_tmp[i];
		}
		// if (n % var_save_invl == 0) CheckVar("TIP_cutoff_", cpts, tips);
	}
	if (n % var_save_invl == 0) CheckVar("TIP_cutoff_", cpts, tips);

}

vector<pair<Vertex3D, int>> NeuronGrowth::FindClosestVerticesWithIndices(const vector<Vertex3D>& vertices, const Vertex3D& inputVertex, int k) {
	// Custom comparator that prioritizes larger squared distances and considers the vertex index
	auto comp = [&inputVertex](const pair<Vertex3D, int>& a, const pair<Vertex3D, int>& b) {
		return SquaredDistance(inputVertex, a.first) < SquaredDistance(inputVertex, b.first);
	};

	// Initialize the priority queue with the custom comparator
	priority_queue<pair<Vertex3D, int>, vector<pair<Vertex3D, int>>, decltype(comp)> pq(comp);

	for (int i = 0; i < vertices.size(); ++i) {
		pq.emplace(vertices[i], i);
		if (pq.size() > k) pq.pop(); // Remove the vertex with the largest distance
	}

	vector<pair<Vertex3D, int>> closestVertices;
	while (!pq.empty()) {
		closestVertices.push_back(pq.top()); // Collect the closest vertices and their indices
		pq.pop();
	}
	// Reverse to correct the order from farthest of the closest to nearest
	reverse(closestVertices.begin(), closestVertices.end());

	return closestVertices;
}

vector<tuple<Vertex3D, int, float>> NeuronGrowth::FindClosestVerticesWithIndicesAndDistances(
    const KDTree& kdTree, const Vertex3DCloud& cloud, const Vertex3D& inputVertex, int k) 
{
    // Ensure k does not exceed the number of points in the KDTree
    k = min(k, static_cast<int>(cloud.pts.size()));

    vector<size_t> closestIndices(k);
    vector<float> squaredDistances(k);

    // Query point: extract coordinates from inputVertex
    float queryPoint[3] = {inputVertex.coor[0], inputVertex.coor[1], inputVertex.coor[2]};

    // Initialize result set for k nearest neighbors
    nanoflann::KNNResultSet<float> resultSet(k);
    resultSet.init(closestIndices.data(), squaredDistances.data());

    // Perform the nearest neighbor search
    kdTree.findNeighbors(resultSet, queryPoint, nanoflann::SearchParameters());

    // Collect results: vertices, indices, and distances
    vector<tuple<Vertex3D, int, float>> closestVertices;
    for (size_t i = 0; i < k; ++i) {
        const Vertex3D& vertex = cloud.pts[closestIndices[i]];
        float distance = sqrt(squaredDistances[i]); // Convert squared distance to actual distance
        // float distance = squaredDistances[i]; // remove sqrt to improve computational efficiency
        closestVertices.emplace_back(vertex, static_cast<int>(closestIndices[i]), distance);
    }

    return closestVertices;
}

vector<float> NeuronGrowth::ComputeRefine(
    const vector<float>& phi_in, 
	const int& NX, const int& NY, const int& NZ,
	const int& originX, const int& originY, const int& originZ,
    const KDTree& kdTree, const Vertex3DCloud& cloud) 
{
	// CheckVar("PHI_", cpts, phi_in); // check control points phi for debugging

    // Initialize the refined elements vector
    vector<float> ele_refine(NX * NY * NZ, 0.0);
	const float phi_ceil = 0.95f, phi_floor = 0.0f;

    // Loop through the 3D grid to compute the refinement flags
	// (some limitations in THS3D when it comes to locally refining boundary elements - see Xiaodong's code)
    for (int i = 1; i < NX - 1; ++i) {       // Avoid boundaries in x
        for (int j = 1; j < NY - 1; ++j) {   // Avoid boundaries in y
            for (int k = 1; k < NZ - 1; ++k) { // Avoid boundaries in z

                // Define the current query point
                Vertex3D queryPoint;
				float spacing = 4.0f;
				queryPoint.coor[0] = i * spacing + originX;
				queryPoint.coor[1] = j * spacing + originY;
				queryPoint.coor[2] = k * spacing + originZ;

                // Find the k closest vertices and their distances
                int numNeighbors = 6; // Number of neighbors to consider
                vector<tuple<Vertex3D, int, float>> closestVertices = 
                    FindClosestVerticesWithIndicesAndDistances(kdTree, cloud, queryPoint, numNeighbors);
									
				// float phi_average(0);
				float phi_max = 0.0f; // Initialize phi_max to a default value
				if (!closestVertices.empty()) {
					// // Compute weighted average for phi
					// float weightedSum = 0.0;
					// float weightTotal = 0.0;
					// const float epsilon = 1e-6f; // Avoid division by zero

					// for (const auto& neighbor : closestVertices) {
					// 	int idx = get<1>(neighbor);     // Index of the neighbor
					// 	float distance = get<2>(neighbor); // Distance to the neighbor
					// 	float weight = 1.0f / (distance + epsilon); // Weight based on distance
					// 	weightedSum += CellBoundary(phi_in[idx], 0.1) * weight;
					// 	weightTotal += weight;
					// }

					// phi_average = weightedSum / weightTotal; // Compute weighted average

					// Iterate over the closest vertices to find the maximum phi value
					for (const auto& neighbor : closestVertices) {
						int idx = get<1>(neighbor);         // Index of the neighbor
						float phi_value = CellBoundary(phi_in[idx], 0.1); // Adjust phi value using CellBoundary
						phi_max = max(phi_max, phi_value); // Update phi_max if current phi_value is greater
					}
				}

                // Compute indices for output grid (ele_refine)
				int index_out = k * NX * NY + j * NX + i;

                // Apply refinement criteria based on phi thresholds
                const float phi_detected_threshold = 0.001f;
				if (phi_max > phi_detected_threshold) {
					ele_refine[index_out] = 1.0f;
				} else {
                    ele_refine[index_out] = 0.0f; // No refinement
                    // ele_refine[index_out] = 1.0f; // for debugging
                }
            }
        }
    }
    return ele_refine;
}

// Function to perform Breadth-First Search (BFS) for clustering in 3D
void NeuronGrowth::BFS3D(const vector<float>& matrix, int depth, int rows, int cols, int dep, int row, int col,
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
					BFS3D(matrix, depth, rows, cols, d, i, j, visited, cluster);

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

// void NeuronGrowth::FloodFill3DWithKDTree(vector<vector<vector<int>>>& image,
//                                          int x, int y, int z, int newColor, int originalColor,
//                                          const KDTree& kdTree, const Vertex3DCloud& cloud) 
// {
//     if (x < 0 || x >= image.size() || 
//         y < 0 || y >= image[0].size() || 
//         z < 0 || z >= image[0][0].size() || 
//         image[x][y][z] != originalColor || 
//         image[x][y][z] == newColor) {
//         return;
//     }

//     image[x][y][z] = newColor;

//     // Search for neighbors using KD_SearchPair in six possible directions
//     int dx[] = {0, 0, -1, 1, 0, 0};
//     int dy[] = {-1, 1, 0, 0, 0, 0};
//     int dz[] = {0, 0, 0, 0, -1, 1};

//     for (int i = 0; i < 6; ++i) {
//         int nx = x + dx[i];
//         int ny = y + dy[i];
//         int nz = z + dz[i];

//         // Validate bounds
//         if (nx >= 0 && nx < image.size() && ny >= 0 && ny < image[0].size() && nz >= 0 && nz < image[0][0].size()) {
//             int index;
//             if (KD_SearchPair(cloud.pts, kdTree, nx, ny, nz, index)) {
//                 FloodFill3DWithKDTree(image, nx, ny, nz, newColor, originalColor, kdTree, cloud);
//             }
//         }
//     }
// }

void NeuronGrowth::FloodFill3DWithKDTree(std::vector<std::vector<std::vector<int>>>& image,
	int x, int y, int z, float newColor, float originalColor,
	const KDTree& kdTree, const Vertex3DCloud& cloud)
{
	std::stack<std::array<int, 3>> stack;
	stack.push({x, y, z});

	// Get the correct dimensions.
	int imageNX = image.size();
	int imageNY = image[0].size();
	int imageNZ = image[0][0].size();

	while (!stack.empty()) {
		auto [cx, cy, cz] = stack.top();
		stack.pop();

		// Correct bounds checking here!
		if (cx < 0 || cx >= imageNX ||
		cy < 0 || cy >= imageNY ||
		cz < 0 || cz >= imageNZ ||
		image[cx][cy][cz] != 1) {
			continue;
		}

		image[cx][cy][cz] = newColor;

		int dx[] = {0, 0, -1, 1, 0, 0};
		int dy[] = {-1, 1, 0, 0, 0, 0};
		int dz[] = {0, 0, 0, 0, -1, 1};

		for (int i = 0; i < 6; ++i) {
			int nx = cx + dx[i];
			int ny = cy + dy[i];
			int nz = cz + dz[i];

			if (nx >= 0 && nx < imageNX && ny >= 0 && ny < imageNY && nz >= 0 && nz < imageNZ) {
				stack.push({nx, ny, nz});
			}
		}
	}
}

void NeuronGrowth::IdentifyNeurons3DWithKDTree(vector<float> phi_fine,
                                               vector<vector<vector<int>>>& neurons, 
                                               const vector<array<float, 3>>& seed,
                                               int NX, int NY, int NZ, 
                                               int originX, int originY, int originZ,
                                               const KDTree& kdTree, const Vertex3DCloud& cloud) 
{
    neurons = ConvertTo3DIntVector(phi_fine, NX * 4 + 1, NY * 4 + 1, NZ * 4 + 1);
	// PetscPrintf(PETSC_COMM_WORLD, "phi_fine %d, %d %d %d\n", phi_fine.size(), NX, NY, NZ);

    for (size_t i = 0; i < seed.size(); ++i) {
        int startZ = static_cast<int>((seed[i][0] - originX));
        int startY = static_cast<int>((seed[i][1] - originY));
        int startX = static_cast<int>((seed[i][2] - originZ));

        PetscPrintf(PETSC_COMM_WORLD, "Cluster %zu:\n", i + 1);
        PetscPrintf(PETSC_COMM_WORLD, "Seed: (%f, %f, %f)\n", seed[i][0], seed[i][1], seed[i][2]);
        PetscPrintf(PETSC_COMM_WORLD, "Origin: (%d, %d, %d)\n", originX, originY, originZ);
        PetscPrintf(PETSC_COMM_WORLD, "Start: (%d, %d, %d)\n", startX, startY, startZ);

        if (startX < 0 || startX >= neurons.size() || 
            startY < 0 || startY >= neurons[0].size() || 
            startZ < 0 || startZ >= neurons[0][0].size()) {
            continue;
        }

        float newColor = static_cast<float>(i + 2);
        float originalColor = static_cast<float>(neurons[startX][startY][startZ]);
        // PetscPrintf(PETSC_COMM_WORLD, "Colors - neuron: %f, Original: %f, New: %f\n", neurons[startX][startY][startZ], originalColor, newColor);

        if (originalColor != newColor) {
            FloodFill3DWithKDTree(neurons, startX, startY, startZ, newColor, originalColor, kdTree, cloud);
        }
    }
}

bool NeuronGrowth::IsValid(const int& x, const int& y, const int& z, 
						const int& rows, const int& cols, 
						const int& depth) 
{
    return x >= 0 && x < rows && y >= 0 && y < cols && z >= 0 && z < depth;
}

vector<vector<vector<int>>> NeuronGrowth::CalculateGeodesicDistanceFromPoint3D(
    vector<vector<vector<int>>> neurons, const vector<array<float, 3>>& seed,
    int originX, int originY, int originZ) 
{
    int rows = neurons.size();
    int cols = (rows > 0) ? neurons[0].size() : 0;
    int depth = (cols > 0) ? neurons[0][0].size() : 0;

    // Initialize distances with INF and visited flags as false.
    vector<vector<vector<int>>> distances(rows, vector<vector<int>>(cols, vector<int>(depth, INF)));
    vector<vector<vector<bool>>> visited(rows, vector<vector<bool>>(cols, vector<bool>(depth, false)));

    // 6 possible movements in 3D (up, down, left, right, front, back).
    static const int dx[] = {-1, 1, 0, 0, 0, 0};
    static const int dy[] = {0, 0, -1, 1, 0, 0};
    static const int dz[] = {0, 0, 0, 0, -1, 1};

    // Perform BFS for each seed point.
    for (const auto& point : seed) {
        int startX = point[0] - originX;
        int startY = point[1] - originY;
        int startZ = point[2] - originZ;

        // Ensure the starting point is valid and within bounds.
        if (!IsValid(startX, startY, startZ, rows, cols, depth) || neurons[startX][startY][startZ] == 0) {
            continue;
        }

        distances[startX][startY][startZ] = 0;
        visited[startX][startY][startZ] = true;

        queue<array<int, 3>> q;
        q.push({startX, startY, startZ});

        // BFS traversal to calculate distances.
        while (!q.empty()) {
            auto [x, y, z] = q.front();
            q.pop();

            for (int i = 0; i < 6; ++i) {
                int newX = x + dx[i];
                int newY = y + dy[i];
                int newZ = z + dz[i];

                // Check bounds, neuron presence, and visitation status.
                if (IsValid(newX, newY, newZ, rows, cols, depth) && neurons[newX][newY][newZ] == 1 && !visited[newX][newY][newZ]) {
                    visited[newX][newY][newZ] = true;
                    distances[newX][newY][newZ] = distances[x][y][z] + 1;
                    q.push({newX, newY, newZ});
                }
            }
        }
    }
    return distances;
}

void NeuronGrowth::SaveNGvars(const vector<vector<float>>& NGvars, int NX, int NY, const string& fn) {
	bool visualization = true; // Flag to enable/disable visualization

	// Array of variable names corresponding to NGvars indices
	const char* varNames[] = {"phi", "syn", "tub", "theta", "phi_0", "tub_0"};

	// Loop through NGvars and save each to a text file with a corresponding name
	for (size_t i = 0; i < NGvars.size() && i < sizeof(varNames)/sizeof(varNames[0]); ++i) {
		string filepath = fn + "/" + varNames[i] + "_" + to_string(n) + ".txt";
		PrintVec2TXT(NGvars[i], filepath, visualization);
	}
}

void NeuronGrowth::PrintOutNeurons3D(const vector<vector<vector<int>>>& neurons) 
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

void NeuronGrowth::PrintStatus(int n, int end_iter, int reason_phi, int its_phi, double t_phi,
                               int reason_syn, int reason_tub, int its_syn, int its_tub, double t_syn_tub, 
                               int n_bzmesh) {
    ostringstream oss;

    // Use setw for aligned, compact output
    oss << "Step: " << setw(4) << n << "/" << setw(4) << end_iter
        << " | Phi: " << setw(2) << reason_phi << "[" << setw(3) << its_phi << "] " 
        << fixed << setprecision(2) << setw(6) << t_phi << "s"
        << " | Syn: " << setw(2) << reason_syn << "[" << setw(3) << its_syn << "] " 
        << "Tub: " << setw(2) << reason_tub << "[" << setw(3) << its_tub << "] " 
        << fixed << setprecision(2) << setw(6) << t_syn_tub << "s"
        << " | #bzElem: " << setw(4) << n_bzmesh;

    // MPI-safe output
    PetscPrintf(PETSC_COMM_WORLD, "%s\n", oss.str().c_str());
}

PetscErrorCode SetupSNES(SNES &snes, const char *solverType, void *ctx,
                         PetscErrorCode (*formFunction)(SNES, Vec, Vec, void *),
                         PetscErrorCode (*formJacobian)(SNES, Vec, Mat, Mat, void *),
                         PetscReal rtol = 1e-5, PetscReal atol = 1e-7, PetscReal dtol = 1e-9,
                         PetscInt maxIters = 100, PetscInt maxFails = 1000) {
    // Create and set SNES type
    CHKERRQ(SNESCreate(PETSC_COMM_WORLD, &snes));
    CHKERRQ(SNESSetType(snes, solverType));

    // Set solver tolerances
    CHKERRQ(SNESSetTolerances(snes, rtol, atol, dtol, maxIters, maxFails));

    // Configure line search
    SNESLineSearch linesearch;
    CHKERRQ(SNESGetLineSearch(snes, &linesearch));
    CHKERRQ(SNESLineSearchSetType(linesearch, SNESLINESEARCHCP));
    CHKERRQ(SNESLineSearchSetDamping(linesearch, 0.8));

    // Validate function pointers
    if (!formFunction || !formJacobian) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Function or Jacobian is null");
    }

    // Set function and Jacobian
    CHKERRQ(SNESSetFunction(snes, NULL, formFunction, ctx));
    CHKERRQ(SNESSetJacobian(snes, NULL, NULL, formJacobian, ctx));

    return PETSC_SUCCESS; // Explicitly return PETSC_SUCCESS on success
}

PetscErrorCode SetupKSP(KSP &ksp, Mat &A, const char *kspType, const char *pcType,
                        PetscReal rtol = 1.e-8, PetscReal atol = PETSC_DEFAULT,
                        PetscReal dtol = PETSC_DEFAULT, PetscInt maxIters = 100000, PetscInt restart = 100) {
    // Validate matrix
    if (!A) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Matrix A is not initialized");
    }

    // Create and configure KSP solver
    CHKERRQ(KSPCreate(PETSC_COMM_WORLD, &ksp));
    CHKERRQ(KSPSetOperators(ksp, A, A));
    CHKERRQ(KSPSetType(ksp, kspType));

    // Configure preconditioner
    PC pc;
    CHKERRQ(KSPGetPC(ksp, &pc));
    CHKERRQ(PCSetType(pc, pcType));
    // Additional preconditioner options can be set here if needed
    // Example: CHKERRQ(PCSetType(pc, PCFIELDSPLIT)); or CHKERRQ(PCSetType(pc, PCGAMG));

    // Set GMRES-specific options
    if (string(kspType) == "KSPGMRES") {
        CHKERRQ(KSPGMRESSetRestart(ksp, restart));
    }

    // Configure solver tolerances and options
    CHKERRQ(KSPSetInitialGuessNonzero(ksp, PETSC_TRUE));
    CHKERRQ(KSPSetTolerances(ksp, rtol, atol, dtol, maxIters));
    CHKERRQ(KSPSetFromOptions(ksp));
    CHKERRQ(KSPSetUp(ksp));

    return PETSC_SUCCESS; // Explicitly return success
}

PetscErrorCode ScatterVector(Vec src, vector<float>& target, PetscInt size, 
							bool applyBoundary = false, NeuronGrowth* NG = nullptr) {
							PetscErrorCode ierr;
							Vec temp_seq;             // Sequential vector for scattered data
							VecScatter ctx;           // Scatter context
							const PetscScalar* array; // Array pointer for data in `temp_seq`

    // Create scatter context and scatter data from distributed vector to sequential vector
    ierr = VecScatterCreateToAll(src, &ctx, &temp_seq); CHKERRQ(ierr);
    ierr = VecScatterBegin(ctx, src, temp_seq, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterEnd(ctx, src, temp_seq, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);

    // Access the array in the sequential vector
    ierr = VecGetArrayRead(temp_seq, &array); CHKERRQ(ierr);

    // Resize target vector to match the size
    target.resize(size);

    // Populate the target vector
    for (PetscInt i = 0; i < size; ++i) {
        if (applyBoundary && NG) {
            // Apply boundary function and clamp values as needed
            target[i] = PetscMax(PetscRealPart(array[i]), 0.0f) * NG->CellBoundary(NG->phi[i], 0.5);
            target[i] = isnan(target[i]) || target[i] > 1 ? 0 : target[i];
        } else {
            target[i] = PetscRealPart(array[i]);
        }
    }

    // Restore the array and clean up resources
    ierr = VecRestoreArrayRead(temp_seq, &array); CHKERRQ(ierr);
    ierr = VecScatterDestroy(&ctx); CHKERRQ(ierr);
    ierr = VecDestroy(&temp_seq); CHKERRQ(ierr);

    return ierr;
}

PetscErrorCode FormFunction_phi(SNES snes, Vec x, Vec F, void *ctx)
{
    PetscErrorCode ierr;
    NeuronGrowth *user = (NeuronGrowth *)ctx;
    PetscInt e, i, j, k;

    Vec P_seq;
    VecScatter scatter_ctx1;
    PetscScalar *Parray;
    ierr = VecScatterCreateToAll(x, &scatter_ctx1, &P_seq); CHKERRQ(ierr);
    ierr = VecScatterBegin(scatter_ctx1, x, P_seq, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterEnd(scatter_ctx1, x, P_seq, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGetArray(P_seq, &Parray); CHKERRQ(ierr);
    ierr = VecSet(F, 0.0); CHKERRQ(ierr);

    int ind = 0; 
    const size_t nel = user->bzmesh_process.size();
    const int gptSize = (int)user->Gpt.size();

    for (e = 0; e < (int)nel; e++) {
        int nen = (int)user->bzmesh_process[e].IEN.size();

		// Reset EVectorSolve and EMatrixSolve without resizing
		for (int m = 0; m < nen; m++) {
			user->EVectorSolve[m] = 0.0f; // Reset vector values

			// Reset each row of the matrix up to nen
			fill(user->EMatrixSolve[m].begin(), user->EMatrixSolve[m].begin() + nen, 0.0f);
		}
		fill(user->eleVal[0].begin(), user->eleVal[0].end(), 0.0f); // Reset values

        // Extract nodal values just once
        const auto &IEN = user->bzmesh_process[e].IEN;
        for (int ii = 0; ii < nen; ii++) {
            int A = IEN[ii];
            user->pre_eleVal[e][0][ii] = (float)Parray[A];           // phiGuess
        }

        for (i = 0; i < gptSize; i++) {
            for (j = 0; j < gptSize; j++) {
                for (k = 0; k < gptSize; k++) {
					float eleAniso(0), dA_dPdx(0), dA_dPdy(0), dA_dPdz(0);
					if (user->n > 0)
						user->EvaluateOrientation(nen, user->pre_Nx[ind], user->pre_dNdx[ind], user->pre_eleVal[e][1], user->pre_eleVal[e][4], eleAniso, dA_dPdx, dA_dPdy, dA_dPdz);
					user->ElementEvaluationAll_phi(nen, user->pre_Nx[ind], user->pre_dNdx[ind], user->pre_eleVal[e], user->vars);

					float eleMp;
					eleMp = 1;

					// adjust rg (assembly rate) and sg (disassembly rate) based on detected tips
					if (user->n < 0) {
						user->vars[8] = user->alphaOverPi*atan(user->gamma * (1 - user->vars[6]));
					} else {
						if (user->vars[9] > 0) {
							user->vars[8] = user->alphaOverPi*atan(user->gamma * 1 * (1 - user->vars[6]));
						} else {
							user->vars[8] = user->alphaOverPi*atan(user->gamma * 0.01 * (1 - user->vars[6]));
						}
					}

					// calculate C1 variable for phase field energy term
					user->vars[0] = user->vars[8] - user->pre_C0[ind];
					// loop through control points
					for (int m = 0; m < nen; m++) {
						user->EVectorSolve[m] += (user->vars[2] * user->pre_Nx[ind][m] - user->dt * eleMp * (
							(- eleAniso * eleAniso * (user->vars[3] * user->pre_dNdx[ind][m][0] + user->vars[4] * user->pre_dNdx[ind][m][1] + user->vars[18] * user->pre_dNdx[ind][m][2]))
							+ (- user->pre_dNdx[ind][m][0] * eleAniso * dA_dPdx * ( user->vars[3] * user->vars[3] + user->vars[4] * user->vars[4] + user->vars[18] * user->vars[18] ) )
							+ (- user->pre_dNdx[ind][m][1] * eleAniso * dA_dPdy * ( user->vars[3] * user->vars[3] + user->vars[4] * user->vars[4] + user->vars[18] * user->vars[18] ) )
							+ (- user->pre_dNdx[ind][m][2] * eleAniso * dA_dPdz * ( user->vars[3] * user->vars[3] + user->vars[4] * user->vars[4] + user->vars[18] * user->vars[18] ) )
							+ (- user->vars[2] * user->vars[2] * user->vars[2] + (1 - user->vars[0]) * user->vars[2] * user->vars[2] + user->vars[0] * user->vars[2]) * user->pre_Nx[ind][m]
							) - user->vars[5] * user->pre_Nx[ind][m]
							) * user->pre_detJ[ind];

						// loop through 16 control points
						for (int n = 0; n < nen; n++) {
							user->EMatrixSolve[m][n] += (user->pre_Nx[ind][m] * user->pre_Nx[ind][n] - user->dt * eleMp * (
								(- eleAniso * eleAniso * (user->pre_dNdx[ind][m][0] * user->pre_dNdx[ind][n][0] + user->pre_dNdx[ind][m][1] * user->pre_dNdx[ind][n][1] + user->pre_dNdx[ind][m][2] * user->pre_dNdx[ind][n][2])) // terma2
								+ (- user->pre_dNdx[ind][m][0] * eleAniso * dA_dPdx * ( 2 * user->vars[3] + 2 * user->vars[4] + 2 * user->vars[18]))
								+ (- user->pre_dNdx[ind][m][1] * eleAniso * dA_dPdy * ( 2 * user->vars[3] + 2 * user->vars[4] + 2 * user->vars[18]))
								+ (- user->pre_dNdx[ind][m][2] * eleAniso * dA_dPdz * ( 2 * user->vars[3] + 2 * user->vars[4] + 2 * user->vars[18]))
								+ (- 3 * user->vars[2] * user->vars[2] + 2 * (1 - user->vars[0]) * user->vars[2] + user->vars[0] * user->pre_Nx[ind][m]) * user->pre_Nx[ind][n] // termdbl
								)) * user->pre_detJ[ind];

						}
					}
					ind += 1; // incrementing index for extracting pre-calculated variables
                }
            }
        }

        // Apply Boundary Condition
        for (int ii = 0; ii < nen; ii++) {
            int A = user->bzmesh_process[e].IEN[ii];
            if (user->cpts[A].label == 1) {
                user->ApplyBoundaryCondition(0, ii, 0, user->EMatrixSolve, user->EVectorSolve);
            }
        }

        user->ResidualAssembly(user->EVectorSolve, user->bzmesh_process[e].IEN, F);
        user->MatrixAssembly(user->EMatrixSolve, user->bzmesh_process[e].IEN, user->J);
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

/**
 * @brief FormFunction_phi_preComputed:
 *        Computes the residual and Jacobian for the phase-field equation
 *        using arrays that were precomputed in PreparePhaseField_SNES_preComputed().
 *
 * @details
 *  - We gather the current solution guess (phi) from the global vector \a x
 *    into \a P_seq, then copy those values into user->pre_eleVal[e][0].
 *  - We then loop over each element and each Gauss point, computing the
 *    element-level contributions to the residual vector (F) and the global
 *    Jacobian matrix (user->J).
 *  - The pre_eleMp, pre_eleAniso, pre_C1, etc. arrays have been computed
 *    beforehand, so the only new computations here are those strictly
 *    dependent on the updated phi guess.
 *
 * @param[in]  snes  The SNES (nonlinear solver) object
 * @param[in]  x     Current solution vector
 * @param[out] F     Residual vector
 * @param[in]  ctx   Pointer to NeuronGrowth data structure
 *
 * @return PETSc error code
 */
PetscErrorCode FormFunction_phi_preComputed(SNES snes, Vec x, Vec F, void *ctx)
{
    PetscErrorCode ierr;
    NeuronGrowth *user = (NeuronGrowth *)ctx;

    //============================================
    // 1) Scatter 'x' into a sequential vector
    //    so we can read off the solution values.
    //============================================
    Vec P_seq;
    VecScatter scatter_ctx1;
    PetscScalar *Parray;
    ierr = VecScatterCreateToAll(x, &scatter_ctx1, &P_seq); CHKERRQ(ierr);
    ierr = VecScatterBegin(scatter_ctx1, x, P_seq, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterEnd(scatter_ctx1, x, P_seq, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGetArray(P_seq, &Parray); CHKERRQ(ierr);

    //============================================
    // 2) Initialize the residual vector F to zero
    //============================================
    ierr = VecSet(F, 0.0); CHKERRQ(ierr);

    //============================================
    // 3) Basic sizing & indexing
    //============================================
    const size_t nel = user->bzmesh_process.size();  // # of elements
    const int gptSize = (int)user->Gpt.size();        // Gauss points per dim
    int ind = 0;                            // Index into precomputed arrays

    //============================================
    // 4) Loop over elements to assemble residual
    //    and Jacobian contributions
    //============================================
    for (int e = 0; e < (int)nel; e++) {

        // (a) Get # of nodes in this element
        const int nen = (int)user->bzmesh_process[e].IEN.size();

        // (b) Reset EVectorSolve and EMatrixSolve for local assembly
        for (int m = 0; m < nen; m++) {
            user->EVectorSolve[m] = 0.0f;
            // Only reset up to 'nen' columns in each row
            fill(user->EMatrixSolve[m].begin(), user->EMatrixSolve[m].begin() + nen, 0.0f);
        }
        // Also reset user->eleVal[0] if it’s reused as a scratch
        fill(user->eleVal[0].begin(), user->eleVal[0].end(), 0.0f);

        // (c) Copy the current phi-guess from Parray into pre_eleVal[e][0]
        const auto &IEN = user->bzmesh_process[e].IEN;
        for (int ii = 0; ii < nen; ii++) {
            int A = IEN[ii];
            user->pre_eleVal[e][0][ii] = static_cast<float>(Parray[A]);
        }

        // (d) Triple nested loop over Gauss points
        for (int i = 0; i < gptSize; i++) {
            for (int j = 0; j < gptSize; j++) {
                for (int k = 0; k < gptSize; k++) {

                    // i) Evaluate phi at this Gauss point (elePG) and its gradient
                    float elePG = 0.0f;   // phi guess at Gauss point
                    float dPdx  = 0.0f;
                    float dPdy  = 0.0f;
                    float dPdz  = 0.0f;

                    user->ElementValue(user->pre_Nx[ind], user->pre_eleVal[e][0], elePG); // new phi guess 
                    user->ElementDeriv(nen, user->pre_dNdx[ind], user->pre_eleVal[e][0], dPdx, dPdy, dPdz);

                    // ii) Compute local PDE contributions for each node m
                    for (int m = 0; m < nen; m++) {

                        //-------------------------------------------
                        // * PDE Residual:
                        //-------------------------------------------
                        user->EVectorSolve[m] += (
                            elePG * user->pre_Nx[ind][m]
                            - user->dt * user->pre_eleMp[ind] * (
                                // 1) Anisotropy term: - (A^2)(grad(phi)·grad(Nx[m]))
                                - user->pre_eleAniso[ind] * user->pre_eleAniso[ind]
                                  * (  dPdx * user->pre_dNdx[ind][m][0]
                                     + dPdy * user->pre_dNdx[ind][m][1]
                                     + dPdz * user->pre_dNdx[ind][m][2]
                                    )
                                // 2) Additional terms: derivs wrt A in dA_dPdx, etc.
                                + ( - user->pre_dNdx[ind][m][0] * user->pre_eleAniso[ind]
                                    * user->pre_dA_dPdx[ind]
                                    * ( dPdx*dPdx + dPdy*dPdy + dPdz*dPdz )
                                  )
                                + ( - user->pre_dNdx[ind][m][1] * user->pre_eleAniso[ind]
                                    * user->pre_dA_dPdy[ind]
                                    * ( dPdx*dPdx + dPdy*dPdy + dPdz*dPdz )
                                  )
                                + ( - user->pre_dNdx[ind][m][2] * user->pre_eleAniso[ind]
                                    * user->pre_dA_dPdz[ind]
                                    * ( dPdx*dPdx + dPdy*dPdy + dPdz*dPdz )
                                  )
                                // 3) Polynomial in phi
                                + ( - elePG*elePG*elePG
                                    + (1 - user->pre_C1[ind])*elePG*elePG
                                    + user->pre_C1[ind]*elePG
                                  ) * user->pre_Nx[ind][m]
                            )
                            // Subtract old-phase at this Gauss point (user->pre_eleP[ind])
                            - user->pre_eleP[ind] * user->pre_Nx[ind][m]
                        ) * user->pre_detJ[ind];

                        //-------------------------------------------
                        // * PDE Jacobian:
                        //-------------------------------------------
                        for (int n = 0; n < nen; n++) {
                            user->EMatrixSolve[m][n] += (
                                user->pre_Nx[ind][m] * user->pre_Nx[ind][n]
                                - user->dt * user->pre_eleMp[ind] * (
                                    // a) Anisotropy term
                                    - user->pre_eleAniso[ind] * user->pre_eleAniso[ind]
                                      * (  user->pre_dNdx[ind][m][0] * user->pre_dNdx[ind][n][0]
                                         + user->pre_dNdx[ind][m][1] * user->pre_dNdx[ind][n][1]
                                         + user->pre_dNdx[ind][m][2] * user->pre_dNdx[ind][n][2]
                                        )
                                    // b) Derivs wrt orientation
                                    + ( - user->pre_dNdx[ind][m][0] * user->pre_eleAniso[ind]
                                        * user->pre_dA_dPdx[ind]
                                        * (2 * dPdx + 2 * dPdy + 2 * dPdz)
                                      )
                                    + ( - user->pre_dNdx[ind][m][1] * user->pre_eleAniso[ind]
                                        * user->pre_dA_dPdy[ind]
                                        * (2 * dPdx + 2 * dPdy + 2 * dPdz)
                                      )
                                    + ( - user->pre_dNdx[ind][m][2] * user->pre_eleAniso[ind]
                                        * user->pre_dA_dPdz[ind]
                                        * (2 * dPdx + 2 * dPdy + 2 * dPdz)
                                      )
                                    // c) Polynomial terms in phi
                                    + ( -3 * elePG*elePG
                                        + 2 * (1 - user->pre_C1[ind])*elePG
                                        + user->pre_C1[ind]*user->pre_Nx[ind][m]
                                      ) * user->pre_Nx[ind][n]
                                )
                            ) * user->pre_detJ[ind];
                        } // end n-loop
                    } // end m-loop

                    // (e) Move on to next Gauss point
                    ind++;
                }
            }
        }

        // (f) Apply boundary conditions (e.g. Dirichlet) if cpts[A].label == 1
        for (int ii = 0; ii < nen; ii++) {
            int A = user->bzmesh_process[e].IEN[ii];
            if (user->cpts[A].label == 1) {
                user->ApplyBoundaryCondition(0, ii, 0, user->EMatrixSolve, user->EVectorSolve);
            }
        }

        // (g) Assemble local residual & matrix into global F, user->J
        user->ResidualAssembly(user->EVectorSolve, user->bzmesh_process[e].IEN, F);
        user->MatrixAssembly(user->EMatrixSolve, user->bzmesh_process[e].IEN, user->J);
    }

    //============================================
    // 5) Finalize assembly for PETSc
    //============================================
    ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
    ierr = MatAssemblyBegin(user->J, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(user->J, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

    //============================================
    // 6) Cleanup
    //============================================
    ierr = VecRestoreArray(P_seq, &Parray); CHKERRQ(ierr);
    ierr = VecScatterDestroy(&scatter_ctx1); CHKERRQ(ierr);
    ierr = VecDestroy(&P_seq); CHKERRQ(ierr);

    return 0;
}

// FormJacobian for the phase-field equation
PetscErrorCode FormJacobian_phi(SNES snes, Vec x, Mat J, Mat P, void *ctx) {
    NeuronGrowth *user = static_cast<NeuronGrowth *>(ctx);
    PetscErrorCode ierr;

    // Copy precomputed Jacobian to PETSc matrix
    ierr = MatCopy(user->J, J, SAME_NONZERO_PATTERN); CHKERRQ(ierr);

    return PETSC_SUCCESS;
}

// Monitor function for SNES solver iterations
PetscErrorCode MySNESMonitor(SNES snes, PetscInt its, PetscReal fnorm, PetscViewerAndFormat *vf) {
    PetscFunctionBeginUser; // PETSc macro to mark the start of a user-defined function

    // Call PETSc's default short summary monitor for SNES iterations
    SNESMonitorDefaultShort(snes, its, fnorm, vf);

    // Retrieve the associated KSP solver from the SNES solver
    KSP ksp;
    PetscErrorCode ierr = SNESGetKSP(snes, &ksp); CHKERRQ(ierr);

    // Uncomment the desired KSP monitoring option
    // PetscOptionsSetValue(NULL, "-ksp_monitor", "");               // Monitor each KSP iteration
    PetscOptionsSetValue(NULL, "-ksp_monitor_singular_value", ""); // Monitor singular values of the KSP operator
    // PetscOptionsSetValue(NULL, "-ksp_monitor_solution", "");      // Monitor the KSP solution

    // Get the total number of KSP iterations for the current SNES solve
    PetscInt numIterations;
    ierr = KSPGetTotalIterations(ksp, &numIterations); CHKERRQ(ierr);

    // Print the total number of KSP iterations for debugging or analysis
    PetscPrintf(PETSC_COMM_WORLD, "     - KSP Iterations: %d\n", numIterations);

    PetscFunctionReturn(PETSC_SUCCESS); // Indicate successful execution
}

// Cleans up solvers and associated resources in the NeuronGrowth object
PetscErrorCode CleanUpSolvers(NeuronGrowth &NG) {
	if (NG.phi_solver == "snes") {
		// Safely destroy SNES solver for phi
		CHKERRQ(SNESDestroy(&NG.snes_phi));
		CHKERRQ(MatDestroy(&NG.J));
	} else {
		// Safely destroy KSP solver and resources for synaptogenesis (syn)
		CHKERRQ(KSPDestroy(&NG.ksp_phi));
		CHKERRQ(MatDestroy(&NG.GK_phi));
		CHKERRQ(VecDestroy(&NG.GR_phi));
	}
    CHKERRQ(VecDestroy(&NG.temp_phi));

    // Safely destroy KSP solver and resources for synaptogenesis (syn)
    CHKERRQ(KSPDestroy(&NG.ksp_syn));
    CHKERRQ(MatDestroy(&NG.GK_syn));
    CHKERRQ(VecDestroy(&NG.GR_syn));
    CHKERRQ(VecDestroy(&NG.temp_syn));

    // Safely destroy KSP solver and resources for tubules (tub)
    CHKERRQ(KSPDestroy(&NG.ksp_tub));
    CHKERRQ(MatDestroy(&NG.GK_tub));
    CHKERRQ(VecDestroy(&NG.GR_tub));
    CHKERRQ(VecDestroy(&NG.temp_tub));

	PetscFunctionReturn(PETSC_SUCCESS); // Indicate successful execution
}

int RunNG(
    const int n_bzmesh, vector<vector<int>> ele_process_in,
    vector<Vertex3D> &cpts_initial, vector<Vertex3D> &cpts, vector<Vertex3D>& prev_cpts, vector<Vertex3D>& cpts_fine,
    string path_in, string path_out,
    int &iter, int end_iter,
    vector<vector<float>> &NGvars,
    int &NX, int &NY, int &NZ,
    vector<array<float, 3>> &seed,
	int &originX, int &originY, int &originZ,
    bool &localRefine,
	const string& phi_solver,
	double& t_global,
	bool& restart, int& tmp_restart_check)
{
	/*========================================================*/
	// Initializations
	NeuronGrowth NG(phi_solver, iter, seed.size(), end_iter);
	NG.path_out = path_out;
	NG.SetVariables("simulation_parameters.txt"); // optional variable loading, for quick/batch simulation testing

	// Initialize vertex clouds for the current, fine, and previous configurations
	Vertex3DCloud cloud_initial(cpts_initial); 	// Cloud for initial points
	Vertex3DCloud cloud(cpts); 					// Cloud for current points
	Vertex3DCloud cloud_prev(prev_cpts);  		// Cloud for previous points
	Vertex3DCloud cloud_fine(cpts_fine);  		// Cloud for fine points
	// Initialize KD-Trees for the current, fine, and previous vertex clouds
	KDTree kdTree_initial(3 /* dim */, cloud_initial, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
	KDTree kdTree(3 /* dim */, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
	KDTree kdTree_prev(3 /* dim */, cloud_prev, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
	KDTree kdTree_fine(3 /* dim */, cloud_fine, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
	// Build indexes for KD-Trees to optimize search operations
	kdTree_initial.buildIndex();
	kdTree.buildIndex();
	kdTree_prev.buildIndex();
	kdTree_fine.buildIndex();

	// Call InitializeProblemNG with proper arguments
	NG.InitializeProblemNG(n_bzmesh, cpts, cloud, kdTree, prev_cpts, cloud_prev, kdTree_prev, NGvars, seed);
	// NG.ToPETScVec(NG.phi, NG.temp_phi); // initial guess for SNES (optional)
	PetscPrintf(PETSC_COMM_WORLD, "Set initial guess!-----------------------------------------------------------\n");	
	NG.AssignProcessor(ele_process_in);
	// Check MPI element assignments, and print out in orders
	for (int i = 0; i < NG.nProcess; i++) {
		if (i == NG.comRank) {
			cout << "[Rank " << NG.comRank << "/" << NG.nProcess << "] "
				<< "Element process size: " << NG.ele_process.size() << endl;
		}
		// Ensure barrier synchronization for clean output order
		CHKERRQ(MPI_Barrier(PETSC_COMM_WORLD));
	}
	// Read bezier mesh and prepare SNES initial guess
	NG.ReadBezierElementProcess(path_in);
	PetscPrintf(PETSC_COMM_WORLD, "Read bzmesh!-----------------------------------------------------------------\n");	

	/*========================================================*/
	// Write initial variables
	string varName;	
	if (NG.n == 0 && restart == false) {
		NG.VisualizeVTK_PhysicalDomain_All(0, NG.path_out);
		PetscPrintf(PETSC_COMM_WORLD, "Saving all variables!--------------------------------------------------------\n");		
	}

	/*========================================================*/
	// calculate some variable in advance to save computational cost
	NG.PrepareBasis();
	PetscPrintf(PETSC_COMM_WORLD, "Prepared basis!--------------------------------------------------------------\n");	
	CHKERRQ(MPI_Allreduce(&NG.sum_grad_phi0_local, &NG.sum_grad_phi0_global, 1, MPI_FLOAT, MPI_SUM, PETSC_COMM_WORLD));
	PetscPrintf(PETSC_COMM_WORLD, "Calculated sum grad phi0!----------------------------------------------------\n");
	NG.PrepareTermSource();
	PetscPrintf(PETSC_COMM_WORLD, "Prepared variables!----------------------------------------------------------\n");

	PetscPrintf(PETSC_COMM_WORLD, "*****************************************************************************************\n");
	PetscPrintf(PETSC_COMM_WORLD, "Running simmulations ... \n");	
	tic();
	while (iter <= NG.end_iter) {
		NG.n = iter;

		/*========================================================*/
		// if we want to restart the simulation
		if (restart == true) {
			restart = false;
			tmp_restart_check++;

			string restartVTK = FindLatestVTK(path_out);
			if (restartVTK == "") {
				cerr << "Failed to read the latest VTK file.\n";
			} else {
				// cout << "Read " << cpts.size() << " points from the largest-step file.\n";
				cout << "Reading " << cpts.size() << " points from:" << restartVTK << endl;;
			}
			NG.ReadVTK(restartVTK);
			
			cpts = NG.cpts;
			
			NGvars = {NG.phi, NG.syn, NG.tub, NG.theta, NG.phi_0, NG.tub_0};
			// Clean up solvers and synchronize processes
			// CHKERRQ(CleanUpSolvers(NG)); // some potential memory issues here

			CHKERRQ(MPI_Barrier(PETSC_COMM_WORLD));
			return 2;
		}
		
		/*--------------------------------------------------------*/
		// Domain expansion and variable passing
		if (NG.n % NG.expandCK_invl == 0 && NG.n >= 10) {
			// localRefine = true;

			NG.HandleExpansion(NG.phi, NX, NY, NZ, originX, originY, originZ);
			// Store NG variables
			NGvars = {NG.phi, NG.syn, NG.tub, NG.theta, NG.phi_0, NG.tub_0};

			// Clean up solvers and synchronize processes
			if (restart == false && tmp_restart_check == 1) {
			} else {
				CHKERRQ(CleanUpSolvers(NG));
			}
			// // Clean up solvers and synchronize processes
			// CHKERRQ(CleanUpSolvers(NG));
			NG.VisualizeVTK_ControlMesh(cpts, iter, path_out);

			if (NG.comRank == 0) {
				// Compute refinement values and save to file
				vector<float> ele_refine = NG.ComputeRefine(NG.phi, NX, NY, NZ, originX, originY, originZ, kdTree, cloud);
				writeVectorToFile(ele_refine, path_in + "phi.txt", false);
				vector<float> tmp = {NX, NY, NZ, originX, originY, originZ};
				writeVectorToFile(tmp, path_in + "domain_size.txt", false);
			}
			CHKERRQ(MPI_Barrier(PETSC_COMM_WORLD));

			iter++;
			return 2;
		}

		/*--------------------------------------------------------*/
		// Update local refinement information (2 scenarios, initial local refinement and later interval based local refinements)
		// if ((NG.n == 10 && !localRefine) || (NG.n % NG.refine_invl == 0 && NG.n != 0)) {
		if (NG.n == 10 && !localRefine) {
		// if (!localRefine) {
			PetscPrintf(PETSC_COMM_WORLD, "-----------------------------------------------------------------------------------------\n");
			PetscPrintf(PETSC_COMM_WORLD, "Checking and compuing local refinement information\n");
			// Store NG variables
			NGvars = {NG.phi, NG.syn, NG.tub, NG.theta, NG.phi_0, NG.tub_0};

			// Clean up solvers and synchronize across processes
			CleanUpSolvers(NG);
			CHKERRQ(MPI_Barrier(PETSC_COMM_WORLD));

			// Perform local refinement if applicable
			if (NG.comRank == 0) {
				vector<float> ele_refine = NG.ComputeRefine(NG.phi, NX, NY, NZ, originX, originY, originZ, kdTree, cloud);
				writeVectorToFile(ele_refine, path_in + "phi.txt", false);
			}
			CHKERRQ(MPI_Barrier(PETSC_COMM_WORLD));

			// Special case for simulation reset
			if (NG.n == 10 && !localRefine) {
				iter = 0;           // Reset iteration counter
				localRefine = true; // Enable local refinement
			} else {
				iter++;
			}
			return 1;
		}

		/*--------------------------------------------------------*/
		// Neuron identification and tip detection
		if ((NG.n % NG.tip_detect_invl == 0) || (NG.n == 0) || (NG.tips.size() != NG.phi.size()) || (NG.n == NG.end_iter)) {			
			// Detect tips and save intermediate results
			PetscPrintf(PETSC_COMM_WORLD, "-----------------------------------------------------------------------------------------\n");
			PetscPrintf(PETSC_COMM_WORLD, "Detecting tips\n");

			float tip_intensity_sz = 8.0f; // box size for calculating tip intensity
			// NG.DetectTips(cpts_fine, cloud_fine, kdTree_fine, tip_intensity_sz, cpts, cloud, kdTree, seed, NX, NY, NZ, originX, originY, originZ);
			NG.DetectTips_multi(cpts_fine, cloud_fine, kdTree_fine, tip_intensity_sz, cpts, cloud, kdTree, seed, NX, NY, NZ, originX, originY, originZ);
			PetscPrintf(PETSC_COMM_WORLD, "-----------------------------------------------------------------------------------------\n");
		}

		// Write physical domain results to file
		if (NG.n != 0 && NG.n % NG.var_save_invl == 0) {
			PetscPrintf(PETSC_COMM_WORLD, "-----------------------------------------------------------------------------------------\n");
			NG.VisualizeVTK_PhysicalDomain_All(NG.n, path_out);
			NG.VisualizeVTK_ControlMesh(cpts, iter, path_out);
			PetscPrintf(PETSC_COMM_WORLD, 
						"Step: %d/%d | Wrote Physical Domain! | Average time %fs | Total time: %f |\n", 
						NG.n, NG.end_iter, t_write / NG.var_save_invl, t_global);
			PetscPrintf(PETSC_COMM_WORLD, "-----------------------------------------------------------------------------------------\n");
		}

		UpdateSimulationTimers(t_collect, t_write, t_global);

		/*==============================================================================*/
		/* Implicit Nonlinear Solver for Phase Field Equation */
		PetscInt its_phi;               // KSP iteration count
		KSPConvergedReason ksp_reason_phi;  // KSP convergence reason
		SNESConvergedReason snes_reason_phi; // SNES convergence reason
		if (NG.phi_solver == "ksp") { // Unstable, susceptable to different parameter adjustments, use SNES
			/*--------------------------------------------------------*/
			/* Implicit Nonlinear Newton-Raphson Solver with Adaptive Time-Stepping for Phase Field Equation */

			// Time-stepping parameters
			double dt_min = 1e-6;      // Minimum time step
			double dt_max = 1e-2;      // Maximum time step
			double factor_grow = 1.5;  // Time step growth factor
			double factor_shrink = 0.5; // Time step shrink factor
			double tol_increase = 1e-3; // Residual tolerance for increasing dt
			double tol_decrease = 1e-1; // Residual tolerance for decreasing dt

			NG.phi_prev = NG.phi; // Save previous phi state for comparison
			NG.PreparePhaseField(); // Prepare necessary data structures for solving

			// Initialize residuals and other variables
			double tol = 1e-4;   // Convergence tolerance
			double residual = INF;  // Initialize residual
			int NR_itr = 0;      // Newton-Raphson iteration counter
			int max_NR_iters = 50;  // Maximum NR iterations allowed per time step

			while (residual > tol) {
				// === Step 1: Reset Linear System ===
				NG.ierr = MatZeroEntries(NG.GK_phi); CHKERRQ(NG.ierr);
				NG.ierr = VecSet(NG.GR_phi, 0); CHKERRQ(NG.ierr);

				// Build the linear system for the current iteration
				NG.BuildLinearSystemProcessNG_phi(cpts);

				// Finalize assembly of matrix and residual vector
				NG.ierr = VecAssemblyEnd(NG.GR_phi); CHKERRQ(NG.ierr);
				NG.ierr = MatAssemblyEnd(NG.GK_phi, MAT_FINAL_ASSEMBLY); CHKERRQ(NG.ierr);

				// === Step 2: Initialize KSP Solver (if not done) ===
				if (NG.judge_phi == 0) {
					CHKERRQ(SetupKSP(NG.ksp_phi, NG.GK_phi, KSPGMRES, PCBJACOBI, 1.e-6));
					if (NG.n == 0) {
						CHKERRQ(KSPView(NG.ksp_phi, PETSC_VIEWER_STDOUT_WORLD));
					}
					NG.judge_phi = 1; // Mark solver as initialized
				}

				// === Step 3: Solve Linear System ===
				NG.ierr = KSPSolve(NG.ksp_phi, NG.GR_phi, NG.temp_phi); CHKERRQ(NG.ierr);

				// Retrieve solver statistics
				NG.ierr = KSPGetIterationNumber(NG.ksp_phi, &its_phi); CHKERRQ(NG.ierr);
				NG.ierr = KSPGetConvergedReason(NG.ksp_phi, &ksp_reason_phi); CHKERRQ(NG.ierr);

				// Check if KSP solver converged
				if (ksp_reason_phi < 0) {
					PetscPrintf(PETSC_COMM_WORLD, 
								"[ERROR] KSP did not converge. Reason: %d | Iterations: %d\n", 
								ksp_reason_phi, its_phi);

					// Reduce time step and restart if solver fails
					NG.dt = PetscMax(NG.dt * factor_shrink, dt_min);
					// PetscPrintf(PETSC_COMM_WORLD, "[INFO] Reducing time step to %.5e and restarting.\n", NG.dt);
					NR_itr = 0;  // Reset NR iteration counter
					continue;    // Restart NR solver with reduced time step
				}

				// === Step 4: Update Solution and Compute Residual ===
				Vec temp_phi_seq;
				VecScatter ctx_phi;
				PetscScalar *_p;
				// Scatter solution to all processors
				NG.ierr = VecScatterCreateToAll(NG.temp_phi, &ctx_phi, &temp_phi_seq); CHKERRQ(NG.ierr);
				NG.ierr = VecScatterBegin(ctx_phi, NG.temp_phi, temp_phi_seq, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(NG.ierr);
				NG.ierr = VecScatterEnd(ctx_phi, NG.temp_phi, temp_phi_seq, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(NG.ierr);
				// Retrieve array for computation
				NG.ierr = VecGetArray(temp_phi_seq, &_p); CHKERRQ(NG.ierr);

				// Compute residual and update phi
				double prev_residual = residual; // Save previous residual
				residual = 0;
				for (size_t i = 0; i < NG.phi.size(); i++) {
					double delta_phi = PetscRealPart(_p[i]);
					NG.phi[i] -= delta_phi; // Update phi
					residual = PetscMax(residual, fabs(delta_phi)); // Update residual

					// Check for divergence
					if (fabs(NG.phi[i]) > 3) {
						PetscPrintf(PETSC_COMM_WORLD, 
									"[ERROR] Diverging phi detected! Stopping simulation.\n");
						return 3;
					}
				}

				// Restore and clean up scattered vector
				NG.ierr = VecRestoreArray(temp_phi_seq, &_p); CHKERRQ(NG.ierr);
				NG.ierr = VecScatterDestroy(&ctx_phi); CHKERRQ(NG.ierr);
				NG.ierr = VecDestroy(&temp_phi_seq); CHKERRQ(NG.ierr);
				// Increment NR iteration counter
				NR_itr++;

				// === Step 5: Adjust Time Step ===
				if (residual < tol_increase * prev_residual) {
					// Increase time step if residual decreases significantly
					NG.dt = PetscMin(NG.dt * factor_grow, dt_max);
					// PetscPrintf(PETSC_COMM_WORLD, "[INFO] Increasing time step to %.5e.\n", NG.dt);
				} else if (residual > tol_decrease * prev_residual) {
					// Reduce time step if residual decreases insufficiently
					NG.dt = PetscMax(NG.dt * factor_shrink, dt_min);
					// PetscPrintf(PETSC_COMM_WORLD, "[INFO] Reducing time step to %.5e.\n", NG.dt);
				}

				// Check for maximum NR iterations per step
				if (NR_itr > max_NR_iters) {
					PetscPrintf(PETSC_COMM_WORLD, 
								"[ERROR] Maximum Newton-Raphson iterations exceeded. Stopping.\n");
					return 3;
				}
			}
		} else if (NG.phi_solver == "snes") {
		/*--------------------------------------------------------*/
			/*Implcit Non-liear SNES solver for Phase field equation*/
			NG.phi_prev = NG.phi;		
			// NG.PreparePhaseField_SNES();
			NG.PreparePhaseField_SNES_preComputed();
			// Phase Field Equation Solver (SNES)
			if (NG.judge_phi == 0) {
				// SetupSNES(NG.snes_phi, SNESNEWTONLS, &NG, FormFunction_phi, FormJacobian_phi);
				SetupSNES(NG.snes_phi, SNESNEWTONLS, &NG, FormFunction_phi_preComputed, FormJacobian_phi);

				if (NG.n == 0) {
					CHKERRQ(SNESView(NG.snes_phi, PETSC_VIEWER_STDOUT_WORLD));
				}
				NG.judge_phi = 1;
			}

			CHKERRQ(SNESSolve(NG.snes_phi, NULL, NG.temp_phi));
			CHKERRQ(SNESGetIterationNumber(NG.snes_phi, &its_phi));
			CHKERRQ(SNESGetConvergedReason(NG.snes_phi, &snes_reason_phi));
			if (snes_reason_phi < 0) {
				PetscPrintf(PETSC_COMM_WORLD, "SNES phi not converging ........ | SNESreason: %d | SNESiter: %d \n", snes_reason_phi, its_phi);	
				return 3; // divering, ending simulation
			}
			/*Collecting scattered Phi variable from all processors*/
			// Scatter phi variable
			CHKERRQ(ScatterVector(NG.temp_phi, NG.phi, NG.phi.size()));
		}
		
		// Timing updates for Phi
		UpdateSimulationTimers(t_phi, t_write, t_global);
		
		/*==============================================================================*/
		/* Synaptogenesis and Tubulin Equations Solver */

		// Synchronize and update global values
		CHKERRQ(MPI_Allreduce(&NG.sum_grad_phi0_local, &NG.sum_grad_phi0_global, 1, MPI_FLOAT, MPI_SUM, PETSC_COMM_WORLD));
		// Prepare source terms
		NG.PrepareTermSource();
		// Build Linear Systems
		if (NG.judge_syn == 0) {
			CHKERRQ(MatZeroEntries(NG.GK_syn)); // Reset Synaptogenesis matrix if needed
		}
		CHKERRQ(VecSet(NG.GR_syn, 0));          // Reset Synaptogenesis RHS
		CHKERRQ(MatZeroEntries(NG.GK_tub));    // Reset Tubulin matrix
		CHKERRQ(VecSet(NG.GR_tub, 0));         // Reset Tubulin RHS
		// Build both Synaptogenesis and Tubulin systems in one loop
		NG.BuildLinearSystemProcessNG_syn_tub(cpts); // Setup both equations together to same time and computational costs
		// Assemble RHS vectors and matrices
		CHKERRQ(VecAssemblyEnd(NG.GR_syn)); // Synaptogenesis RHS
		if (NG.judge_syn == 0) {
			CHKERRQ(MatAssemblyEnd(NG.GK_syn, MAT_FINAL_ASSEMBLY)); // Synaptogenesis matrix
		}
		CHKERRQ(VecAssemblyEnd(NG.GR_tub)); // Tubulin RHS
		CHKERRQ(MatAssemblyEnd(NG.GK_tub, MAT_FINAL_ASSEMBLY));    // Tubulin matrix

		/*--------------------------------------------------------*/
		/* Synaptogenesis Solver (KSP) */
		// Initialize KSP solver for Synaptogenesis if needed
		if (NG.judge_syn == 0) {
			CHKERRQ(SetupKSP(NG.ksp_syn, NG.GK_syn, KSPGMRES, PCBJACOBI, 1.e-8)); // GMRES with Jacobi preconditioner
			if (NG.n == 0) {
				CHKERRQ(KSPView(NG.ksp_syn, PETSC_VIEWER_STDOUT_WORLD)); // View KSP configuration for debugging
			}
			NG.judge_syn = 1; // Mark as initialized
		}
		// Solve Synaptogenesis equation
		CHKERRQ(KSPSolve(NG.ksp_syn, NG.GR_syn, NG.temp_syn)); // Solve system
		PetscInt its_syn;
		KSPConvergedReason reason_syn;
		CHKERRQ(KSPGetIterationNumber(NG.ksp_syn, &its_syn));    // Get iteration count
		CHKERRQ(KSPGetConvergedReason(NG.ksp_syn, &reason_syn)); // Get convergence reason
		// Check for divergence
		if (reason_syn < 0) {
			PetscPrintf(PETSC_COMM_WORLD, "KSP Synaptogenesis not converging  | KSPreason: %d | KSPiter: %d \n", reason_syn, its_syn);
			return 3; // Divergence detected
		}

		/*--------------------------------------------------------*/
		/* Tubulin Solver (KSP) */
		// Initialize KSP solver for Tubulin if needed
		if (NG.judge_tub == 0) {
			CHKERRQ(SetupKSP(NG.ksp_tub, NG.GK_tub, KSPGMRES, PCBJACOBI, 1.e-8)); // GMRES with Jacobi preconditioner
			if (NG.n == 0) {
				CHKERRQ(KSPView(NG.ksp_tub, PETSC_VIEWER_STDOUT_WORLD)); // View KSP configuration for debugging
			}
			NG.judge_tub = 1; // Mark as initialized
		}
		// Solve Tubulin equation
		CHKERRQ(KSPSolve(NG.ksp_tub, NG.GR_tub, NG.temp_tub)); // Solve system
		PetscInt its_tub;
		KSPConvergedReason reason_tub;
		CHKERRQ(KSPGetIterationNumber(NG.ksp_tub, &its_tub));    // Get iteration count
		CHKERRQ(KSPGetConvergedReason(NG.ksp_tub, &reason_tub)); // Get convergence reason
		// Check for divergence
		if (reason_tub < 0) {
			PetscPrintf(PETSC_COMM_WORLD, "KSP Tubulin not converging ..... | KSPreason: %d | KSPiter: %d \n", reason_tub, its_tub);
			return 3; // Divergence detected
		}

		// Timing updates for Tubulin
		UpdateSimulationTimers(t_syn_tub, t_write, t_global);

		/*Collecting scattered Synaptogenesis and Tubulin variables from all processors*/
		// Scatter synaptogenesis and tubulin variable
		CHKERRQ(ScatterVector(NG.temp_syn, NG.syn, NG.syn.size()));
		CHKERRQ(ScatterVector(NG.temp_tub, NG.tub, NG.tub.size(), true, &NG));

		/*==============================================================================*/
		/*Iteration summary printout*/	
		auto reason_phi = (NG.phi_solver == "snes") ? snes_reason_phi : ksp_reason_phi;
		NG.PrintStatus(NG.n, NG.end_iter, reason_phi, its_phi, t_phi, 
            reason_syn, reason_tub, its_syn, its_tub, t_syn_tub, 
            n_bzmesh);

		// Increment iteration counter if no expansion
		iter++;
	}

		/*==============================================================================*/
	CHKERRQ(CleanUpSolvers(NG)); // Destroy solvers
	CHKERRQ(MPI_Barrier(PETSC_COMM_WORLD));

	return 0; // ending simulation
}