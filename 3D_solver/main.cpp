#include <iostream>
#include <vector>
#include <array>
#include <sstream>
#include <iomanip>
#include "NeuronGrowth.h"
#include "utils.h"

using namespace std;

static char help[] = "Solve 3DNG\n";

int main(int argc, char** argv) {
    int rank, nProcs;
    CHKERRQ(initializeMPI(rank, nProcs, argc, argv, help));

    // Validate input arguments
    if (argc < 4) {
        cerr << "Usage: <numNeuron> <end_iter> <solver> <path_in>" << endl;
        return 1;
    }

    // User inputs
    int numNeuron = atoi(argv[1]);   // Number of neurons
    int end_iter = atoi(argv[2]);   // Number of iterations
    string phi_solver = string(argv[3]);
    if (phi_solver != "ksp" && phi_solver != "snes") {
        PetscPrintf(PETSC_COMM_WORLD, "Please specify phi solver type: snes or ksp.\n");
        return 0;
    }
    string path_in = argv[4];       // Working directory path
    string path_out = path_in + "outputs/";

    // Simulation parameters
    int NX, NY, NZ, originX(0), originY(0), originZ(0);
    vector<array<float, 3>> seed;
    InitializeSoma(numNeuron, seed, NX, NY, NZ);

    // Data structures for mesh and simulation
    int n_bzmesh;
    vector<vector<float>> vertices;
    vector<vector<int>> elements, ele_process(nProcs);
    vector<Vertex3D> cpts_initial, cpts, prev_cpts;
    vector<vector<float>> NGvars(6); // Stores neuron growth variables

    bool localRefine = false; // Flag for local refinement
    int iter = 0, state = 1;  // Simulation state: 0-end, 1-running, 2-expanding, 3-diverging

    PetscPrintf(PETSC_COMM_WORLD, "Starting Simulation\n");
    double t_global = 0;

    // Main simulation loop
    while (iter <= end_iter) {
        // Reset simulation state
        prev_cpts = move(cpts); // Efficiently transfer ownership instead of copying
        cpts.clear();
        cpts_initial.clear();

        // Clear and resize `ele_process` for parallel processing
        ele_process.assign(nProcs, {}); // Clear and resize in one step

        if (rank == 0) {
            // Setup simulation files, generate mesh, and partition if needed
            SetupSimulationFiles(nProcs, path_in, localRefine, vertices, elements, NX, NY, NZ, originX, originY, originZ);
        }

        CHKERRQ(MPI_Barrier(PETSC_COMM_WORLD)); // Synchronize processes

        // File paths for reading control points
        string fn_mesh_initial = path_in + "controlmesh_initial.vtk";
        string fn_mesh = localRefine ? path_in + "controlPoints.vtk" : path_in + "controlmesh.vtk";
        string fn_bz = path_in + "bzmeshinfo.txt.epart." + to_string(nProcs);

        // Read control points and assign processors
        ReadControlPoints(fn_mesh_initial, cpts_initial);
        ReadControlPoints(fn_mesh, cpts);
        AssignProcessor(fn_bz, n_bzmesh, ele_process);

        PetscPrintf(PETSC_COMM_WORLD, "Processor Assigned!\n");

        // Run neuron growth simulation for the current iteration
        state = RunNG(
            n_bzmesh, ele_process, 
            cpts_initial, cpts, prev_cpts, 
            path_in, path_out,
            iter, end_iter,
            NGvars,
            NX, NY, NZ,
            seed,
            originX, originY, originZ,
            localRefine,
            phi_solver,
            t_global);

        // Exit if simulation diverges
        if (state == 3) {
            PetscPrintf(PETSC_COMM_WORLD, "Simulation diverging, ending program.\n");
            break;
        }
    }

    // Finalize simulation
    PetscPrintf(PETSC_COMM_WORLD, "Simulation Complete!\n");
    CHKERRQ(PetscFinalize());
    return 0;
}