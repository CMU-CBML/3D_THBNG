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
    PetscErrorCode ierr = initializeMPI(rank, nProcs, argc, argv, help); 
    if (ierr) return ierr;

    // Validate input arguments
    if (argc < 4) {
        cerr << "Usage: <numNeuron> <end_iter> <path_in>" << endl;
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
    vector<Element3D> tmesh_initial, tmesh;
    vector<vector<float>> NGvars(6); // Stores neuron growth variables

    bool localRefine = false; // Flag for local refinement
    int iter = 0, state = 1;  // Simulation state: 0-end, 1-running, 2-expanding, 3-diverging

    PetscPrintf(PETSC_COMM_WORLD, "Starting Simulation\n");

    // Main simulation loop
    while (iter <= end_iter) {
        // Reset previous simulation state
        prev_cpts = cpts;
        cpts_initial.clear();
        tmesh_initial.clear();
        cpts.clear();
        tmesh.clear();
        ele_process.clear();
        ele_process.resize(nProcs);

        if (rank == 0) {
            // Setup simulation files, generate mesh, and partition if needed
            setupSimulationFiles(nProcs, path_in, localRefine, vertices, elements, NX, NY, NZ, originX, originY, originZ);
        }

        ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(ierr); // Synchronize processes

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
            seed, originX, originY, originZ,
            localRefine,
            phi_solver);

        // Exit if simulation diverges
        if (state == 3) {
            PetscPrintf(PETSC_COMM_WORLD, "Simulation diverging, ending program.\n");
            break;
        }
    }

    // Finalize simulation
    PetscPrintf(PETSC_COMM_WORLD, "Simulation Complete!\n");
    ierr = PetscFinalize(); CHKERRQ(ierr);
    return 0;
}

// #include <iostream>
// #include <vector>
// #include <array>
// #include "NeuronGrowth.h"
// #include <sstream>
// #include <iomanip>
// // #include "time.h"
// #include "utils.h"

// using namespace std;

// static char help[] = "Solve 3DNG\n";

// int main(int argc, char **argv)
// {
// 	// // int numNeuron = atoi(argv[1]);		// user input number of neurons
// 	// // int end_iter = atoi(argv[2]);		// user input number of iterations
// 	// // string path_in = argv[3];		// user input working directory
// 	// // string fn_mesh_initial(path_in + "controlmesh_initial.vtk");
// 	// string fn_mesh_initial("../io3D/controlmesh_initial.vtk");
// 	// vector<vector<float>> vertices;
// 	// vector<vector<int>> elements;
// 	// gen3Dmesh(0,0,0, 10,10,10, vertices, elements);
// 	// write_hex_toVTK(fn_mesh_initial.c_str(), vertices, elements);

// 	int rank, nProcs;
// 	PetscErrorCode ierr;
// 	/// start up petsc
// 	ierr = PetscInitialize(&argc, &argv, (char*)0, help); if (ierr) return ierr;
// 	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
// 	ierr = MPI_Comm_size(PETSC_COMM_WORLD, &nProcs); CHKERRQ(ierr);

// 	int numNeuron = atoi(argv[1]);		// user input number of neurons
// 	int end_iter = atoi(argv[2]);		// user input number of iterations
// 	string path_in = argv[3];		// user input working directory

// 	int NX, NY, NZ, originX(0), originY(0), originZ(0);
// 	vector<array<float, 3>> seed; 
// 	InitializeSoma(numNeuron, seed, NX, NY, NZ);

// 	/// Set simulation parameters and mesh
// 	string fn_mesh_initial(path_in + "controlmesh_initial.vtk");
// 	string fn_mesh(path_in + "controlmesh.vtk");
// 	// string fn_mesh(path_in + "controlPoints.vtk");
// 	string fn_bz(path_in + "bzmeshinfo.txt.epart." + to_string(nProcs));
// 	string path_out(path_in + "outputs/");

// 	int n_bzmesh;
// 	vector<vector<float>> vertices;
// 	vector<vector<int>> elements, ele_process;
// 	vector<Vertex3D> cpts_initial, cpts, prev_cpts;
// 	vector<Element3D> tmesh_initial, tmesh;
// 	ele_process.resize(nProcs);
// 	PetscPrintf(PETSC_COMM_WORLD, "Initial element process size: %ld\n", ele_process.size()); CHKERRQ(ierr);

// 	// // vector<int> rfid, rftype; // list of elements to refine
// 	vector<vector<float>> NGvars; // for passing neuron growth variables in-between domain expansion
// 	NGvars.clear(); NGvars.resize(6);

// 	bool localRefine = false;
// 	// UserSetting *NGuser = new UserSetting;

// 	int iter(0), state(1); // 0-end, 1-running, 2-expanding, 3-diverging simulation
// 	while (iter <= end_iter) {
// 		prev_cpts = cpts; // back up old control points for later NGvars interpolations (old cpts to new cpts)
// 		cpts_initial.clear(); tmesh_initial.clear(); 
// 		cpts.clear(); tmesh.clear(); 
// 		ele_process.clear();
// 		ele_process.resize(nProcs);

// 		if (rank == 0) {
// 			// to make sure correct files are generated and then read later on
// 			std::remove("../io3D/controlmesh.vtk");
// 			std::remove("../io3D/controlPoints.vtk");
// 			std::remove("../io3D/controlmesh_initial.vtk");
// 			std::remove("../io3D/bzpt.txt");
// 			std::remove("../io3D/cmat.txt");
// 			std::remove("../io3D/bzmesh.vtk");
// 			std::remove("../io3D/bzmeshinfo.txt");
// 			string epart("../io3D/bzmeshinfo.txt.epart." + std::to_string(nProcs));
// 			string npart("../io3D/bzmeshinfo.txt.npart." + std::to_string(nProcs));
// 			std::remove(epart.c_str());
// 			std::remove(npart.c_str());
// 			// std::remove("../io3D/phi.txt");

// 			gen3Dmesh(originX, originY, originZ, NX, NY, NZ, vertices, elements); // Generating 3D hex mesh
// 			write_hex_toVTK(fn_mesh_initial.c_str(), vertices, elements);
// 			// if (iter == 0) {
// 			// 	vector<float> ele_refine(elements.size(), 0);
// 			// 	writeVectorToFile(ele_refine, path_in + "phi.txt", false);
// 			// }

// 			if (localRefine != true) {
// 				write_hex_toVTK(fn_mesh.c_str(), vertices, elements);
// 				bzmesh3D(path_in); // generating 3D bezier mesh information (bzmeshinfo.txt)
// 			} 
// 			else {
// 				THS3D(path_in);
// 			}
// 			// write_hex_toVTK(fn_mesh.c_str(), vertices, elements);
// 			// THS3D(path_in);

// 			mpmetis(nProcs, path_in); // partitioning bzmeshinfo using mpmetis for parallelization
// 		}
// 		ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(ierr);

// 		// ReadMesh(fn_mesh_initial, cpts_initial, tmesh_initial);
// 		// ReadMesh(fn_mesh, cpts, tmesh);
	
// 		if (localRefine == true) {
// 			fn_mesh = path_in + "controlPoints.vtk";
// 		}

// 		ReadControlPoints(fn_mesh_initial, cpts_initial);
// 		ReadControlPoints(fn_mesh, cpts);

// 		AssignProcessor(fn_bz, n_bzmesh, ele_process);
// 		PetscPrintf(PETSC_COMM_WORLD, "Processor Assigned!----------------------------------------------------------\n");

// 		// state = RunNG(n_bzmesh, ele_process, cpts_initial, tmesh_initial, cpts, prev_cpts, tmesh, path_in, path_out,
// 		// 	iter, end_iter, NGvars, NX, NY, NZ, seed, originX, originY, originZ, localRefine);
// 		state = RunNG(n_bzmesh, ele_process, cpts_initial, cpts, prev_cpts, path_in, path_out,
// 			iter, end_iter, NGvars, NX, NY, NZ, seed, originX, originY, originZ, localRefine);
// 		// return 0;

// 		if (state == 3) {
// 			PetscPrintf(PETSC_COMM_WORLD, "Simulation divering, ending program.\n");
// 			return 0;
// 		}
// 	}

// 	PetscPrintf(PETSC_COMM_WORLD, "Done - Main reached end iteration!\n");
// 	ierr = PetscFinalize(); CHKERRQ(ierr);
// 	return 0;
// }
