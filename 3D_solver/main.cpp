#include <iostream>
#include <vector>
#include <array>
#include <sstream>
#include <iomanip>
#include "NeuronGrowth.h"
#include "utils.h"

using namespace std;

static int printUsageAndExit(){
    PetscPrintf(PETSC_COMM_WORLD,
        "\nUsage:\n"
        "  mpirun -np <N> ./3DNG \\\n"
        "    --numNeuron=<int> \\\n"
        "    --end_iter=<int> \\\n"
        "    --solver=<snes|ksp> \\\n"
        "    --path_in=<dir> \\\n"
        "    [--restart=<yes|no|true|false|1|0>]\n\n"
        "Example:\n"
        "/xxx/xxx/xxx/mpirun -np 14 ./3DNG --numNeuron=1 --end_iter=20000 --solver=snes "
        "--path_in=../io3D/ --restart=yes\n\n"
    );
    PetscFinalize();
    return EXIT_FAILURE;
}

PetscErrorCode initializeMPI(int &rank, int &nProcs, int argc, char** argv, const char help[] = "");

int main(int argc, char** argv)
{
    int rank = 0, nProcs = 1;
    const char help[] = "3D Neuron Growth Solver.\n";
    // If using PETSc initialization:
    CHKERRQ(initializeMPI(rank, nProcs, argc, argv, help));

    //--------------------------------------------------------------------------
    // 1) Defaults (set them, or leave them uninitialized if you want them required)
    //--------------------------------------------------------------------------
    bool haveNumNeuron = false;
    bool haveEndIter   = false;
    bool haveSolver    = false;
    bool havePathIn    = false;

    int numNeuron  = 0;
    int end_iter   = 0;
    string solver;
    string path_in;
    bool restart   = false; // optional, default = false

    //--------------------------------------------------------------------------
    // 2) Parse arguments of the form: --key=value
    //    We skip argv[0], which is the program name ./3DNG
    //--------------------------------------------------------------------------
    if (argc == 1) {
        // No arguments provided
        cerr << "Error: No arguments provided.\n";
        return printUsageAndExit();
    }

    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg.rfind("--", 0) != 0) {
            // Doesn't start with '--', error
            cerr << "Error: Argument '" << arg << "' must start with '--'.\n";
            return printUsageAndExit();
        }

        // Remove the leading '--'
        arg.erase(0, 2); // now arg is "key=value"

        // Find '='
        size_t eqPos = arg.find('=');
        if (eqPos == string::npos) {
            // No '=' found
            cerr << "Error: Argument '--" << arg
                      << "' is missing '=value'.\n";
            return printUsageAndExit();
        }

        // Split key and value
        string key   = arg.substr(0, eqPos);
        string value = arg.substr(eqPos + 1);

        // Convert key to lowercase if you want case-insensitivity
        transform(key.begin(), key.end(), key.begin(),
                       [](unsigned char c){ return tolower(c); });

        //--------------------------------------------------------------------------
        // 3) Match recognized keys and parse values
        //--------------------------------------------------------------------------
        if (key == "numneuron") {
            numNeuron = atoi(value.c_str());
            haveNumNeuron = true;
        }
        else if (key == "end_iter") {
            end_iter = atoi(value.c_str());
            haveEndIter = true;
        }
        else if (key == "solver") {
            solver = value;
            // Validate
            if (solver != "snes" && solver != "ksp") {
                cerr << "Error: solver must be 'snes' or 'ksp'.\n";
                return printUsageAndExit();
            }
            haveSolver = true;
        }
        else if (key == "path_in") {
            path_in = value;
            havePathIn = true;
        }
        else if (key == "restart") {
            // Accept "yes", "true", "1" => true; else false
            string valLower = value;
            transform(valLower.begin(), valLower.end(), valLower.begin(),
                           [](unsigned char c){ return tolower(c); });
            if (valLower == "yes" || valLower == "true" || valLower == "1") {
                restart = true;
            } else {
                restart = false;
            }
        }
        else {
            // Unknown key => warn or error out
            cerr << "Warning: Unrecognized argument key '--" << key << "'. Ignored.\n";
        }
    }

    //--------------------------------------------------------------------------
    // 4) Check if all required options were specified
    //--------------------------------------------------------------------------
    if (!haveNumNeuron || !haveEndIter || !haveSolver || !havePathIn) {
        cerr << "Error: Missing one or more required arguments.\n\n";
        return printUsageAndExit();
    }

    // Path out could be derived from path_in, e.g.:
    string path_out = path_in + "outputs/";

    //--------------------------------------------------------------------------
    // 5) Print final configuration (optional)
    //--------------------------------------------------------------------------
    PetscPrintf(PETSC_COMM_WORLD,
                "Configuration:\n"
                "  numNeuron  = %d\n"
                "  end_iter   = %d\n"
                "  solver     = %s\n"
                "  path_in    = %s\n"
                "  path_out   = %s\n"
                "  restart    = %s\n",
                numNeuron, end_iter, solver.c_str(),
                path_in.c_str(), path_out.c_str(),
                restart ? "true" : "false");
                
    // Simulation parameters
    bool localRefine = false; // Flag for local refinement
    int iter = 0, state = 1;  // Simulation state: 0-end, 1-running, 2-expanding, 3-diverging
    int NX, NY, NZ, originX(0), originY(0), originZ(0);
    vector<array<float, 3>> seed;
    if (restart == false) {
        InitializeSoma(numNeuron, seed, NX, NY, NZ);
    } else {
        localRefine = true;
        vector<float> domain_size = readVectorFromFile(path_in + "domain_size.txt", false);
        NX = (int)domain_size[0];
        NY = (int)domain_size[1];
        NZ = (int)domain_size[2];
        originX = (int)domain_size[3];
        originY = (int)domain_size[4];
        originZ = (int)domain_size[5];
        cout << "Read origin: " << originX << " " << originY << " " << originZ << endl;
        
        string latestVTK = FindLatestVTK(path_out);
        iter = getStepFromVTK(latestVTK.substr(path_out.size(), latestVTK.size() - path_out.size()));
    }

    // Data structures for mesh and simulation
    int n_bzmesh;
    vector<vector<float>> vertices;
    vector<vector<int>> elements, ele_process(nProcs);
    vector<Vertex3D> cpts_initial, cpts, prev_cpts, cpts_fine;
    vector<vector<float>> NGvars(6); // Stores neuron growth variables
    // vector<Element3D> mesh;

    PetscPrintf(PETSC_COMM_WORLD, "Starting Simulation\n");
    double t_global = 0;

    int tmp_restart_check = 0;
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
        string fn_mesh_fine = path_in + "controlmesh_fine.vtk";
        string fn_bz = path_in + "bzmeshinfo.txt.epart." + to_string(nProcs);

        // Read control points and assign processors
        ReadControlPoints(fn_mesh_initial, cpts_initial);
        ReadControlPoints(fn_mesh, cpts);
        // ReadMesh(fn_mesh, cpts, mesh);
        // cout << "Finishing reading -===================" << endl;
        ReadControlPoints(fn_mesh_fine, cpts_fine);
        AssignProcessor(fn_bz, n_bzmesh, ele_process);

        PetscPrintf(PETSC_COMM_WORLD, "Processor Assigned!\n");

        // Run neuron growth simulation for the current iteration
        state = RunNG(
            n_bzmesh, ele_process, 
            cpts_initial, cpts, prev_cpts, cpts_fine,
            path_in, path_out,
            iter, end_iter,
            NGvars,
            NX, NY, NZ,
            seed,
            originX, originY, originZ,
            localRefine,
            solver,
            t_global,
            restart, tmp_restart_check);

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
