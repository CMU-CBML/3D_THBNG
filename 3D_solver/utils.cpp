#include "utils.h"
#include <iostream>
#include <algorithm>
#include "BasicDataStructure.h"
#include <cmath>

// Removes a list of files from the filesystem
void removeFiles(const vector<string>& files) {
    for (const auto& file : files) {
        remove(file.c_str());
    }
}

// Initializes MPI environment and retrieves rank and size
PetscErrorCode initializeMPI(int& rank, int& nProcs, int argc, char** argv, const char help[]) {
    PetscErrorCode ierr;
    ierr = PetscInitialize(&argc, &argv, nullptr, help); CHKERRQ(ierr);
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &nProcs); CHKERRQ(ierr);
    return ierr;
}

// Sets up simulation files, generates mesh, and partitions if needed
void SetupSimulationFiles(const int nProcs, const string& path_in, bool localRefine,
                          vector<vector<float>>& vertices, vector<vector<int>>& elements,
                          int NX, int NY, int NZ, int originX, int originY, int originZ) {

    // Remove old files to ensure fresh outputs
    const vector<string> filesToRemove = {
        path_in + "controlmesh.vtk",
        path_in + "controlPoints.vtk",
        path_in + "controlmesh_initial.vtk",
        path_in + "bzpt.txt",
        path_in + "cmat.txt",
        path_in + "bzmesh.vtk",
        path_in + "bzmeshinfo.txt",
        path_in + "bzmeshinfo.txt.epart." + to_string(nProcs),
        path_in + "bzmeshinfo.txt.npart." + to_string(nProcs),
        path_in + "phi_10.txt." + to_string(nProcs),
        path_in + "phi_refine.txt." + to_string(nProcs),
        // path_in + "phi.txt." + to_string(nProcs)
    };
    removeFiles(filesToRemove);

    // File paths for the meshes
    string fn_mesh_initial = path_in + "controlmesh_initial.vtk";
    string fn_mesh = path_in + "controlmesh.vtk";

    // Generate the initial 3D mesh and write it to file
    gen3Dmesh(originX, originY, originZ, NX, NY, NZ, 2, 2, 2, vertices, elements);
    write_hex_toVTK(fn_mesh_initial.c_str(), vertices, elements);

	int level = 3;
    // Handle local refinement or default processing
    if (!localRefine) {
        write_hex_toVTK(fn_mesh.c_str(), vertices, elements);
        bzmesh3D(path_in); // Generate Bezier mesh info
    } else {
        THS3D(path_in, level); // Perform local refinement
        // THS3D(path_in); // Perform local refinement
    }

    // Partition the mesh for parallel processing
    mpmetis(nProcs, path_in);
}

void gen3Dmesh(int originX, int originY, int originZ, int Nx, int Ny, int Nz, 
               float dx, float dy, float dz, 
               vector<vector<float>>& vertices, vector<vector<int>>& elements) {
	cout << "******************************************************************************\n";
    cout << "                        Generating 3D Structured Initial Mesh                 \n";
	cout << "------------------------------------------------------------------------------\n";
    cout << "Mesh Configuration:\n";
    cout << "  Dimensions (Nx x Ny x Nz):   " << Nx << " x " << Ny << " x " << Nz << "\n";
    cout << "  Origin:                      (" << originX << ", " << originY << ", " << originZ << ")\n";
    cout << "  Spacing (dx, dy, dz):        (" << dx << ", " << dy << ", " << dz << ")\n";
    cout << "-----------------------------------------------------------------------------\n";

    // Clear previous vertices and elements
    vertices.clear();
    elements.clear();

    vector<float> tmp_vtx; // Temporary storage for vertex
    vector<int> tmp_ele;   // Temporary storage for element

    // Generate vertices
    for (int k = 0; k <= Nz; ++k) {
        for (int j = 0; j <= Ny; ++j) {
            for (int i = 0; i <= Nx; ++i) {
                tmp_vtx.clear();
                tmp_vtx.push_back(originX + i * dx); // X-coordinate
                tmp_vtx.push_back(originY + j * dy); // Y-coordinate
                tmp_vtx.push_back(originZ + k * dz); // Z-coordinate
                vertices.push_back(tmp_vtx);
            }
        }
    }

    // Generate elements
    int tl_pt; // Top-left point of the current element
    for (int k = 0; k < Nz; ++k) {
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                tl_pt = (k * (Ny + 1) + j) * (Nx + 1) + i;
                tmp_ele.clear();
                tmp_ele.push_back(tl_pt);                                   // Bottom-left-front
                tmp_ele.push_back(tl_pt + 1);                               // Bottom-right-front
                tmp_ele.push_back(tl_pt + Nx + 2);                          // Bottom-right-back
                tmp_ele.push_back(tl_pt + Nx + 1);                          // Bottom-left-back
                tmp_ele.push_back(tl_pt + (Nx + 1) * (Ny + 1));             // Top-left-front
                tmp_ele.push_back(tl_pt + (Nx + 1) * (Ny + 1) + 1);         // Top-right-front
                tmp_ele.push_back(tl_pt + (Nx + 1) * (Ny + 1) + Nx + 2);    // Top-right-back
                tmp_ele.push_back(tl_pt + (Nx + 1) * (Ny + 1) + Nx + 1);    // Top-left-back
                elements.push_back(tmp_ele);
            }
        }
    }

    cout << "Mesh generation complete. Total vertices: " << vertices.size() 
         << ", Total elements: " << elements.size() << "\n";
	cout << "------------------------------------------------------------------------------\n";
}

// generating 3D bezier mesh using spline_src
void bzmesh3D(string path_in){
	cout << "******************************************************************************\n";
	string spline_cmd_tmd("../spline3D_src/spline " + path_in);
	const char* spline_cmd = spline_cmd_tmd.c_str();
	system(spline_cmd);
} 

// partitioning mesh using mpmetis
void mpmetis(int n_process, string path_in){
	string mpmetis_cmd_tmp("mpmetis " + path_in + "bzmeshinfo.txt " + to_string(n_process));
	const char* mpmetis_cmd = mpmetis_cmd_tmp.c_str();
	system(mpmetis_cmd);
}

void write_hex_toVTK(const char* qs, vector<vector<float>>& vertices, vector<vector<int>>& elements)
{
    FILE* fp;
    fp = fopen(qs, "w");
    int nv, nhex, i, j;

    nv = vertices.size();
    nhex = elements.size();

    fprintf(fp, "# vtk DataFile Version 2.0\n");
    fprintf(fp, "3DmeshGen\n");
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

    fprintf(fp, "POINTS %d float\n", nv);
    for (i = 0; i < nv; i++) {
        fprintf(fp, "%.4f %.4f %.4f\n", vertices[i][0], vertices[i][1], vertices[i][2]);
    }

    fprintf(fp, "\nCELLS %d %d\n", nhex, nhex * 9);

    for (i = 0; i < nhex; i++) {
        fprintf(fp, "8 %d %d %d %d %d %d %d %d\n", elements[i][0], elements[i][1], elements[i][2], elements[i][3],
                elements[i][4], elements[i][5], elements[i][6], elements[i][7]);
    }

    fprintf(fp, "\nCELL_TYPES %d\n", nhex);
    for (i = 0; i < nhex; i++) {
        fprintf(fp, "12\n");
    }

    fclose(fp);
}

void PrintVec2TXT(const vector<float>& v, const string& fn, bool visualization)
{
    ofstream fout(fn);
    if (!fout.is_open()) {
        cerr << "Failed to open file: " << fn << endl;
        return; // Exit if file cannot be opened
    }

    fout << setprecision(2) << fixed;

    if (!visualization) {
        // Print each element on a new line for non-visualization mode
        for (size_t i = 0; i < v.size(); i++) {
            fout << v[i] << endl;
        }
    } else {
        // Visualization mode assumes a square layout
        int sq_sz = static_cast<int>(sqrt(v.size()));
        for (int i = 0; i < sq_sz; i++) {
            for (int j = 0; j < sq_sz; j++) {
                // Print with alignment, ensure spacing for zero and non-zero values
                fout << setw(5);
                if (v[i * sq_sz + j] == 0) {
                    fout << " "; // Use a single space for zero values for better visibility
                } else {
                    fout << v[i * sq_sz + j];
                }
            }
            fout << endl;
        }
    }

    fout.close();
}

void writeVectorToFile(const vector<float>& data, const string& filename, bool binary) {
	ofstream outfile;

	if (binary) {
		outfile.open(filename, ios::out | ios::binary);
	} else {
		outfile.open(filename);
	}

	if (!outfile) {
		cerr << "Error opening file: " << filename << endl;
		return;
	}

	if (binary) {
		outfile.write(reinterpret_cast<const char*>(data.data()), data.size() * sizeof(float));
	} else {
		for (const auto& value : data) {
			outfile << value << " ";
		}
	}

	cout << "Vector successfully written to " << filename << endl;
	outfile.close();
}

vector<float> readVectorFromFile(const string& filename, bool binary) {
	ifstream infile;

	if (binary) {
		infile.open(filename, ios::in | ios::binary);
	} else {
		infile.open(filename);
	}

	if (!infile) {
		cerr << "Error opening file: " << filename << endl;
		return {};
	}

	vector<float> data;

	if (binary) {
		infile.seekg(0, ios::end);
		size_t fileSize = infile.tellg();
		infile.seekg(0, ios::beg);

		data.resize(fileSize / sizeof(float));
		infile.read(reinterpret_cast<char*>(data.data()), fileSize);
	} else {
		float value;

		while (infile >> value) {
			data.push_back(value);
		}
	}

	cout << "Vector successfully read from " << filename << endl;
	infile.close();

	return data;
}

void InitializeSoma(int numNeuron, vector<array<float, 3>>& seed, int& NX, int& NY, int& NZ) {
    // Resize the seed vector to accommodate the number of neurons
    seed.resize(numNeuron);

    // Initialize neuron soma based on the number of neurons
    switch (numNeuron) {
        case 1:
            // // Single neuron case
            // NX = 10;
            // NY = 10;
            // NZ = 10;  // Assumes 3D initialization
            // seed[0] = {10.0f, 10.0f, 10.0f};
            NX = 10;
            NY = 10;
            NZ = 10;  // Assumes 3D initialization
            seed[0] = {10.0f, 10.0f, 10.0f};
            break;

        case 2:
            // Two neurons in a 2D plane
            NX = 140;
            NY = 70;
            NZ = 1;  // Flat plane
            seed[0] = {35.0f, 35.0f, 0.0f};
            seed[1] = {105.0f, 35.0f, 0.0f};
            break;

        case 3:
            // Three neurons in a 2D plane
            NX = 140;
            NY = 130;
            NZ = 1;  // Flat plane
            seed[0] = {35.0f, 35.0f, 0.0f};
            seed[1] = {105.0f, 35.0f, 0.0f};
            seed[2] = {70.0f, 95.0f, 0.0f};
            break;

        case 4:
            // Four neurons in a 2D grid
            NX = 140;
            NY = 140;
            NZ = 1;  // Flat plane
            seed[0] = {35.0f, 35.0f, 0.0f};
            seed[1] = {105.0f, 35.0f, 0.0f};
            seed[2] = {35.0f, 105.0f, 0.0f};
            seed[3] = {105.0f, 105.0f, 0.0f};
            break;

        case 5:
            // Five neurons in a 2D grid with a center neuron
            NX = 140;
            NY = 140;
            NZ = 1;  // Flat plane
            seed[0] = {35.0f, 35.0f, 0.0f};
            seed[1] = {105.0f, 35.0f, 0.0f};
            seed[2] = {35.0f, 105.0f, 0.0f};
            seed[3] = {105.0f, 105.0f, 0.0f};
            seed[4] = {70.0f, 70.0f, 0.0f};
            break;

        default:
            cerr << "Unsupported number of neurons: " << numNeuron << endl;
            throw invalid_argument("Number of neurons must be between 1 and 5.");
    }
}

void ReadMesh(string fn, vector<Vertex3D>& pts, vector<Element3D>& mesh)//need vtk file with point label
{
	string fname(fn), stmp;
	int npts, neles, itmp;
	ifstream fin;
	fin.open(fname);
	if (fin.is_open())
	{
		for (int i = 0; i < 4; i++) getline(fin, stmp);//skip lines
		fin >> stmp >> npts >> stmp;
		pts.resize(npts);
		for (int i = 0; i < npts; i++)
		{
			fin >> pts[i].coor[0] >> pts[i].coor[1] >> pts[i].coor[2];
		}
		getline(fin, stmp);
		fin >> stmp >> neles >> itmp;
		mesh.resize(neles);
		for (int i = 0; i < neles; i++)
		{
			fin >> itmp >> mesh[i].IEN[0] >> mesh[i].IEN[1] >> mesh[i].IEN[2] >> mesh[i].IEN[3] >>
				mesh[i].IEN[4] >> mesh[i].IEN[5] >> mesh[i].IEN[6] >> mesh[i].IEN[7];
			for (int j = 0; j < 8; j++)
			{
				mesh[i].pts[j][0] = pts[mesh[i].IEN[j]].coor[0];
				mesh[i].pts[j][1] = pts[mesh[i].IEN[j]].coor[1];
				mesh[i].pts[j][2] = pts[mesh[i].IEN[j]].coor[2];
			}

		}
		for (int i = 0; i < neles + 5; i++) getline(fin, stmp);//skip lines
		for (int i = 0; i < npts; i++)	fin >> pts[i].label;
		fin.close();
		PetscPrintf(PETSC_COMM_WORLD, "Mesh Loaded!\n");
	}
	else
	{
		PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s!\n", fname.c_str());
	}
}

void ReadControlPoints(string fn, vector<Vertex3D>& pts)
{
	string fname(fn), stmp;
	int npts, neles, itmp;
	ifstream fin;
	fin.open(fname);
	if (fin.is_open())
	{
		for (int i = 0; i < 4; i++) getline(fin, stmp);//skip lines
		fin >> stmp >> npts >> stmp;
		pts.resize(npts);
		for (int i = 0; i < npts; i++)
		{
			fin >> pts[i].coor[0] >> pts[i].coor[1] >> pts[i].coor[2];
		}
		getline(fin, stmp);
		// fin >> stmp >> neles >> itmp;
		// mesh.resize(neles);
		// for (int i = 0; i < neles; i++)
		// {
		// 	fin >> itmp >> mesh[i].IEN[0] >> mesh[i].IEN[1] >> mesh[i].IEN[2] >> mesh[i].IEN[3] >>
		// 		mesh[i].IEN[4] >> mesh[i].IEN[5] >> mesh[i].IEN[6] >> mesh[i].IEN[7];
		// 	for (int j = 0; j < 8; j++)
		// 	{
		// 		mesh[i].pts[j][0] = pts[mesh[i].IEN[j]].coor[0];
		// 		mesh[i].pts[j][1] = pts[mesh[i].IEN[j]].coor[1];
		// 		mesh[i].pts[j][2] = pts[mesh[i].IEN[j]].coor[2];
		// 	}

		// }
		// for (int i = 0; i < neles + 5; i++) getline(fin, stmp);//skip lines
		// for (int i = 0; i < npts; i++)	fin >> pts[i].label;
		fin.close();
		PetscPrintf(PETSC_COMM_WORLD, "Control Points Loaded!\n");
	}
	else
	{
		PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s!\n", fname.c_str());
	}
}

void AssignProcessor(const string& fn, int& n_bzmesh, vector<vector<int>>& ele_process)
{
    // Open the input file
    ifstream fin(fn);
    if (!fin.is_open()) {
        PetscPrintf(PETSC_COMM_WORLD, "Error: Cannot open file %s\n", fn.c_str());
        return;
    }

    // Read and distribute elements to processes
    int tmp;
    n_bzmesh = 0; // Initialize mesh counter
    while (fin >> tmp) {
        ele_process[tmp].push_back(n_bzmesh++); // Increment element counter
        // Skip any non-integer characters (e.g., newlines or spaces)
        fin.ignore(numeric_limits<streamsize>::max(), '\n');
    }

    PetscPrintf(PETSC_COMM_WORLD, "Mesh partition finished. Total elements: %d\n", n_bzmesh);

    // Close the file
    fin.close();
}

// Converts a 3D vector of integers to a 1D vector of floats
vector<float> Convert3DIntTo1DFloatVector(const vector<vector<vector<int>>>& input) 
{
    vector<float> output;
    output.reserve(input.size() * input[0].size() * input[0][0].size()); // Preallocate memory for performance

    for (const auto& matrix : input) {
        for (const auto& row : matrix) {
            for (int value : row) {
                output.emplace_back(static_cast<float>(value)); // Cast int to float and store in the 1D vector
            }
        }
    }

    return output;
}

// Converts a 3D vector of floats to a 1D vector of floats
vector<float> Convert3DFloatTo1DFloatVector(const vector<vector<vector<float>>>& input) 
{
    vector<float> output;
    output.reserve(input.size() * input[0].size() * input[0][0].size()); // Preallocate memory for performance

    for (const auto& matrix : input) {
        for (const auto& row : matrix) {
            for (float value : row) {
                output.emplace_back(value); // Store float value in the 1D vector
            }
        }
    }

    return output;
}

// Function to search for a particular vertex in the vector of Vertex3D
bool SearchVertex(const vector<Vertex3D>& vertices, float targetX, float targetY, float targetZ, int& ind) {
	for (int i = 0; i < vertices.size(); i++) {
		if (vertices[i].coor[0] == targetX && vertices[i].coor[1] == targetY && vertices[i].coor[2] == targetZ) {
			ind = i;
			return true; // Found the vertex (targetX, targetY, targetZ) in the vector
		}
	}
	return false; // Vertex not found in the vector
}

// Function to search for a particular x, y, and z in the vector of Vertex3D
bool SearchPair3D(const vector<Vertex3D> prev_cpts, float targetX, float targetY, float targetZ, int &ind) {
	for (int i = 0; i < prev_cpts.size(); i++) {
		if (prev_cpts[i].coor[0] == targetX && prev_cpts[i].coor[1] == targetY && prev_cpts[i].coor[2] == targetZ) {
			ind = i;
			return true; // Found the triplet (targetX, targetY, targetZ) in the vector
		}
	}
	return false; // Triplet not found in the vector
}

// Function to interpolate values for a new mesh based on coordinates
vector<float> InterpolateValues3D(const vector<Vertex3D>& cpts_initial, const vector<float>& input,
                                      const vector<Vertex3D>& cpts_new) {
	vector<float> output;
	output.resize(cpts_new.size());

	for (int i = 0; i < cpts_new.size(); i++) {
		// cout << i << endl;
		float x = cpts_new[i].coor[0];
		float y = cpts_new[i].coor[1];
		float z = cpts_new[i].coor[2];

		int ind;
		if (SearchVertex(cpts_initial, x, y, z, ind)) {
			// Exact match found, no need for interpolation
			output[i] = input[ind];
			// cout << ind << " ";
		} 
		else {
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

			// cout << ind << " ";
		}
	}
	return output;
}

vector<float> InterpolateVars3D(vector<vector<vector<int>>> input, vector<Vertex3D> cpts_initial, vector<Vertex3D> cpts, int type) 
{   
	vector<float> output;
	output.resize(cpts.size());
	vector<float> tmp = Convert3DIntTo1DFloatVector(input);

	for (int i = 0; i < cpts.size(); i++) {
		float x = cpts[i].coor[0];
		if (abs(remainder(x, 1)) != 0.5) {
			x = round(x);
		}
		float y = cpts[i].coor[1];
		if (abs(remainder(y, 1)) != 0.5) {
			y = round(y);
		}
		float z = cpts[i].coor[2];
		if (abs(remainder(z, 1)) != 0.5) {
			z = round(z);
		}

		int ind;
		if (SearchPair3D(cpts_initial, x, y, z, ind)) {
			output[i] = tmp[ind];
		} 
		else {
			int indDownX, indUpX, indDownY, indUpY, indDownZ, indUpZ;
			float weightDownX, weightUpX, weightDownY, weightUpY, weightDownZ, weightUpZ;

			// Calculate interpolation weights
			weightDownX = abs(x - floor(x));
			weightUpX = 1.0 - weightDownX;
			weightDownY = abs(y - floor(y));
			weightUpY = 1.0 - weightDownY;
			weightDownZ = abs(z - floor(z));
			weightUpZ = 1.0 - weightDownZ;

			// Find neighboring points
			indDownX = SearchNeighbor(cpts_initial, floor(x), round(y), round(z));
			indUpX = SearchNeighbor(cpts_initial, floor(x) + 1, round(y), round(z));
			indDownY = SearchNeighbor(cpts_initial, round(x), floor(y), round(z));
			indUpY = SearchNeighbor(cpts_initial, round(x), floor(y) + 1, round(z));
			indDownZ = SearchNeighbor(cpts_initial, round(x), round(y), floor(z));
			indUpZ = SearchNeighbor(cpts_initial, round(x), round(y), floor(z) + 1);

			// Perform linear interpolation
			float interpolatedValue = 0.0;

			interpolatedValue += weightDownX * weightDownY * weightDownZ * tmp[indDownX];
			interpolatedValue += weightDownX * weightDownY * weightUpZ * tmp[indUpZ];
			interpolatedValue += weightDownX * weightUpY * weightDownZ * tmp[indDownY];
			interpolatedValue += weightDownX * weightUpY * weightUpZ * tmp[indUpY];
			interpolatedValue += weightUpX * weightDownY * weightDownZ * tmp[indDownX];
			interpolatedValue += weightUpX * weightDownY * weightUpZ * tmp[indUpZ];
			interpolatedValue += weightUpX * weightUpY * weightDownZ * tmp[indDownY];
			interpolatedValue += weightUpX * weightUpY * weightUpZ * tmp[indUpY];

			output[i] = interpolatedValue;
		}
	}
	return output;
}

int SearchNeighbor(const vector<Vertex3D>& cpts, float targetX, float targetY, float targetZ) {
	float epsilon = 1e-5;  // A small value to handle floating-point precision issues

	for (int i = 0; i < cpts.size(); i++) {
		float diffX = abs(cpts[i].coor[0] - targetX);
		float diffY = abs(cpts[i].coor[1] - targetY);
		float diffZ = abs(cpts[i].coor[2] - targetZ);

		if (diffX < epsilon && diffY < epsilon && diffZ < epsilon) {
			return i;  // Found the neighbor with matching coordinates
		}
	}

	// If no exact match is found, you may need to handle this case based on your requirements
	// You might consider more sophisticated search algorithms or handle interpolation differently
	return -1;  // Return -1 to indicate that no exact match is found
}

float round5(float value) {
	// Extract the fractional part and the whole part of the value
	float wholePart = floor(value);
	float fractionalPart = value - wholePart;

	// Check if the fractional part is exactly 0.5
	if (fractionalPart == 0.5 || fractionalPart == -0.5) {
		return wholePart; // Round down to 0 for 0.5
	} else {
		return round(value); // Use standard rounding for other cases
	}
}

// Function to perform linear interpolation between two values
float Lerp(float a, float b, float t) {
	return a + t * (b - a);
}

float ComputeDistance(const Vertex3D& a, const Vertex3D& b) {
    float dx = a.coor[0] - b.coor[0];
    float dy = a.coor[1] - b.coor[1];
    float dz = a.coor[2] - b.coor[2];
    return sqrt(dx * dx + dy * dy + dz * dz);
}

void THS3D(const string &path_in, int level) {
// void THS3D(const string &path_in) {
    cout << "******************************************************************************" << endl;
    cout << "Local refinement based on Xiaodong's THS3D code ..." << endl;
    cout << "------------------------------------------------------------------------------" << endl;

    // Construct the command
    string ths3d_cmd = "../THS3D/THS3D " + path_in + " " + to_string(level);
    // string ths3d_cmd = "../THS3D/THS3D " + path_in;

    // Log the command for debugging
    cout << "Executing THS3D command: " << ths3d_cmd << endl;

    // Execute the command
    int ret_code = system(ths3d_cmd.c_str());
    if (ret_code != 0) {
        cerr << "Error: THS3D command failed with return code " << ret_code << endl;
    } else {
        cout << "THS3D completed successfully." << endl;
    }
}