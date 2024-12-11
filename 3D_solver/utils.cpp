#include "utils.h"
#include <iostream>
#include <algorithm>
#include "BasicDataStructure.h"
#include <cmath>

void gen3Dmesh(int originX, int originY, int originZ, int Nx, int Ny, int Nz, vector<vector<float>>& vertices, vector<vector<int>>& elements)
{
    cout << "******************************************************************************" << endl;
    cout << "Generating 3D structured initial mesh" << endl;
    cout << "-----------------------------------------------------------------------------" << endl;
    cout << "Nx by Ny by Nz: " << Nx << " x " << Ny << " x " << Nz << " | origin: " << originX << "," << originY << "," << originZ << endl;
    vertices.clear(); elements.clear();
    vector<float> tmp_vtx;
    vector<int> tmp_ele;
    
    for (int k = originZ; k <= (originZ + Nz); k++) {
        for (int j = originY; j <= (originY + Ny); j++) {
            for (int i = originX; i <= (originX + Nx); i++) {
                tmp_vtx.clear();
		// tmp_vtx.push_back((float)i/4);
                // tmp_vtx.push_back((float)j/4);
                // tmp_vtx.push_back((float)k/4);
		tmp_vtx.push_back((float)i);
                tmp_vtx.push_back((float)j);
                tmp_vtx.push_back((float)k);
                vertices.push_back(tmp_vtx);
            }
        }
    }

    int tl_pt;
    for (int k = 0; k < Nz; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                tl_pt = (k * (Ny + 1) + j) * (Nx + 1) + i;
                tmp_ele.clear();
                tmp_ele.push_back(tl_pt);
                tmp_ele.push_back(tl_pt + 1);
                tmp_ele.push_back(tl_pt + Nx + 2);
                tmp_ele.push_back(tl_pt + Nx + 1);
                tmp_ele.push_back(tl_pt + (Nx + 1) * (Ny + 1));
                tmp_ele.push_back(tl_pt + (Nx + 1) * (Ny + 1) + 1);
                tmp_ele.push_back(tl_pt + (Nx + 1) * (Ny + 1) + Nx + 2);
                tmp_ele.push_back(tl_pt + (Nx + 1) * (Ny + 1) + Nx + 1);
                elements.push_back(tmp_ele);
            }
        }
    }
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

// // Creating 2D mesh (incrementing from lo to hi)
// void gen3Dmesh(int originX, int originY, int originZ, int Nx, int Ny, int Nz, vector<vector<float>>& vertices, vector<vector<int>>& elements)
// {
//     std::cout << "******************************************************************************" << std::endl;
//     std::cout << "Generating 3D structured initial mesh" << std::endl;
//     std::cout << "-----------------------------------------------------------------------------" << std::endl;
//     std::cout << "Nx by Ny: " << Nx << " x " << Ny << " z " << Nz << " | origin: " << originX << "," << originY << "," << originY << std::endl;
//     vertices.clear(); elements.clear();
//     vector<float> tmp_vtx;
//     vector<int> tmp_ele;
//     for(int i = 0; i < Nx; i++)
//     {
//         for(int j = 0; j < Ny; j++)
//         {
//             for(int k = 0; k < Nz; k++)
//             {
//                 tmp_vtx.clear();
//                 tmp_vtx.push_back(i);
//                 tmp_vtx.push_back(j);
//                 tmp_vtx.push_back(k);
//                 vertices.push_back(tmp_vtx);
//             }
//         }
//     }

//     int tl_pt;
//     // std::cout << "Testing tl_pt" << std::endl;
//     for(int i = 0; i < Nx-1; i++)
//     {
//         for(int j = 0; j < Ny-1; j++)
//         {
//             for(int k = 0; k < Nz-1; k++)
//             {
//                 tl_pt = i*(Nx)+j*(Ny)+k;
//                 // std::cout<< tl_pt << std::endl;
//                 tmp_ele.clear();
// 		tmp_ele.push_back(tl_pt);
// 		tmp_ele.push_back(tl_pt+1);
// 		tmp_ele.push_back(tl_pt+Ny+2);
// 		tmp_ele.push_back(tl_pt+Ny+1);

// 		tmp_ele.push_back(tl_pt + Nx*Ny);
// 		tmp_ele.push_back(tl_pt+1 + Nx*Ny);
// 		tmp_ele.push_back(tl_pt+Ny+2 + Nx*Ny);
// 		tmp_ele.push_back(tl_pt+Ny+1 + Nx*Ny);

//                 // tmp_ele.push_back(tl_pt);
//                 // tmp_ele.push_back(tl_pt+1);
//                 // tmp_ele.push_back(tl_pt+Nx+1);
//                 // tmp_ele.push_back(tl_pt+Nx);
//                 // tmp_ele.push_back(tl_pt+Ny);
//                 // tmp_ele.push_back(tl_pt+Ny+1);
//                 // tmp_ele.push_back(tl_pt+Ny+Nx+1);
//                 // tmp_ele.push_back(tl_pt+Ny+Nx);
//                 elements.push_back(tmp_ele);
//             }
//         }
//     }
// 
// }

void PrintVec2TXT(const std::vector<float>& v, const std::string& fn, bool visualization)
{
    std::ofstream fout(fn);
    if (!fout.is_open()) {
        std::cerr << "Failed to open file: " << fn << std::endl;
        return; // Exit if file cannot be opened
    }

    fout << std::setprecision(2) << std::fixed;

    if (!visualization) {
        // Print each element on a new line for non-visualization mode
        for (size_t i = 0; i < v.size(); i++) {
            fout << v[i] << std::endl;
        }
    } else {
        // Visualization mode assumes a square layout
        int sq_sz = static_cast<int>(std::sqrt(v.size()));
        for (int i = 0; i < sq_sz; i++) {
            for (int j = 0; j < sq_sz; j++) {
                // Print with alignment, ensure spacing for zero and non-zero values
                fout << std::setw(5);
                if (v[i * sq_sz + j] == 0) {
                    fout << " "; // Use a single space for zero values for better visibility
                } else {
                    fout << v[i * sq_sz + j];
                }
            }
            fout << std::endl;
        }
    }

    fout.close();
}

// // Export hex mesh to vtk for visualization
// void write_hex_toVTK(const char* qs, vector<vector<float>>& vertices, vector<vector<int>>& elements)
// {
// 	FILE* fp;
// 	fp = fopen(qs, "w");
// 	int nv, nhex, i, j;

// 	nv = vertices.size();
// 	nhex = elements.size();

// 	fprintf(fp, "# vtk DataFile Version 2.0\n");
// 	fprintf(fp, "2DmeshGen\n");
// 	fprintf(fp, "ASCII\n");
// 	fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

// 	fprintf(fp, "POINTS %d float\n", nv);
// 	for (i = 0; i < nv; i++) {
// 		fprintf(fp, "%.2f %.2f %.2f\n", vertices[i][0], vertices[i][1], vertices[i][2]);
// 	}

// 	fprintf(fp, "\nCELLS %d %d\n", nhex, nhex * 9);

// 	for (i = 0; i < nhex; i++) {
// 		fprintf(fp, "8 %d %d %d %d %d %d %d %d\n", elements[i][0], elements[i][1], elements[i][2], elements[i][3], elements[i][4], elements[i][5], elements[i][6], elements[i][7]);
// 	}

// 	fprintf(fp, "\nCELL_TYPES %d\n", nhex);
// 	for (i = 0; i < nhex; i++) {
// 		fprintf(fp, "12\n");

// 	}
// 	fclose(fp);
// }

// generating 2D bezier mesh using spline2D_src
void bzmesh2D(string path_in){
    std::cout << "******************************************************************************" << std::endl;

    string spline_cmd_tmd("../spline2D_src/spline -i " + path_in);
    const char* spline_cmd = spline_cmd_tmd.c_str();
    system(spline_cmd);
} 

// generating 3D bezier mesh using spline_src
void bzmesh3D(string path_in){
	std::cout << "******************************************************************************" << std::endl;
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

void THS3D(const std::string &path_in) {
    std::cout << "******************************************************************************" << std::endl;
    std::cout << "Local refinement based on Xiaodong's THS3D code ..." << std::endl;
    std::cout << "******************************************************************************" << std::endl;

    // Construct the command
    std::string ths3d_cmd = "../THS3D/THS3D " + path_in;

    // Log the command for debugging
    std::cout << "Executing THS3D command: " << ths3d_cmd << std::endl;

    // Execute the command
    int ret_code = std::system(ths3d_cmd.c_str());
    if (ret_code != 0) {
        std::cerr << "Error: THS3D command failed with return code " << ret_code << std::endl;
    } else {
        std::cout << "THS3D completed successfully." << std::endl;
    }
}

void InitializeSoma(int numNeuron, vector<array<float, 3>>& seed, int& NX, int& NY, int& NZ) {
    // Resize the seed vector to accommodate the number of neurons
    seed.resize(numNeuron);

    // Initialize neuron soma based on the number of neurons
    switch (numNeuron) {
        case 1:
            // Single neuron case
            NX = 20;
            NY = 20;
            NZ = 20;  // Assumes 3D initialization
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

void AssignProcessor(string fn, int &n_bzmesh, vector<vector<int>> &ele_process)
{
	int tmp;
	int i = 0;
	string fname(fn);
	ifstream fin;
	fin.open(fname, ios::in);
	if (fin.is_open())
	{
		while (!fin.eof() && fin.peek() != EOF)
		{
			fin >> tmp;
			ele_process[tmp].push_back(i);
			i++;
			fin.get();
		}
		n_bzmesh = i;
		//PetscPrintf(PETSC_COMM_WORLD, "Mesh partition finished!\n");
		fin.close();
	}
	else
	{
		PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s!\n", fname.c_str());
	}
}

vector<float> Convert3DIntTo1DFloatVector(const vector<vector<vector<int>>> input) 
{
	vector<float> output;

	for (const auto& matrix : input) {
		for (const auto& row : matrix) {
			for (int value : row) {
				output.push_back(static_cast<float>(value)); // Cast int value to float and add to the 1D vector
			}
		}
	}

	return output;
}

vector<float> Convert3DFloatTo1DFloatVector(const vector<vector<vector<float>>> input) 
{
	vector<float> output;

	for (const auto& matrix : input) {
		for (const auto& row : matrix) {
			for (float value : row) {
				output.push_back(value); // Add float value to the 1D vector
			}
		}
	}

	return output;
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

// vector<float> InterpolateVars3D(vector<vector<vector<int>>> input, vector<Vertex2D> cpts_initial, vector<Vertex2D> cpts, int type) 
// {   
//     vector<float> output;
//     output.resize(cpts.size());
//     vector<float> tmp = Convert3DTo1DFloatVector(input);

//     for (int i = 0; i < cpts.size(); i++) {
//         float x = cpts[i].coor[0];
//         if (abs(remainder(x, 1)) != 0.5) {
//             x = round(x);
//         }
//         float y = cpts[i].coor[1];
//         if (abs(remainder(y, 1)) != 0.5) {
//             y = round(y);
//         }
//         float z = cpts[i].coor[2];
//         if (abs(remainder(z, 1)) != 0.5) {
//             z = round(z);
//         }

//         int ind;
//         if (SearchPair3D(cpts_initial, x, y, z, ind)) {
//             output[i] = tmp[ind];
//         } 
//         else {
//             int indDownX, indUpX, indDownY, indUpY, indDownZ, indUpZ;
//             if ((abs(remainder(x, 1)) == 0.5) && (abs(remainder(y, 1)) != 0.5) && (abs(remainder(z, 1)) != 0.5)) {
//                 if (SearchPair3D(cpts_initial, floorf(x), round(y), round(z), indDownX) &&
//                     SearchPair3D(cpts_initial, floorf(x) + 1, round(y), round(z), indUpX)) {
//                     if (type == 0) {
//                         output[i] = max(tmp[indDownX], tmp[indUpX]);
//                     } else if (type == 1) {
//                         output[i] = (tmp[indDownX] + tmp[indUpX]) / 2;
//                     } else if (type == 2) {
//                         output[i] = 0;
//                     }
//                 }
//                 else {
//                     PetscPrintf(PETSC_COMM_WORLD, "Failed to find pts (ck0)!\n");
//                 }
//             } 
//             else if ((abs(remainder(x, 1)) != 0.5) && (abs(remainder(y, 1)) == 0.5) && (abs(remainder(z, 1)) != 0.5)) {
//                 if (SearchPair3D(cpts_initial, round(x), floorf(y), round(z), indDownY) &&
//                     SearchPair3D(cpts_initial, round(x), floorf(y) + 1, round(z), indUpY)) {
//                     if (type == 0) {
//                         output[i] = max(tmp[indDownY], tmp[indUpY]);
//                     } else if (type == 1) {
//                         output[i] = (tmp[indDownY] + tmp[indUpY]) / 2;
//                     } else if (type == 2) {
//                         output[i] = 0;
//                     }
//                 } 
//                 else {
//                     PetscPrintf(PETSC_COMM_WORLD, "Failed to find pts (ck1)!\n");
//                 }
//             } 
//             else if ((abs(remainder(x, 1)) != 0.5) && (abs(remainder(y, 1)) != 0.5) && (abs(remainder(z, 1)) == 0.5)) {
//                 if (SearchPair3D(cpts_initial, round(x), round(y), floorf(z), indDownZ) &&
//                     SearchPair3D(cpts_initial, round(x), round(y), floorf(z) + 1, indUpZ)) {
//                     if (type == 0) {
//                         output[i] = max(tmp[indDownZ], tmp[indUpZ]);
//                     } else if (type == 1) {
//                         output[i] = (tmp[indDownZ] + tmp[indUpZ]) / 2;
//                     } else if (type == 2) {
//                         output[i] = 0;
//                     }
//                 } 
//                 else {
//                     PetscPrintf(PETSC_COMM_WORLD, "Failed to find pts (ck2)!\n");
//                 }
//             }
//         }
//     }
//     return output;
// }

// // vector<float> InterpolateVars3D(vector<vector<int>> input, vector<Vertex2D> cpts_initial, vector<Vertex2D> cpts, int type) 
// // {	
// // 	vector<float> output;
// // 	output.resize(cpts.size());
// // 	vector<float> tmp = ConvertTo1DFloatVector(input);

// // 	for (int i = 0; i < cpts.size(); i++) {
// // 		// float x = cpts[i].coor[0];
// // 		// float y = cpts[i].coor[1];

// // 		float x = cpts[i].coor[0];
// // 		if (abs(remainder(x,1)) != 0.5) {
// // 			x = round(x);
// // 		}
// // 		float y = cpts[i].coor[1];
// // 		if (abs(remainder(y,1)) != 0.5) {
// // 			y = round(y);
// // 		}

// // 		int ind;
// // 		// if (SearchPair(cpts_initial, round(x), round(y), ind)) {
// // 		if (SearchPair(cpts_initial, x, y, ind)) {
// // 			output[i] = tmp[ind];
// // 		} 
// // 		else {
// // 			int indDown, indUp, indLeft, indRight;
// // 			if ((abs(remainder(x,1)) == 0.5) && (abs(remainder(y,1)) != 0.5)) {
// // 				if (SearchPair(cpts_initial, floorf(x), round(y), indDown) &&
// // 					SearchPair(cpts_initial, floorf(x)+1, round(y), indUp)) {
// // 					if (type == 0) {
// // 						output[i] = max(tmp[indDown], tmp[indUp]);
// // 					} else if (type == 1) {
// // 						output[i] = (tmp[indDown] + tmp[indUp])/2;
// // 					} else if (type == 2) {
// // 						output[i] = 0;
// // 					}
// // 				} else {
// // 					PetscPrintf(PETSC_COMM_WORLD, "Failed to find pts (ck0)!\n");
// // 				}
// // 			} else if ((abs(remainder(x,1)) != 0.5) && (abs(remainder(y,1)) == 0.5)) {
// // 				if (SearchPair(cpts_initial, round(x), floorf(y), indLeft) &&
// // 					SearchPair(cpts_initial, round(x), floorf(y)+1, indRight)) {
// // 					if (type == 0) {
// // 						output[i] = max(tmp[indLeft], tmp[indRight]);
// // 					} else if (type == 1) {
// // 						output[i] = (tmp[indLeft] + tmp[indRight])/2;
// // 					} else if (type == 2) {
// // 						output[i] = 0;
// // 					}
// // 				} else {
// // 					PetscPrintf(PETSC_COMM_WORLD, "Failed to find pts (ck1)!\n");
// // 				}
// // 			} else if ((abs(remainder(x,1)) == 0.5) && (abs(remainder(y,1)) == 0.5)) {
// // 				if (SearchPair(cpts_initial, floorf(x), floorf(y), indDown) &&
// // 					SearchPair(cpts_initial, floorf(x)+1, floorf(y), indUp) &&
// // 					SearchPair(cpts_initial, floorf(x), floorf(y), indLeft) &&
// // 					SearchPair(cpts_initial, floorf(x), floorf(y)+1, indRight)) {
// // 					if (type == 0) {
// // 						output[i] = max(max(tmp[indDown], tmp[indUp]), max(tmp[indLeft], tmp[indRight]));
// // 					} else if (type == 1) {
// // 						output[i] = (tmp[indDown] + tmp[indUp] + tmp[indLeft] + tmp[indRight])/4;
// // 					} else if (type == 2) {
// // 						output[i] = 0;
// // 					}
// // 				} else {
// // 					PetscPrintf(PETSC_COMM_WORLD, "Failed to find pts (ck2)!\n");
// // 				}
// // 			}
// // 		}			
// // 	}
// // 	return output;
// // }

// vector<float> InterpolateVars3D(vector<float> input, vector<Vertex2D> cpts_initial, vector<Vertex2D> cpts, int type) 
// {	
// 	vector<float> output;
// 	output.resize(cpts.size());

// 	for (int i = 0; i < cpts_initial.size(); i++) {
// 		if (abs(remainder(cpts_initial[i].coor[0],1)) != 0.5)
// 			cpts_initial[i].coor[0] = round(cpts_initial[i].coor[0]);
// 		if (abs(remainder(cpts_initial[i].coor[1],1)) != 0.5)
// 			cpts_initial[i].coor[1] = round(cpts_initial[i].coor[1]);
// 	}
	
// 	for (int i = 0; i < cpts.size(); i++) {
// 		// float x = cpts[i].coor[0];
// 		// float y = cpts[i].coor[1];

// 		float x = cpts[i].coor[0];
// 		if (abs(remainder(x,1)) != 0.5) {
// 			x = round(x);
// 		}
// 		float y = cpts[i].coor[1];
// 		if (abs(remainder(y,1)) != 0.5) {
// 			y = round(y);
// 		}

// 		int ind;
// 		// if (SearchPair(cpts_initial, round(x), round(y), ind)) {
// 		if (SearchPair(cpts_initial, x, y, ind)) {
// 			output[i] = input[ind];
// 		} 
// 		else {
// 			int indDown, indUp, indLeft, indRight;
// 			if ((abs(remainder(x,1)) == 0.5) && (abs(remainder(y,1)) != 0.5)) {
// 				if (SearchPair(cpts_initial, floorf(x), round(y), indDown) &&
// 					SearchPair(cpts_initial, floorf(x)+1, round(y), indUp)) {
// 					if (type == 0) {
// 						output[i] = max(input[indDown], input[indUp]);
// 					} else if (type == 1) {
// 						output[i] = (input[indDown] + input[indUp])/2;
// 					} else if (type == 2) {
// 						output[i] = 0;
// 					}
// 				} else {
// 					PetscPrintf(PETSC_COMM_WORLD, "Failed to find pts (ck0)!\n");
// 				}
// 			} else if ((abs(remainder(x,1)) != 0.5) && (abs(remainder(y,1)) == 0.5)) {
// 				if (SearchPair(cpts_initial, round(x), floorf(y), indLeft) &&
// 					SearchPair(cpts_initial, round(x), floorf(y)+1, indRight)) {
// 					if (type == 0) {
// 						output[i] = max(input[indLeft], input[indRight]);
// 					} else if (type == 1) {
// 						output[i] = (input[indLeft] + input[indRight])/2;
// 					} else if (type == 2) {
// 						output[i] = 0;
// 					}
// 				} else {
// 					PetscPrintf(PETSC_COMM_WORLD, "Failed to find pts (ck1)!\n");
// 				}
// 			} else if ((abs(remainder(x,1)) == 0.5) && (abs(remainder(y,1)) == 0.5)) {
// 				if (SearchPair(cpts_initial, floorf(x), floorf(y), indDown) &&
// 					SearchPair(cpts_initial, floorf(x)+1, floorf(y), indUp) &&
// 					SearchPair(cpts_initial, floorf(x), floorf(y), indLeft) &&
// 					SearchPair(cpts_initial, floorf(x), floorf(y)+1, indRight)) {
// 					if (type == 0) {
// 						output[i] = max(max(input[indDown], input[indUp]), max(input[indLeft], input[indRight]));
// 					} else if (type == 1) {
// 						output[i] = (input[indDown] + input[indUp] + input[indLeft] + input[indRight])/4;
// 					} else if (type == 2) {
// 						output[i] = 0;
// 					}
// 				} else {
// 					PetscPrintf(PETSC_COMM_WORLD, "Failed to find pts (ck2)!\n");
// 				}
// 			}
// 		}			
// 	}
// 	return output;
// }

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

// Function to perform linear interpolation between two values
float Lerp(float a, float b, float t) {
	return a + t * (b - a);
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

// Function to interpolate values for a new mesh based on coordinates
vector<float> InterpolateValues3D(const vector<Vertex3D>& cpts_initial, const vector<float>& input,
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

// Function to compute the average of surrounding points
std::vector<float> ComputeRefine(const std::vector<float>& phi, int NX, int NY, int NZ) 
{
	std::vector<float> ele_refine(NX * NY * NZ, 0.0);
	float maxPhi(0), phi_average;
	for (int i = 0; i < phi.size(); i++) {
		maxPhi = max(maxPhi, phi[i]);
	}
	
	// std::cout << maxPhi << std::endl;
	for (int i = 0; i < NX; i++) {
		for (int j = 0; j < NY; j++) {
			for (int k = 0; k < NZ; k++) {

				int index_in = i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k;
				int index_out = i * NY * NZ + j * NZ + k;

				// Compute the average of surrounding points
				float phi_average = (phi[index_in - 1] + phi[index_in + 1] +
							phi[index_in - (NZ + 1)] + phi[index_in + (NZ + 1)] +
							phi[index_in - (NY + 1) * (NZ + 1)] + phi[index_in + (NY + 1) * (NZ + 1)]) / 6.0;

				// cout << 'check' << endl;
				if ((phi_average < (0.5 * maxPhi)) && (phi_average > (0.001 * maxPhi))) {
					ele_refine[index_out] = 1;
					// cout << 'check' << endl;
				} else {
					ele_refine[index_out] = 0;
				}

				// ele_refine.push_back(0);
			}
		}
	}

	return ele_refine;
}

// // Function to reshape the input 3D vector to a new size
// std::vector<float> Reshape3DGrid(const std::vector<float>& phi, int NX, int NY, int NZ)
// {
//     // Check if the original size matches the expected size
//     if (phi.size() != (NX + 1) * (NY + 1) * (NZ + 1)) {
//         std::cerr << "Error: Incorrect input size!" << std::endl;
//         return std::vector<float>();
//     }

//     std::vector<float> phi_out;
//     phi_out.reserve(NX * NY * NZ);

//     for (int i = 0; i < NX; ++i) {
//         for (int j = 0; j < NY; ++j) {
//             for (int k = 0; k < NZ; ++k) {
//                 int index = i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k;
//                 phi_out.push_back(phi[index]);
//             }
//         }
//     }

//     return phi_out;
// }

void writeVectorToFile(const std::vector<float>& data, const std::string& filename, bool binary) {
	std::ofstream outfile;

	if (binary) {
		outfile.open(filename, std::ios::out | std::ios::binary);
	} else {
		outfile.open(filename);
	}

	if (!outfile) {
		std::cerr << "Error opening file: " << filename << std::endl;
		return;
	}

	if (binary) {
		outfile.write(reinterpret_cast<const char*>(data.data()), data.size() * sizeof(float));
	} else {
		for (const auto& value : data) {
			outfile << value << " ";
		}
	}

	std::cout << "Vector successfully written to " << filename << std::endl;
	// std::cout << "ckck0" << std::endl;
	outfile.close();
	// std::cout << "ckck1" << std::endl;
}

std::vector<float> readVectorFromFile(const std::string& filename, bool binary) {
	std::ifstream infile;

	if (binary) {
		infile.open(filename, std::ios::in | std::ios::binary);
	} else {
		infile.open(filename);
	}

	if (!infile) {
		std::cerr << "Error opening file: " << filename << std::endl;
		return {};
	}

	std::vector<float> data;

	if (binary) {
		infile.seekg(0, std::ios::end);
		size_t fileSize = infile.tellg();
		infile.seekg(0, std::ios::beg);

		data.resize(fileSize / sizeof(float));
		infile.read(reinterpret_cast<char*>(data.data()), fileSize);
	} else {
		float value;

		while (infile >> value) {
			data.push_back(value);
		}
	}

	std::cout << "Vector successfully read from " << filename << std::endl;
	infile.close();

	return data;
}