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
                tmp_vtx.push_back(i);
                tmp_vtx.push_back(j);
                tmp_vtx.push_back(k);
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

void PrintVec2TXT(std::vector<float>& v, std::string fn, bool visualization)
{
	std::ofstream fout;

	fout.open(fn);

	fout << std::setprecision(2) << std::fixed;

	int v_size = v.size();
	if (visualization == 0) {
		for (int i = 0; i < v.size(); i++) {
			fout << v[i] << std::endl;
		}
	} else {
		int sq_sz = (int)sqrt(v_size);
		int ind = 0;
		for (int i = 0; i < sq_sz; i++) {
			for (int j = 0; j < sq_sz; j++) {
				if (v[ind] == 0) {
					fout << "     ";
				} else {
					fout << v[ind] << " ";
				}
				ind += 1;
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
	string spline_cmd_tmd("../spline3D_src/spline -i " + path_in);
	const char* spline_cmd = spline_cmd_tmd.c_str();
	system(spline_cmd);
} 

// partitioning mesh using mpmetis
void mpmetis(int n_process, string path_in){
	string mpmetis_cmd_tmp("mpmetis " + path_in + "bzmeshinfo.txt " + to_string(n_process));
	const char* mpmetis_cmd = mpmetis_cmd_tmp.c_str();
	system(mpmetis_cmd);
}

// // partitioning mesh using mpmetis
// void THS3D(string path_in, vector<int> rfid, vector<int> rftype){
// 	std::cout << "******************************************************************************" << std::endl;
// 	std::cout << "Local refinement based on Xiaodong's THS3D code ... " << std::endl;
// 	std::cout << "  - see: Truncated T-splines: Fundamentals and methods (2017)" << std::endl << std::endl;
// 	std::cout << "-----------------------------------------------------------------------------" << std::endl;
// 	std::cout << "Calling command | input mesh directory | refine ID | refine element type" << std::endl << std::endl;
// 	string ths2d_cmd_tmp("../THS3D/THS3D " + path_in + " ");
// 	for (int i = 0; i < rfid.size(); i++) {
// 		ths2d_cmd_tmp = ths2d_cmd_tmp + std::to_string(rfid[i]) + " ";
// 	}
// 	for (int i = 0; i < rftype.size(); i++) {
// 		ths2d_cmd_tmp = ths2d_cmd_tmp + std::to_string(rftype[i]) + " ";
// 	}
// 	std::cout << ths2d_cmd_tmp << std::endl;
// 	const char* ths2d_cmd = ths2d_cmd_tmp.c_str();
// 	system(ths2d_cmd);
// }

//
void InitializeSoma(int numNeuron, vector<array<int, 3>> &seed, int &NX, int &NY, int &NZ){
    seed.resize(numNeuron);
    // 2D neuron soma initialization
    switch (numNeuron) {
        case 1:
            NX = 40;
            NY = 40;
	    NZ = 40;
            seed[0][0] = 20;	seed[0][1] = 20;	seed[0][2] = 20;
            break;
        case 2:
            NX = 140;
            NY = 70;
            seed[0][0] = 35;    seed[0][1] = 35;
            seed[1][0] = 105;   seed[1][1] = 35;
            break;
        case 3:
            NX = 140;
            NY = 130;
            seed[0][0] = 35;    seed[0][1] = 35;
            seed[1][0] = 105;   seed[1][1] = 35;
            seed[2][0] = 70;    seed[2][1] = 95;
            break;
        case 4:
            NX = 140;
            NY = 140;
            seed[0][0] = 35;    seed[0][1] = 35;
            seed[1][0] = 105;   seed[1][1] = 35;
            seed[2][0] = 35;    seed[2][1] = 105;
            seed[3][0] = 105;   seed[3][1] = 105;
            break;
        case 5:
            NX = 140;
            NY = 140;
            seed[0][0] = 35;    seed[0][1] = 35;
            seed[1][0] = 105;   seed[1][1] = 35;
            seed[2][0] = 35;    seed[2][1] = 105;
            seed[3][0] = 105;   seed[3][1] = 105;
            seed[4][0] = 70;    seed[4][1] = 70;
            break;
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