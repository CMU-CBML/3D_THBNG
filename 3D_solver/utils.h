#ifndef UTILS_H
#define UTILS_H

#include <array>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <vector>
#include "BasicDataStructure.h"

using namespace std;

// Creating 3D mesh (incrementing from lo to hi)
void gen3Dmesh(int originX, int originY, int originZ, int Nx, int Ny, int Nz, vector<vector<float>>& vertices, vector<vector<int>>& elements);

// Export hex mesh to vtk for visualization
void write_hex_toVTK(const char* qs, vector<vector<float>>& vertices, vector<vector<int>>& elements);

void PrintVec2TXT(vector<float>& v, string fn, bool visualization); // print out to commandline for debugging

// // generating 3D bezier mesh using spline_src
void bzmesh3D(string path_in);

// partitioning mesh using mpmetis
void mpmetis(int n_process, string path_in);

void THS3D(string path_in, vector<int> rfid, vector<int> rftype);

void InitializeSoma(int numNeuron, vector<array<int, 3>> &seed, int &NX, int &NY, int &NZ);

void ReadMesh(string fn, vector<Vertex3D>& pts, vector<Element3D>& mesh);

void AssignProcessor(string fn, int &n_bzmesh, vector<vector<int>> &ele_process);

#endif