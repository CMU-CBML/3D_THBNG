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

void PrintVec2TXT(const std::vector<float>& v, const std::string& fn, bool visualization);
// void PrintVec2TXT(vector<float>& v, string fn, bool visualization); // print out to commandline for debugging

// // generating 3D bezier mesh using spline_src
void bzmesh3D(string path_in);

// partitioning mesh using mpmetis
void mpmetis(int n_process, string path_in);

void THS3D(const std::string &path_in);

void InitializeSoma(int numNeuron, vector<array<float, 3>> &seed, int &NX, int &NY, int &NZ);

void ReadMesh(string fn, vector<Vertex3D>& pts, vector<Element3D>& mesh);

void ReadControlPoints(string fn, vector<Vertex3D>& pts);

void AssignProcessor(string fn, int &n_bzmesh, vector<vector<int>> &ele_process);

vector<float> Convert3DIntTo1DFloatVector(const vector<vector<vector<int>>> input);
vector<float> Convert3DFloatTo1DFloatVector(const vector<vector<vector<float>>> input);

bool SearchPair3D(const vector<Vertex3D> prev_cpts, float targetX, float targetY, float targetZ, int &ind);
// vector<float> InterpolateVars3D(vector<vector<vector<int>>> input, vector<Vertex2D> cpts_initial, vector<Vertex2D> cpts, int type);
// // vector<float> InterpolateVars3D(vector<vector<int>> input, vector<Vertex2D> cpts_initial, vector<Vertex2D> cpts, int type);
// vector<float> InterpolateVars3D(vector<float> input, vector<Vertex2D> cpts_initial, vector<Vertex2D> cpts, int type);

vector<float> InterpolateVars3D(vector<vector<vector<int>>> input, vector<Vertex3D> cpts_initial, vector<Vertex3D> cpts, int type);
int SearchNeighbor(const vector<Vertex3D>& cpts, float targetX, float targetY, float targetZ);


// Function to perform linear interpolation between two values
float Lerp(float a, float b, float t);
// Function to search for a particular vertex in the vector of Vertex3D
bool SearchVertex(const vector<Vertex3D>& vertices, float targetX, float targetY, float targetZ, int& ind);
vector<float> InterpolateValues3D(const vector<Vertex3D>& cpts_initial, const vector<float>& input, const vector<Vertex3D>& cpts_new);

std::vector<float> ComputeRefine(const std::vector<float>& phi, int NX, int NY, int NZ);

void writeVectorToFile(const std::vector<float>& data, const std::string& filename, bool binary);
std::vector<float> readVectorFromFile(const std::string& filename, bool binary);

#endif