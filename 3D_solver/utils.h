#ifndef UTILS_H
#define UTILS_H

#include <array>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <vector>
#include "BasicDataStructure.h"

using namespace std;

// Utility Functions
void removeFiles(const vector<string>& files);
PetscErrorCode initializeMPI(int& rank, int& nProcs, int argc, char** argv, const char help[]);

// Simulation Setup
void setupSimulationFiles(const int nProcs, const string& path_in, bool localRefine,
                          vector<vector<float>>& vertices, vector<vector<int>>& elements,
                          int NX, int NY, int NZ, int originX, int originY, int originZ);

// Mesh Generation and Manipulation
void gen3Dmesh(int originX, int originY, int originZ, int Nx, int Ny, int Nz, 
               float dx, float dy, float dz, 
               vector<vector<float>>& vertices, vector<vector<int>>& elements);
void bzmesh3D(string path_in);
void mpmetis(int n_process, string path_in);

// Visualization
void write_hex_toVTK(const char* qs, vector<vector<float>>& vertices, vector<vector<int>>& elements);

// Printing and File Operations
void PrintVec2TXT(const std::vector<float>& v, const std::string& fn, bool visualization);
void writeVectorToFile(const std::vector<float>& data, const std::string& filename, bool binary);
std::vector<float> readVectorFromFile(const std::string& filename, bool binary);

// Soma Initialization
void InitializeSoma(int numNeuron, vector<array<float, 3>>& seed, int& NX, int& NY, int& NZ);

// Mesh Processing
void ReadMesh(string fn, vector<Vertex3D>& pts, vector<Element3D>& mesh);
void ReadControlPoints(string fn, vector<Vertex3D>& pts);
void AssignProcessor(const string& fn, int& n_bzmesh, vector<vector<int>>& ele_process);

// 3D Vector Conversion
vector<float> Convert3DIntTo1DFloatVector(const vector<vector<vector<int>>>& input);
vector<float> Convert3DFloatTo1DFloatVector(const vector<vector<vector<float>>>& input);

// Search and Interpolation
bool SearchVertex(const vector<Vertex3D>& vertices, float targetX, float targetY, float targetZ, int& ind);
bool SearchPair3D(const vector<Vertex3D> prev_cpts, float targetX, float targetY, float targetZ, int& ind);
vector<float> InterpolateValues3D(const vector<Vertex3D>& cpts_initial, const vector<float>& input, const vector<Vertex3D>& cpts_new);
vector<float> InterpolateVars3D(vector<vector<vector<int>>> input, vector<Vertex3D> cpts_initial, vector<Vertex3D> cpts, int type);
int SearchNeighbor(const vector<Vertex3D>& cpts, float targetX, float targetY, float targetZ);

// Mathematical Utilities
float round5(float value);
float Lerp(float a, float b, float t);

// Refinement and Clustering
std::vector<float> ComputeRefine(const std::vector<float>& phi, int NX, int NY, int NZ);

// THS3D Integration
void THS3D(const std::string& path_in);

#endif // UTILS_H