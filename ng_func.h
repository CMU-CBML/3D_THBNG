#ifndef NG_FUNC_H    // To make sure you don't declare the function more than once by including the header multiple times.
#define NG_FUNC_H

#include <iostream>     /* cout */
#include <vector>
#include <math.h>       /* tanh */
#include <fstream>

#include <Eigen/Dense>

// using namespace std;

// Setup domain Nx Ny based on number of neurons
void neuron_domain_setup(int numNeuron, int& Nx, int& Ny);

// Creating knot vector (incrementing from lo to hi)
void gen_knotvector(Eigen::VectorXd& knotvector, int lo, int hi, int N, int order);
// void gen_knotvector(vector<float>& knotvector, int lo, int hi, int N, int order);

// Printing the vector
void printVector1D(Eigen::VectorXd& vec);
void printMatrix2D(Eigen::MatrixXd& vec);
void printArray2TXT(Eigen::MatrixXd& v, std::string fn);

// Initialize domains based on given parameters
void initialize_neurite_growth(Eigen::MatrixXd& phi, Eigen::MatrixXd& conct, Eigen::VectorXd& seed_x, Eigen::VectorXd& seed_y, int seed_radius, int lenu, int lenv, int numNeuron);

// Find knot span
int FindSpan(int n, int p, float u, Eigen::VectorXd U);

// Calculate derivatives
void dersbasisfuns(int i, int p, int order_deriv, float u, Eigen::VectorXd U, Eigen::MatrixXd& ders);

// Calculate derivatives based on collocation method
void collocationDers(Eigen::VectorXd knotvectorU,int p, Eigen::VectorXd knotvectorV, int q, int order_deriv, std::vector<Eigen::MatrixXd>& cm);

// Creating a matrix consisting of val
Eigen::MatrixXd matInit(int rows, int cols, float val);

// Creating a vector of zero matrices/arraies
std::vector<Eigen::MatrixXd> vectorMatInit(int n, int rows, int cols);
std::vector<Eigen::ArrayXd> arrayMatInit(int n, int rows, int cols);

// Solve for x from Ax = b
Eigen::MatrixXd linSol(Eigen::MatrixXd A, Eigen::MatrixXd b);

// calculate element-wise atan2, cos, sin, atan
Eigen::MatrixXd getAtan2(Eigen::MatrixXd A, Eigen::MatrixXd b, int lenu, int lenv);
Eigen::MatrixXd getCos(Eigen::MatrixXd input, int lenu, int lenv);
Eigen::MatrixXd getSin(Eigen::MatrixXd input, int lenu, int lenv);
Eigen::MatrixXd getAtan(Eigen::MatrixXd input, int lenu, int lenv);

// calculates epsilon and aap (a*a') based on phi, theta, NuN1v, and N1uNv.
void getEpsilonAndAap(std::vector<Eigen::MatrixXd>& output, float epsilonb, float delta, Eigen::MatrixXd phi, Eigen::MatrixXd xtheta, std::vector<Eigen::MatrixXd>& cm, int lenu, int lenv);

void updateRgSg(Eigen::ArrayXd& rMat, Eigen::ArrayXd& sMat, Eigen::ArrayXd nnT, int lenu, int lenv);

Eigen::ArrayXd regular_Heiviside_fun(Eigen::ArrayXd input, int lenu, int lenv);

void stiffMatSetupBCID(Eigen::MatrixXd& coll_Lhs, Eigen::MatrixXd& coll_Rhs, Eigen::MatrixXd bcid, Eigen::MatrixXd N, int lenu, int lenv);

#endif