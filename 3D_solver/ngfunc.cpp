#include <iostream>     /* cout */
#include <vector>
#include <math.h>       /* tanh */
#include <fstream>
#include <iomanip>     /* setprecision */

#include "ngfunc.h"
#include <Eigen/Dense>

// using namespace std;

NGfunc::NGfunc()
{
}

void NGfunc::neuron_domain_setup(int numNeuron, int& Nx, int& Ny)
{
    switch(numNeuron) {
        case 1: // 1 neuron
            Nx = 60;
            Ny = 60;   
            break;  
        case 2: // 2 neurons
            Nx = 150;
            Ny = 150;   
            break;
        case 3: // 3 neurons
            Nx = 170;
            Ny = 170;   
            break;
        case 4: // 4 neurons
            Nx = 200;
            Ny = 200;   
            break;
        case 5: // 5 neurons
            Nx = 220;
            Ny = 220;   
            break;
        case 6: // 6 neurons
            Nx = 260;
            Ny = 260;   
            break;
        case 7: // 7 neurons
            Nx = 260;
            Ny = 260;   
            break;
    }
    std::cout << "Domain setup test!" << std::endl;
}

// Creating knot vector (incrementing from lo to hi)
void NGfunc::gen_knotvector(Eigen::VectorXd& knotvector, int lo, int hi, int N, int order)
{
    float incr = (float)(hi - lo) / (float)(N-1);

    for(int i = 0; i < N; i++)
    {
        knotvector(i+3) = lo + i * incr;
    }
    
    // to be updated based on q and p
    for(int i = 0; i < order; i++)
    {
        knotvector(0+i) = 0;
        knotvector(N+i+3) = hi;
    }
}

// Printing the vector
void NGfunc::printVector1D(Eigen::VectorXd& vec)
{
    std::cout << "**********************************************************" << std::endl;
    std::cout << "Vector testing | Printing vector..." << std::endl;
    std::cout << "Vector size: " << vec.rows() << std::endl;
    std::cout << "Vector content:";
    for(int i = 0; i < vec.rows(); i++)
    {
        std::cout << vec(i) << " ";
    }
    std::cout << std::endl;
}

void NGfunc::printMatrix2D(Eigen::MatrixXd& vec)
{
    std::cout << "**********************************************************" << std::endl;
    std::cout << "Matrix testing | Printing Matrix..." << std::endl;
    std::cout << "Matrix size: " << vec.size() << std::endl;
    std::cout << "Matrix content:";
    for(int i = 0; i < vec.rows(); i++)
    {
        for(int j = 0; j < vec.cols(); j++)
        {
            std::cout << vec(i,j) << " ";
        }
    }
    std::cout << std::endl;
}

void NGfunc::printArray2TXT(Eigen::MatrixXd& v, std::string fn)
{
    std::ofstream fout;
    
	fout.open(fn);

    fout << std::setprecision(2) << std::fixed;

	for (int i = 0; i < v.rows(); i++)
    {
        // fout << i << "\t";
        for (int j = 0; j < v.cols(); j++)
        {
            // fout << v(i,j);
            if (v(i,j) == 0)
            // if (abs(v(i,j)) <= 1e-1)
            {
                fout << "     ";
            } else{
                fout << v(i,j) << " ";
            }

            // if (v(i,j) == 0)
            // {
            //     fout << "  ";
            // } else{
            //     fout << "##";
            // }
        }
        fout << std::endl;
    }

	fout.close();
}

void NGfunc::initialize_neurite_growth(Eigen::MatrixXd& phi, Eigen::MatrixXd& conct, Eigen::VectorXd& seed_x, Eigen::VectorXd& seed_y, int seed_radius, int lenu, int lenv, int numNeuron)
{        
    if (numNeuron == 1)
    {
        seed_x << 0;
        seed_y << 0;
    } else if (numNeuron == 2)
    {
        seed_x << 45, -45;
        seed_y << 45, -45;
    } else if (numNeuron == 3)
    {
        seed_x << -60,60,60;
        seed_y << 0,40,-40;
    } else if (numNeuron == 4)
    {        
        seed_x << -60,60,60,-60;
        seed_y << -60,60,-60,60;
    } else if (numNeuron == 5)
    {
        seed_x << -75,75,75,-75,0;
        seed_y << -75,75,-75,75,0;
    } else if (numNeuron == 6)
    {
        seed_x << -85,-10,85,-85,10,85;
        seed_y << -15,-85,-85,85,85,15;
    } else if (numNeuron == 7)
    {
        seed_x << -95,-45,45,-45,45,95,0;
        seed_y << 0,-95,-95,95,95,0,0;
    }
 
    int seed = pow(seed_radius,2);
    for (int i = 0; i < lenu; ++i)
    {
        // phi.push_back(vector<float>());
        // conct.push_back(vector<float>());
        for (int j = 0; j < lenv; ++j)
        {
            if (((i-(float)lenu/2)*(i-(float)lenu/2)+(j-(float)lenv/2)*(j-(float)lenv/2)) < seed)
            {
                float r = sqrt((i-lenu/2)*(i-lenu/2)+(j-lenv/2)*(j-lenv/2));
                for (int l = 0; l < numNeuron; ++l)
                {
                    phi(i,j) = 1;
                    conct(i,j) = 0.5+0.5*tanh((sqrt(seed)-r)/2);
                }
            }
        }
    }

    // printArray2TXT(phi, "./phi_initial.txt");
    // printArray2TXT(conct, "./conct_initial.txt");
}

// Find knot span
int NGfunc::FindSpan(int n, int p, float u, Eigen::VectorXd U)
{
    // int knotSpanIndex;
    if (u == U(n+1))
    {
        return(n);
    }
    int low = p;
    int high = n+1;
    int mid = floor((float)(low + high)/2);
    while (u <U(mid) || u >= U(mid+1) )
    { 
        if( u < U(mid+1))
        // while (u <U[mid] || u >= U[mid+1] )
        //    if( u < U[mid])
        {
            high = mid;
        } else {
            low = mid;
        }
        mid = floor((float)(low+high)/2);
    }
    return(mid);
}

// Compute nonzero basis functions and their derivatives.
void NGfunc::dersbasisfuns(int i, int p, int order_deriv, float u, Eigen::VectorXd U, Eigen::MatrixXd& ders)
{
    Eigen::VectorXd left = Eigen::VectorXd::Zero(p+1);
    Eigen::VectorXd right = Eigen::VectorXd::Zero(p+1);
    Eigen::MatrixXd ndu = Eigen::MatrixXd::Zero(p+1, p+1);

    ndu(0,0) = 1;
    float saved, temp;
    for (int j = 1; j <= p; j++)
    {
        left(j) = u - U(i+1-j);
        right(j) = U(i+j) - u;
        saved = 0;
        for (int r = 0; r < j; r++)
        {    
            ndu(j,r) = right(r+1) + left(j-r);
            temp = ndu(r,j-1)/ndu(j,r);

            ndu(r,j) = saved + right(r+1)*temp;
            saved = left(j-r)*temp;
        } 
        ndu(j,j) = saved;
    }

    for (int j = 0; j <= p; j++) // load basis functions
    {    
        ders(0,j) = ndu(j,p);
    }

    // compute derivatives
    int s1, s2, d, rk, pk, j1, j2, j;
    Eigen::MatrixXd a  = Eigen::MatrixXd::Zero(order_deriv+1, p+1);
    // gen2Dzeros_float(a, order_deriv+1, p+1);
    for (int r = 0; r <= p; r++) // loop over function index
    {
        s1 = 0; s2 = 1; // alternate rows in array a
        a(0,0) = 1;
        for (int k = 1; k <= order_deriv; k++) // loop to compute kth derivative
        {
            d = 0;
            rk = r-k;
            pk = p-k;
            
            if(r >= k)
            {
                a(s2,0) = a(s1,0)/ndu(pk+1,rk);
                d = a(s2,0)*ndu(rk,pk);     
            }

            if(rk >= -1)
            {
                j1 = 1;
            } else {
                j1 = -rk;
            }
            
            if((r-1) <= pk)
            {
                j2 = k-1;
            } else {
                j2 = p-r;
            }
            
            for (int j = j1; j <= j2; j++)
            {   
                a(s2,j) =(a(s1,j) - a(s1,j))/ndu(pk+1,rk+j);
                d += a(s2,j)*ndu(rk+j,pk);
            }
            
            if(r <= pk)
            {   
                a(s2,k) = -a(s1,k)/ndu(pk+1,r);
                d += a(s2,k)*ndu(r,pk);
            }
            
            ders(k,r) = d;
            j = s1; s1 = s2; s2 = j; // switch rows
        }
    }

    // multiply through by the correct factors;
    float r = p;
    for (int k = 1; k <= order_deriv; k++)
    {   
        for (int j = 0; j <= p; j++)
        {            
            ders(k,j) *= r;
        }
        r *= (p-k);
    }

}

void NGfunc::collocationDers(Eigen::VectorXd knotvectorU,int p, Eigen::VectorXd knotvectorV, int q, int order_deriv, std::vector<Eigen::MatrixXd>& cm)
{

    int lenu = knotvectorU.size()-2*(p-1);
    int lenv = knotvectorV.size()-2*(q-1);
    
    Eigen::MatrixXd coll_p = Eigen::MatrixXd::Zero(lenu*lenv, 2);

    float coordx, coordy;

    printVector1D(knotvectorU);

    int k = 0;
    for (int i = 0; i < knotvectorU.rows()-p-1; ++i)
    {
        coordx = 0;
        for (int l = 0; l < p; ++l)
        {
            coordx = coordx + knotvectorU(i+l+1);
        } 
        coordx = coordx/p;

        for (int j = 0; j < knotvectorV.rows()-q-1; ++j)
        {
            coordy = 0;
            for (int l = 0; l < q; ++l)
            {
                coordy = coordy + knotvectorV(j+l+1);
            } 
            coordy = coordy/q;
            
            coll_p(k,0) = coordx;
            coll_p(k,1) = coordy;
            k += 1;
        }
    }
         
    int size_collpts = lenu;
    int nobu = knotvectorU.rows()-p-2; 
    int nobv = knotvectorV.rows()-q-2;

    Eigen::VectorXd knotSpanU = Eigen::VectorXd::Zero(pow(size_collpts,2));
    Eigen::VectorXd knotSpanV = Eigen::VectorXd::Zero(pow(size_collpts,2));
    Eigen::MatrixXd ksU = Eigen::MatrixXd::Zero(pow(size_collpts,2),p+1);
    Eigen::MatrixXd ksV = Eigen::MatrixXd::Zero(pow(size_collpts,2),p+1);

    std::vector<Eigen::MatrixXd> dersU, dersV;
    for (int i = 0; i < pow(size_collpts,2); ++i)
    {
        dersU.push_back(Eigen::MatrixXd::Zero(3, p+1));
        dersV.push_back(Eigen::MatrixXd::Zero(3, p+1));
    }

    // vector<float> Nu, Nv, N1u, N1v, N2u, N2v;
    Eigen::VectorXd Nu = Eigen::VectorXd::Zero(p+1);
    Eigen::VectorXd Nv = Eigen::VectorXd::Zero(q+1);
    Eigen::VectorXd N1u = Eigen::VectorXd::Zero(p+1);
    Eigen::VectorXd N1v = Eigen::VectorXd::Zero(q+1);
    Eigen::VectorXd N2u = Eigen::VectorXd::Zero(p+1);
    Eigen::VectorXd N2v = Eigen::VectorXd::Zero(q+1);

    float u,v;
    int ind;
    k = 0;
    for (int i = 0; i < size_collpts; ++i)
    {
        for (int j = 0; j < size_collpts; ++j)
        {
            // if (j==0)
            // {                
            //     std::cout << coll_p(k,0) << std::endl;
            // }            
            k = i*size_collpts+j;
            // extract u,v position
            u = coll_p(k,0);
            v = coll_p(k,1);

            // calculating basis function, its 1st & 2nd derivatives based on
            // collocation points and knot vector
            knotSpanU(k) = FindSpan(nobu,p,u,knotvectorU);
            dersbasisfuns(knotSpanU(k),p,order_deriv,u,knotvectorU, dersU[k]);
            knotSpanV(k) = FindSpan(nobv,q,v,knotvectorV);
            dersbasisfuns(knotSpanV(k),q,order_deriv,v,knotvectorV, dersV[k]);

            // extracting basis value and second derivatives from dersU&V
            for (int l = 0; l < p+1; ++l)
            {
                Nu(l) = dersU[k](0,l);
                Nv(l) = dersV[k](0,l);
                N1u(l) = dersU[k](1,l);
                N1v(l) = dersV[k](1,l);
                N2u(l) = dersU[k](2,l);
                N2v(l) = dersV[k](2,l);
            }

            // knot span vector -q/p -> for locating corresponding T
            for (int l = 0; l <= p; ++l)
            {
                ksU(k,l) = knotSpanU(k)+l-p;
            }
            // std::cout << std::endl;
            for (int l = 0; l <= q; ++l)
            {
                ksV(k,l) = knotSpanV(k)+l-q;
            }

            for (int l = 0; l < ksU.cols(); ++l)
            {
                for (int m = 0; m < ksV.cols(); ++m)
                {
                    ind = ksU(k,l)*lenv+ksV(k,m); // calculating corresponding position using knotspan
                    cm[0](k,ind) = Nu(l)*Nv(m); // assigning Ni*Nj value
                    cm[1](k,ind) = N1u(l)*Nv(m); // assigning N1i*Nj value
                    cm[2](k,ind) = Nu(l)*N1v(m); // assigning Ni*N1j value
                    cm[3](k,ind) = N1u(l)*N1v(m); // assigning N1i*N1j value
                    cm[4](k,ind) = N2u(l)*Nv(m); // assigning N2i*Nj value
                    cm[5](k,ind) = Nu(l)*N2v(m); // assigning Ni*N2j value
                    cm[6](k,ind) = N2u(l)*N2v(m); // assigning N2i*N2j value
                }
            }
        }         
    }
}

// Creating a matrix consisting of val
Eigen::MatrixXd NGfunc::matInit(int rows, int cols, float val)
{
    Eigen::MatrixXd out = Eigen::MatrixXd::Zero(rows, cols);
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            out(i,j) = val;
        }
    }
    return out;
}

// Creating a vector of zero matrices
std::vector<Eigen::MatrixXd> NGfunc::vectorMatInit(int n, int rows, int cols)
{
    std::vector<Eigen::MatrixXd> output;
    for (int i = 0; i < n; ++i) // constructing 7 vars: NuNv, N1uNv, NuN1v, N1uN1v, N2uNv, NuN2v, N2uN2v
    {
        output.push_back(Eigen::MatrixXd::Zero(rows, cols));
    }
    return output;
}

// Creating a vector of zero arraies
std::vector<Eigen::ArrayXXd> NGfunc::arrayMatInit(int n, int rows, int cols)
{
    std::vector<Eigen::ArrayXXd> output;
    for (int i = 0; i < n; ++i) 
    {
        std::cout << i << std::endl;
        output.push_back(Eigen::ArrayXXd::Zero(rows, cols));
    }
    return output;
}

// Solve for x from Ax = b
Eigen::MatrixXd NGfunc::linSol(Eigen::MatrixXd A, Eigen::MatrixXd b)
{
    // Eigen::MatrixXd x = A.colPivHouseholderQr().solve(b);
    // Eigen::MatrixXd x = A.householderQr().solve(b);
    Eigen::MatrixXd x = A.lu().solve(b);
    // Eigen::MatrixXd x = A.llt().solve(b);
    return x;
}

// calculate element-wise atan2
Eigen::ArrayXXd NGfunc::getAtan2(Eigen::ArrayXXd A, Eigen::ArrayXXd b, int lenu, int lenv)
{
    Eigen::ArrayXXd output = Eigen::ArrayXXd::Zero(lenu*lenv,1);
    for (int i = 0; i < lenu*lenv; i++)
    {
        output(i) = atan2(A(i), b(i));
    }
    return output;
}

// // calculate element-wise cos
// Eigen::MatrixXd NGfunc::getCos(Eigen::MatrixXd input, int lenu, int lenv)
// {
//     Eigen::MatrixXd output = Eigen::MatrixXd::Zero(lenu*lenv,1);
//     for (int i = 0; i < lenu*lenv; i++)
//     {
//         output(i) = cos(input(i));
//     }
//     return output;
// }

// // calculate element-wise sin
// Eigen::MatrixXd NGfunc::getSin(Eigen::MatrixXd input, int lenu, int lenv)
// {
//     Eigen::MatrixXd output = Eigen::MatrixXd::Zero(lenu*lenv,1);
//     for (int i = 0; i < lenu*lenv; i++)
//     {
//         output(i) = sin(input(i));
//     }
//     return output;
// }

// // calculate element-wise atan
// Eigen::MatrixXd NGfunc::getAtan(Eigen::MatrixXd input, int lenu, int lenv)
// {
//     Eigen::MatrixXd output = Eigen::MatrixXd::Zero(lenu*lenv,1);
//     for (int i = 0; i < lenu*lenv; i++)
//     {
//         output(i) = atan(input(i));
//     }
//     return output;
// }

// calculates epsilon and aap (a*a') based on phi, theta, NuN1v, and N1uNv.
void NGfunc::getEpsilonAndAap(std::vector<Eigen::MatrixXd>& output, float epsilonb, float delta, Eigen::MatrixXd phi, Eigen::MatrixXd xtheta, std::vector<Eigen::MatrixXd>& cm, int lenu, int lenv)
{
    // output: 0 - epsilon, 1 - epsilon_deriv, 2 - aap, 3 - P_dy, 4 - P_dx
    int aniso = 6;

    output[3] = (cm[1]*phi).reshaped(lenu*lenv,1); // P_dx = N1uNv*phi;
    output[4] = (cm[2]*phi).reshaped(lenu*lenv,1); // P_dy = NuN1v*phi;

    Eigen::MatrixXd check = (cm[1]).matrix();
    printArray2TXT(check, "./var_check/cm[1]_check.txt");
    printArray2TXT(phi, "./var_check/phi_in_ap_check.txt");
    check = (cm[1]*phi);
    printArray2TXT(check, "./var_check/cm[1]*phi_check.txt");
    check = output[3].matrix();
    printArray2TXT(check, "./var_check/output[3]_check.txt");
    
    Eigen::ArrayXXd atheta = getAtan2(output[4].array(), output[3].array(), lenu, lenv);

    output[0] = epsilonb*(1+delta*cos((aniso*(atheta-xtheta.array())))).matrix(); // epsilon
    output[1] = -epsilonb*(aniso*(delta*(sin((aniso*(atheta-xtheta.array())))))).matrix(); // epsilon_deriv
    output[2] = output[0].cwiseProduct(output[1]); // aap

    output[0] = output[0].reshaped(lenu*lenv,1);
    output[2] = output[2].reshaped(lenu*lenv,1);
}

// update r g based on nnT
void NGfunc::updateRgSg(Eigen::ArrayXXd& rMat, Eigen::ArrayXXd& sMat, Eigen::ArrayXXd nnT, int lenu, int lenv)
{
    for (int i = 0; i < lenu*lenv; i++)
    {
        if (nnT(i) == 1)
        {
            rMat(i) = 50;
            sMat(i) = 0;
        }
    }
}

Eigen::ArrayXXd NGfunc::regular_Heiviside_fun(Eigen::ArrayXXd input, int lenu, int lenv)
{
    float epsilon = 0.0001; // the number is not fixed.
    // H1E = 0.5*(1+(2/pi)*atan(X./epsilon));
    Eigen::ArrayXXd output = Eigen::ArrayXXd::Zero(lenu*lenv,1);
    for (int i = 0; i < lenu*lenv; i++)
    {
        output(i) =  0.5*(1+(2/M_PI)*atan(input(i)/epsilon));
    }
    return output;
}

void NGfunc::stiffMatSetupBCID(Eigen::MatrixXd& coll_Lhs, Eigen::MatrixXd& coll_Rhs, Eigen::MatrixXd bcid, Eigen::MatrixXd N, int lenu, int lenv)
{
    for (int i = 0; i < lenu*lenv; i++)
    {
        // change stiffness mat where there is non-zero boundary condition
        if (bcid(i)==1)
        {
            // coll_Lhs(i,0,lenu*lenv,1) = Eigen::MatrixXd::Zero(1,bcid.rows());
            coll_Lhs.row(i) = Eigen::MatrixXd::Zero(1,bcid.rows());
            coll_Lhs.col(i) = Eigen::MatrixXd::Zero(bcid.rows(),1);
            coll_Lhs(i,i) = 1;
            coll_Rhs(i) = N(i);
        }    }
}

Eigen::MatrixXd NGfunc::N1mulNN(Eigen::MatrixXd N1, Eigen::MatrixXd NN)
{
	Eigen::MatrixXd output = Eigen::MatrixXd::Zero(NN.rows(), NN.cols());
    for (int i = 0; i < NN.rows(); i++)
    {
        output.row(i) = N1(i) * NN.array().row(i);
    }
    return output;
}

Eigen::MatrixXd NGfunc::N1divNN(Eigen::MatrixXd N1, Eigen::MatrixXd NN)
{
	Eigen::MatrixXd output = Eigen::MatrixXd::Zero(NN.rows(), NN.cols());
    for (int i = 0; i < NN.rows(); i++)
    {
        output.row(i) = N1(i) / NN.array().row(i);
    }
    return output;
}

Eigen::MatrixXd NGfunc::sumFilter(Eigen::MatrixXd phi)
{
    int lenu,lenv;
    lenu = phi.rows(); lenv = phi.cols();
	Eigen::MatrixXd tips = Eigen::MatrixXd::Zero(lenu,lenv);

    return tips;
}
	
void NGfunc::Run(int n_process)
{
    // Start Simulation Model
	std::cout << "**********************************************************" << std::endl;
	std::cout << "2D Phase-field Neuron Growth solver using IGA-Collocation" << std::endl;

	// auto start = std::chrono::system_clock::now();
	// auto stop = std::chrono::system_clock::now();
	// auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
	// std::cout << duration.count() << " milliseconds" << std::endl;

	// Eigen::setNbThreads(n_process);

	// variable and png save frequency
	int var_save_invl = 100;
	int png_save_invl = 100;
	int png_plot_invl = 100;

	// variable setup
	int numNeuron = 1;
	int gc_sz = 2;
	float kappa = 2.0;
	int M_axon = 100;
	int M_neurites = 60;
	int delta = 0.20;

	// Simulation Parameter Initialization
	// time stepping variables
	float dtime = 1e-2;
	int end_iter = 1; // for debuging
	// int end_iter = 100000; // large enough number

	// setup domain size
	int Nx, Ny;
	neuron_domain_setup(numNeuron,Nx,Ny);

	// B-spline curve order (U,V direction)
	int p = 3;
	int q = 3;
	// vector<float> knotvectorU(Nx+1+6), knotvectorV(Nx+1+6);
	Eigen::VectorXd knotvectorU = Eigen::VectorXd::Zero(Nx+1+6);
	Eigen::VectorXd knotvectorV = Eigen::VectorXd::Zero(Nx+1+6);

	gen_knotvector(knotvectorU,0,Nx,Nx+1,q);
	gen_knotvector(knotvectorV,0,Nx,Nx+1,q);

	// setting lenu lenv this way for easier access to ghost nodes later on
	int lenu = knotvectorU.rows()-2*(p-1);
	int lenv = knotvectorV.rows()-2*(p-1);
	std::cout << "lenu: " << lenu << " | lenv: " << lenv << std::endl;

	// neuron growth variables
	int aniso = 6;
	// float kappa = 2;end // kappa= 2;
	float alpha = 0.9;
	float pix=4.0*atan(1.0);
	float alphOverPix = alpha/pix;
	int gamma = 10;
	float tau = 0.3;
	float M_phi = 60;
	float M_theta = 0.5*M_phi;
	float s_coeff = 0.007;
	// float delta = 0.1;
	float epsilonb = 0.04;

	// Tubulin parameters
	int r = 5;
	float g = 0.1;
	float alpha_t = 0.001;
	float beta_t = 0.001;
	int Diff = 4;
	float source_coeff = 15;

	// tolerance for NR iterations in phi equation
	float residual, tol;
	tol = 1e-4;

	// Initializing phi and concentration based on neuron seed position
	int seed_radius = 20; // Seed size
	Nx = 20; Ny = 20;
	seed_radius = 5;

	Eigen::VectorXd seed_x = Eigen::VectorXd::Zero(numNeuron);
	Eigen::VectorXd seed_y = Eigen::VectorXd::Zero(numNeuron);
	Eigen::MatrixXd phi = Eigen::MatrixXd::Zero(lenu, lenv);
	Eigen::MatrixXd conct = Eigen::MatrixXd::Zero(lenu, lenv);
	Eigen::MatrixXd theta = Eigen::MatrixXd::Random(lenu, lenv);
	Eigen::MatrixXd tempr = Eigen::MatrixXd::Zero(lenu, lenv);

	initialize_neurite_growth(phi, conct, seed_x, seed_y, seed_radius, lenu, lenv, numNeuron);

	printArray2TXT(phi, "./var_check/phi_check.txt");
	printArray2TXT(conct, "./var_check/conct_check.txt");
	printArray2TXT(theta, "./var_check/theta_check.txt");
	printArray2TXT(tempr, "./var_check/tempr_check.txt");

	// initializing initial phi,theta,tempr for boundary condition (Dirichlet)
	Eigen::MatrixXd phi_initial = phi;
	Eigen::MatrixXd conct_initial = conct;
	Eigen::MatrixXd theta_initial = theta;
	Eigen::MatrixXd tempr_initial = tempr;
	// std::cout << "test2.5" << std::endl;
	// for (int i = 2; i < lenu-2; i++)
	// {
	//     for (int j = 2; j < lenv-2; j++)
	//     {
	//         std::cout << i << "|" << j << std::endl;
	//         phi_initial(i,j) = 0;
	//         theta_initial(i,j) = 0;
	//         tempr_initial(i,j) = 0;
	//     }
	// }

	phi_initial = phi_initial.reshaped(lenu*lenv,1);
	conct_initial  = conct_initial.reshaped(lenu*lenv,1);
	theta_initial  = theta_initial.reshaped(lenu*lenv,1);
	tempr_initial  = tempr_initial.reshaped(lenu*lenv,1);

	phi = phi.reshaped(lenu*lenv,1);
	conct = phi.reshaped(lenu*lenv,1);
	tempr = tempr.reshaped(lenu*lenv,1);

	// Expanding domain parameters
	int BC_clearance = 25;
	int expd_sz = 10;

	// Iterating Variable Initialization
	// constructing collocation basis
	int order_deriv = 2;    // highest order of derivatives to calculate
	// std::vector<Eigen::MatrixXd> cm;
	// for (int i = 0; i < 7; ++i) // constructing 7 vars: NuNv, N1uNv, NuN1v, N1uN1v, N2uNv, NuN2v, N2uN2v
	// {
	//     cm.push_back(Eigen::MatrixXd::Zero(lenu*lenv, lenu*lenv));
	// }
	std::vector<Eigen::MatrixXd> cm = vectorMatInit(7, lenu*lenv, lenu*lenv);
	collocationDers(knotvectorU, p, knotvectorV, q, order_deriv, cm);

	Eigen::MatrixXd tmp = cm[0].block(0,0,80,80);
	// printArray2TXT(tmp, "./var_check/cm0_check.txt");
	tmp = cm[0].block(300,300,20,20);
	// printArray2TXT(tmp, "./var_check/cm0mid_check.txt");
	tmp = cm[0].block(3950,3950,19,19);
	// printArray2TXT(tmp, "./var_check/cm0end_check.txt");

	Eigen::MatrixXd lap = cm[5]+cm[6];

	std::cout << "phi shape: " << phi.rows() << "|" << phi.cols() << std::endl;
	phi = linSol(cm[0], phi);
	tmp = phi.reshaped(lenu,lenv);
	printArray2TXT(tmp, "./var_check/linSolPhi_check.txt");
}