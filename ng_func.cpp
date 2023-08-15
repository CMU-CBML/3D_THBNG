#include <iostream>     /* cout */
#include <vector>
#include <math.h>       /* tanh */
#include <fstream>

#include <Eigen/Dense>

// using namespace std;

void neuron_domain_setup(int numNeuron, int& Nx, int& Ny)
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
}

// Creating knot vector (incrementing from lo to hi)
void gen_knotvector(Eigen::VectorXd& knotvector, int lo, int hi, int N, int order)
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
void printVector1D(Eigen::VectorXd& vec)
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

void printMatrix2D(Eigen::MatrixXd& vec)
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

void printArray2TXT(Eigen::MatrixXd& v, std::string fn)
{
    std::ofstream fout;
    
	fout.open(fn);

	for (int i = 0; i < v.rows(); i++)
    {
        for (int j = 0; j < v.cols(); j++)
        {
            // fout << v(i,j);
            if (v(i,j) == 0)
            {
                fout << "  ";
            } else{
                fout << "##";
            }
        }
        fout << std::endl;
    }

	fout.close();
}

void initialize_neurite_growth(Eigen::MatrixXd& phi, Eigen::MatrixXd& conct, Eigen::VectorXd& seed_x, Eigen::VectorXd& seed_y, int seed_radius, int lenu, int lenv, int numNeuron)
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

    printArray2TXT(phi, "./phi_initial.txt");
    printArray2TXT(conct, "./conct_initial.txt");
}

// Find knot span
int FindSpan(int n, int p, float u, Eigen::VectorXd U)
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
void dersbasisfuns(int i, int p, int order_deriv, float u, Eigen::VectorXd U, Eigen::MatrixXd& ders)
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

void collocationDers(Eigen::VectorXd knotvectorU,int p, Eigen::VectorXd knotvectorV, int q, int order_deriv, std::vector<Eigen::MatrixXd>& cm)
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
            coordx = coordx + knotvectorU(i+l);
        } 
        coordx = coordx/p;

        for (int j = 0; j < knotvectorV.rows()-q-1; ++j)
        {
            coordy = 0;
            for (int l = 0; l < q; ++l)
            {
                coordy = coordy + knotvectorV(j+l);
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
Eigen::MatrixXd matInit(int rows, int cols, float val)
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
std::vector<Eigen::MatrixXd> vectorMatInit(int n, int rows, int cols)
{
    std::vector<Eigen::MatrixXd> output;
    for (int i = 0; i < n; ++i) // constructing 7 vars: NuNv, N1uNv, NuN1v, N1uN1v, N2uNv, NuN2v, N2uN2v
    {
        output.push_back(Eigen::MatrixXd::Zero(rows, cols));
    }
    return output;
}

// Creating a vector of zero arraies
std::vector<Eigen::ArrayXd> arrayMatInit(int n, int rows, int cols)
{
    std::vector<Eigen::ArrayXd> output;
    for (int i = 0; i < n; ++i) // constructing 7 vars: NuNv, N1uNv, NuN1v, N1uN1v, N2uNv, NuN2v, N2uN2v
    {
        output.push_back(Eigen::ArrayXd::Zero(rows, cols));
    }
    return output;
}

// Solve for x from Ax = b
Eigen::MatrixXd linSol(Eigen::MatrixXd A, Eigen::MatrixXd b)
{
    // Eigen::MatrixXd x = A.colPivHouseholderQr().solve(b);
    // Eigen::MatrixXd x = A.lu().solve(b);
    Eigen::MatrixXd x = A.llt().solve(b);
    return x;
}

// calculate element-wise atan2
Eigen::MatrixXd getAtan2(Eigen::MatrixXd A, Eigen::MatrixXd b, int lenu, int lenv)
{
    Eigen::MatrixXd output = Eigen::MatrixXd::Zero(lenu*lenv,1);
    for (int i = 0; i < lenu*lenv; i++)
    {
        output(i) = atan2(A(i), b(i));
    }
    return output;
}

// calculate element-wise cos
Eigen::MatrixXd getCos(Eigen::MatrixXd input, int lenu, int lenv)
{
    Eigen::MatrixXd output = Eigen::MatrixXd::Zero(lenu*lenv,1);
    for (int i = 0; i < lenu*lenv; i++)
    {
        output(i) = cos(input(i));
    }
    return output;
}

// calculate element-wise sin
Eigen::MatrixXd getSin(Eigen::MatrixXd input, int lenu, int lenv)
{
    Eigen::MatrixXd output = Eigen::MatrixXd::Zero(lenu*lenv,1);
    for (int i = 0; i < lenu*lenv; i++)
    {
        output(i) = sin(input(i));
    }
    return output;
}

// calculate element-wise atan
Eigen::MatrixXd getAtan(Eigen::MatrixXd input, int lenu, int lenv)
{
    Eigen::MatrixXd output = Eigen::MatrixXd::Zero(lenu*lenv,1);
    for (int i = 0; i < lenu*lenv; i++)
    {
        output(i) = atan(input(i));
    }
    return output;
}

// calculates epsilon and aap (a*a') based on phi, theta, NuN1v, and N1uNv.
void getEpsilonAndAap(std::vector<Eigen::MatrixXd>& output, float epsilonb, float delta, Eigen::MatrixXd phi, Eigen::MatrixXd xtheta, std::vector<Eigen::MatrixXd>& cm, int lenu, int lenv)
{
    // output: 0 - epsilon, 1 - epsilon_deriv, 2 - aap, 3 - P_dy, 4 - P_dx

    int aniso = 6;

    output[3] = (cm[1]*phi).reshaped(lenu*lenv,1); // P_dx = N1uNv*phi;
    output[4] = (cm[2]*phi).reshaped(lenu*lenv,1); // P_dy = NuN1v*phi;

    Eigen::MatrixXd atheta = getAtan2(output[4], output[3], lenu, lenv);
    Eigen::MatrixXd epbMat = matInit(lenu*lenv,1,epsilonb);
    Eigen::MatrixXd anisoMat = matInit(lenu*lenv,1,aniso);
    Eigen::MatrixXd deltaMat = matInit(lenu*lenv,1,delta);

    output[0] = epbMat.cwiseProduct(matInit(lenu*lenv,1,1)+deltaMat.cwiseProduct(getCos((anisoMat.cwiseProduct(atheta-xtheta)),lenu,lenv))); // epsilon
    output[1] = -epbMat.cwiseProduct(anisoMat.cwiseProduct(deltaMat.cwiseProduct(getSin((anisoMat.cwiseProduct(atheta-xtheta)),lenu,lenv)))); // epsilon_deriv
    output[2] = output[0].cwiseProduct(output[1]); // aap

    output[0] = output[0].reshaped(lenu*lenv,1);
    output[2] = output[2].reshaped(lenu*lenv,1);
}

// update r g based on nnT
void updateRgSg(Eigen::ArrayXd& rMat, Eigen::ArrayXd& sMat, Eigen::ArrayXd nnT, int lenu, int lenv)
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

Eigen::ArrayXd regular_Heiviside_fun(Eigen::ArrayXd input, int lenu, int lenv)
{
    float epsilon = 0.0001; // the number is not fixed.
    // H1E = 0.5*(1+(2/pi)*atan(X./epsilon));
    Eigen::ArrayXd output = Eigen::ArrayXd::Zero(lenu*lenv,1);
    for (int i = 0; i < lenu*lenv; i++)
    {
        output(i) =  0.5*(1+(2/M_PI)*atan(input(i)/epsilon));
    }
    return output;
}

void stiffMatSetupBCID(Eigen::MatrixXd& coll_Lhs, Eigen::MatrixXd& coll_Rhs, Eigen::MatrixXd bcid, Eigen::MatrixXd N, int lenu, int lenv)
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
        }
    }
}