//ihd/dt(x_a)=((-1/2m)(hbar^2)(Del^2) + g_a(x_a)^2 +g_ab(x_b)^2)x_a + omega(x_b)
//ihd/dt(x_b)=((-1/2m)(hbar^2)(Del^2) + g_b(x_b)^2 +g_ab(x_a)^2)x_b + omega^*(x_a)

#include "gpsolver.h"
#include <iostream>

int main()
{
    int N; //gridsize
    double dx; //grid spacing
    double n_0; //initial density of 2D condensate
    double V; //external potential
    double g, g_ab;
    double G_a, G_b; //constant introduced to write dimensionless discretised GP equation in more convenient way
    double coherentCoupling;
    gpsolver homoGP;

    //input discretisation parameters
    std::cout<<"Input gridsize:"; 
    std::cin>>N;
    std::cout<<"Input grid spacing:"; 
    std::cin>>dx;

    //defining a matrix variable type
    matrix M;

    //input BEC parameters
    std::cout<<"Input initial density of condensate n_0:"; 
    std::cin>>n_0;
    std::cout<<"Input intraspecies (a on a, b on b) contact interaction g:"; 
    std::cin>>g;
    std::cout<<"Input interspecies (a on b, b on a) contact interaction g_ab:"; 
    std::cin>>g_ab;
    std::cout<<"Input coherent coupling:"; 
    std::cin>>coherentCoupling;

    //dimensionalising constant G
    G_a=homoGP.discretisedConstant(g, g_ab, n_0, dx);
    G_b=homoGP.discretisedConstant(g, g_ab, n_0, dx);

    homoGP.spatialDiscretiser(N, M, G_a, G_b, g, n_0, coherentCoupling);
    homoGP.temporalDiscretiser(M);

    return 0;
}