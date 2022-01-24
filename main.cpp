//ihd/dt(x_a)=((-1/2m)(hbar^2)(Del^2) + g_a(x_a)^2 +g_ab(x_b)^2)x_a + omega(x_b)
//ihd/dt(x_b)=((-1/2m)(hbar^2)(Del^2) + g_b(x_b)^2 +g_ab(x_a)^2)x_b + omega^*(x_a)

#include "gpsolver.h"
#include <iostream>

int main()
{
    //declaring parameters
    int N; //gridsize
    double dx; //grid spacing
    double n_0; //initial density of 2D condensate
    double V; //external potential
    double g, g_ab;
    double coherentCoupling;
    gpsolver homoGP;

    //input discretisation parameters
    std::cout<<"Input gridsize:"; 
    std::cin>>N;
        if (N%2!=0){
            std::cout<<"Please enter an even grid size N"<<std::endl;
            exit(1);                    
        }
    std::cout<<"Input grid spacing:";
    std::cin>>dx;

    //input BEC parameters
    std::cout<<"Input initial density of condensate n_0:"; 
    std::cin>>n_0;
    std::cout<<"Input intraspecies (a on a, b on b) contact interaction g:"; 
    std::cin>>g;
    std::cout<<"Input interspecies (a on b, b on a) contact interaction g_ab:"; 
    std::cin>>g_ab;
    std::cout<<"Input coherent coupling:"; 
    std::cin>>coherentCoupling;

    //initial size of wavefunction (has to accomodate both a and b so 2*N)
    int temp=2*N;
    //initial condensate wavefunction
    double psi[temp];

    //evaluating psi in time increments using RK4
    homoGP.RK4(N, psi, g, g_ab, n_0, dx, coherentCoupling);

    //piping results to file

    return 0;
}