//ihd/dt(x_a)=((-1/2m)(hbar^2)(Del^2) + g_a(x_a)^2 +g_ab(x_b)^2)x_a + omega(x_b)
//ihd/dt(x_b)=((-1/2m)(hbar^2)(Del^2) + g_b(x_b)^2 +g_ab(x_a)^2)x_b + omega^*(x_a)
//step time dt?
//generating psi_init?

#include "gpsolver.h"
#include <iostream>

int main(){

//declaring objects to store BEC parameters and functions
BEC_parameters parameters;
GPsolver homoGP;

    //input discretisation parameters
    std::cout<<"Input gridsize:"; 
    std::cin>>parameters.gridsize;
        if (parameters.gridsize%2!=0){
            std::cout<<"Please enter an even grid size N"<<std::endl;
            exit(1);                    
        }
    std::cout<<"Input grid spacing:";
    std::cin>>parameters.dx;

    //input BEC parameters
    std::cout<<"Input initial density of condensate n_0:"; 
    std::cin>>parameters.n_0;
    std::cout<<"Input intraspecies (a on a, b on b) contact interaction g:"; 
    std::cin>>parameters.g;
    std::cout<<"Input interspecies (a on b, b on a) contact interaction g_ab:"; 
    std::cin>>parameters.g_ab;
    std::cout<<"Input coherent coupling:"; 
    std::cin>>parameters.omega;

    //initial size of wavefunction (has to accomodate both a and b so 2*N)
    parameters.N=2*parameters.gridsize;
    //initial condensate wavefunction
    double psi_init[parameters.N];

    //evaluating psi in time increments using RK4 using initial condensate wavefunction psi_init
    homoGP.RK4(psi_init, parameters);

    //piping results to file

    return 0;
}