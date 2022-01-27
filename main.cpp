//ihd/dt(x_a)=((-1/2m)(hbar^2)(Del^2) + g_a(x_a)^2 +g_ab(x_b)^2)x_a + omega(x_b)
//ihd/dt(x_b)=((-1/2m)(hbar^2)(Del^2) + g_b(x_b)^2 +g_ab(x_a)^2)x_b + omega^*(x_a)
//generating psi_init?
//healing factor?

#include <iostream>
#include <complex>
#include "gpsolver.h"

int main(){

//declaring objects to store BEC parameters and functions
    BEC_parameters parameter;
    GPsolver GPsolver; 

    //input discretisation parameters
    std::cout<<"Input gridsize:"; 
    std::cin>>parameter.gridsize;
        if (parameter.gridsize%2!=0){
            std::cout<<"Please enter an even grid size N"<<std::endl;
            exit(1);                    
        }
    std::cout<<"Input grid spacing:";
    std::cin>>parameter.dx;
    std::cout<<"Input runtime:";
    std::cin>>parameter.runtime;
    //std::cout<<"Input timestep:"; set to 1 atm
    //std::cin>>parameter.dt;

    //input BEC parameters
    std::cout<<"Input initial density of condensate n_0:"; 
    std::cin>>parameter.n_0;
    std::cout<<"Input intraspecies (a on a, b on b) contact interaction g:"; 
    std::cin>>parameter.g;
    std::cout<<"Input interspecies (a on b, b on a) contact interaction g_ab:"; 
    std::cin>>parameter.g_ab;
    std::cout<<"Input coherent coupling:"; 
    std::cin>>parameter.omega;

    //generating initial wavefunction of condensate
    dcomp psi_init[parameter.N], psi[parameter.N];

    for (int i=0; i<parameter.N; i++){
        psi[i]=psi_init[i];
    }

    //evaluating psi in time increments using RK4 using initial condensate wavefunction psi_init
    for(int t=0; t<parameter.runtime; t++ ){
        GPsolver.RK4(psi); //iterates psi by one time step
        //pipes iterated psi to results file
    }

    return 0;
}