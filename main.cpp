//to do:
//generating psi_init?
//piping psi results to file

#include <iostream>
#include <complex>
#include "gpsolver.h"

int main(){

    //declaring class objects
    BEC_parameters parameter; //declaring object to store BEC parameters
    GPsolver GPsolver; //declaring object that allows access to GPsolver functions

    //input discretisation parameters and storing in parameter class
    std::cout<<"Input gridsize:"; 
    std::cin>>parameter.gridsize;
        if (parameter.gridsize%2!=0){
            std::cout<<"Please enter an even grid size N"<<std::endl;
            exit(1);                    
        }

    std::cout<<"Input grid spacing:";
    std::cin>>parameter.dx;
    //std::cout<<"Input runtime:";
    //std::cin>>parameter.runtime;

    //input BEC parameters and storing in parameter class
    //std::cout<<"Input intraspecies (a on a, b on b) contact interaction g:"; 
    //std::cin>>parameter.g;
    //std::cout<<"Input interspecies (a on b, b on a) contact interaction g_ab:"; 
    //std::cin>>parameter.g_ab;
    //std::cout<<"Input coherent coupling:"; 
    //std::cin>>parameter.omega;

    parameter.N=2*parameter.gridsize;
    parameter.n_0=1/parameter.dx;
    parameter.V_a=0;
    parameter.V_b=0;
    parameter.omega=0;

    //generating initial condensate wavefunction psi at t=0 and saving to results file
    dcomp psi[parameter.N];
    GPsolver.Init_psi_generator(psi);
    std::cout<<parameter.n_0<<std::endl;

    //evaluating psi in time increments using RK4
    //for(int t=0; t<parameter.runtime; t++ ){
    //    GPsolver.RK4(psi); //iterates psi by one time step dt
    //}

    return 0;
}