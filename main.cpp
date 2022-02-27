//GPE is dimensionless and removed e^-imut/hbar time dependance
#include <iostream>
#include <fstream>
#include "gpsolver.h"

int main(){

    //declaring objects
    BEC_parameters parameter; //declaring object to store BEC parameters
    GPsolver GPsolver; //declaring object that allows access to GPsolver functions

    //declaring parameters
    dcomp psi[parameter.N]; //condensate wavefunction storing components a & b
    double x[parameter.gridsize]= {}; //initialising spatial grid of zeros
    double omega[parameter.N]; //coherent coupling    

    //generating initial condensate system at t=0, ie wavefunction psi, spatial grid x, coherent coupling omega
    GPsolver.System_generator(x);
    GPsolver.Init_psi_generator(psi, true, x); //(orderparameter, excitation true/false, gridspace)

    //evaluating psi in time increments using RK4
    int a=1/parameter.dt; //time step normalisation factor
    for(int t=1; t<a*parameter.runtime; t++ ){
        
        GPsolver.RK4(psi); //iterates psi by one time step dt
        
        if(t%a==0){ //saving every 100th iteration
            
            std::ofstream output; //opening up results file
            output.open("results/results.txt", std::ios_base::app);
            
            for(int i=0; i<parameter.N; i++){
                output<<norm(psi[i]); //saving results to "results.txt"
                if(i<parameter.N-1){
                    output<<",";
                }
            }

            output<<"\r\n"; //new line so that next iteration of psi can be appended correctly
            output.close();
        }
    }
    return 0;
}