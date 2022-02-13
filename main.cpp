//GPE is dimensionless

#include <iostream>
#include <fstream>
#include "gpsolver.h"

int main(){

    //declaring class objects
    BEC_parameters parameter; //declaring object to store BEC parameters
    GPsolver GPsolver; //declaring object that allows access to GPsolver functions

    //generating spatial grid
    double x[parameter.gridsize]={};

    //generating initial condensate wavefunction psi at t=0 and saving to results file
    //pass true if generate excitation on top of GS
    dcomp psi[parameter.N];
    GPsolver.Init_psi_generator(psi, true, x); 

    //evaluating psi in time increments using RK4
    for(int t=1; t<100*parameter.runtime; t++ ){
        
        GPsolver.RK4(psi); //iterates psi by one time step dt
        
        if(t%10==0){ //saving every 10th iteration
            
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