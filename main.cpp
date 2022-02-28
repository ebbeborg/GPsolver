//Spatially discretises the RHS of coupled GP eqns (for 2-component homogenous BEC) in 1D using finite differential method
//Returns RHS of coupled GP eqns in the form of a matrix M and coloumn vector x (the discretised order parameter) 
#include "gpsolver.h"
#include <fstream>

int main(){

    //declaring objects
    BEC_parameters parameter; //stores BEC parameters
    GPsolver GPsolver; //allows access to GPsolver functions

    //declaring parameters
    double x[parameter.gridsize]={}; //initialising spatial grid of zeros
    dcomp psi[parameter.N]; //condensate wavefunction psi(a0, b0, a1, b1, etc.)
    double omega[parameter.N]; //coherent coupling

    //generating gridspace
    GPsolver.Gridspace(x);
    
    //generating initial condensate wavefunction of system (true/false to add/remove excitation)
    GPsolver.Init_psi_generator(psi, true, x); 

    //generating spatial distribution of coherent coupling that will modulate system
    GPsolver.Modulator(omega);

    //evaluating psi in time increments using RK4
    int a=1/parameter.dt; //time step normalisation factor
    for(int t=1; t<a*parameter.runtime; t++ ){
        
        GPsolver.RK4(psi, omega); //iterates psi by one time step dt
        
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