//Source code for the FDM and Runge functions which discretise the 2-component coupled GP eqns 

#include "gpsolver.h"
#include <stdlib.h>
#include <iostream>
#include <stdio.h>

using namespace std;

double GPsolver::C(double k[], BEC_parameters& parameters){

}

double GPsolver::spatialDiscretiser(double k[], BEC_parameters& parameters){

    double C_a, C_b; //constant introduced to write dimensionless discretised GP equation in more convenient way 
    C_a=gpsolver::C(k, parameters);
    C_b=gpsolver::C(k, parameters);

    //calculating Mk for each gridpoint i
    for (int i=0; i<N; i++){
        if (i%2==0){ //even entries are for condensate a

        }else{ //odd entries for condensate b

        }
    }
}

void GPsolver::RK4(double psi_init[], BEC_parameters& parameters){
    
    double k[parameters.N];
    double k_1[parameters.N];//, Mk_1[];
    double k_2[parameters.N];//, Mk_2[];
    double k_3[parameters.N];//, Mk_3[];

    //1st RK4 iteration
    for (int i=0; i<parameters.N; i++){
        k_1[i]=psi_init[i];
    }
    Mk_1=GPsolver::spatialDiscretiser(k_1, parameters);
}





