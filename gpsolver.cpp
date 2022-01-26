//Source code for the FDM and Runge functions which discretise the 2-component coupled GP eqns 

#include "gpsolver.h"
#include <stdlib.h>
#include <iostream>
#include <stdio.h>

using namespace std;

double GPsolver::C(double k[]){
    
}

//spatially discretises system and returns -iMk 
double GPsolver::spatialDiscretiser(double k[]){

    double C_a, C_b; //constant introduced to write dimensionless discretised GP equation in more convenient way 
    C_a=GPsolver::C(k);
    C_b=GPsolver::C(k);

    //calculating -iMk for each gridpoint i
    for (int i=0; i<N; i++){
        if (i%2==0){ //even entries are for condensate a

        }else{ //odd entries for condensate b

        }
    }
}

void GPsolver::RK4(double psi[]){
    
    //declaring variables for RK4
    double k[N];
    double k_1[N], Mk_1[N];
    double k_2[N], Mk_2[N];
    double k_3[N], Mk_3[N];
    double k_4[N], Mk_4[N];

    //1st RK4 iteration
    for (int i=0; i<N; i++){
        k_1[i]=psi[i];
    }
    for (int i=0; i<N; i++){
        Mk_1[i]=GPsolver::spatialDiscretiser(k_1);
    }

    //2nd RK4 iteration
    for (int i=0; i<N; i++){
        k_2[i]=psi[i]+dt*Mk_1[i]/2;
    }
    for (int i=0; i<N; i++){
        Mk_2[i]=GPsolver::spatialDiscretiser(k_2);
    }
    
    //3rd RK4 iteration
    for (int i=0; i<N; i++){
        k_3[i]=psi[i]+dt*Mk_2[i]/2;;
    }
    for (int i=0; i<N; i++){
        Mk_3[i]=GPsolver::spatialDiscretiser(k_3);
    }
    
    //4th RK4 iteration
    for (int i=0; i<N; i++){
        k_4[i]=psi[i]+dt*Mk_3[i];;
    }
    for (int i=0; i<N; i++){
        Mk_4[i]=GPsolver::spatialDiscretiser(k_4);
    }
    
    //new wavefunction after dt time increment
    for (int i=0; i<N; i++){
        psi[i]=dt/6.0*(Mk_1[i]+2*Mk_2[i]+2*Mk_3[i]+Mk_4[i]);
    }    
}





