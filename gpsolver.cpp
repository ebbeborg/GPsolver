//Source code for the FDM and Runge functions which discretise the 2-component coupled GP eqns 

#include "gpsolver.h"
#include <stdlib.h>
#include <iostream>
#include <stdio.h>

using namespace std;

double gpsolver::C(double psi[], double g, double g_ab, double n_0, double dx){

}

void gpsolver::spatialDiscretiser(int N, double psi[], double g, double g_ab, double n_0, double dx , double omega){

    double C_a, C_b; //constant introduced to write dimensionless discretised GP equation in more convenient way
    C_a=gpsolver::C(psi, g, g_ab, n_0, dx);
    C_b=gpsolver::C(psi, g, g_ab, n_0, dx);

    //calculating Mx for each gridpoint i
    for (int i=0; i<N; i++){
        if (i%2==0){ //even entries are for condensate a

        }else{ //odd entries for condensate b

        }
    }
}

void gpsolver::RK4(int N, double psi[], double g, double g_ab, double n_0, double dx, double omega){

    gpsolver::spatialDiscretiser(N, psi, g, g_ab, n_0, dx, omega);
}





