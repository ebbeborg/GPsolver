//Source code for the FDM and Runge functions which discretise the 2-component coupled GP eqns 

#include "gpsolver.h"
#include <stdlib.h>
#include <iostream>
#include <stdio.h>

using namespace std;

double gpsolver::discretisedConstant(double g, double g_ab, double n_0, double dx){

}

void gpsolver::spatialDiscretiser(int N, matrix & M, double g, double n_0, double G_a, double G_b, double omega){
    
    //check if gridsize is even
    if (N%2!=0){
        cout<<"Please enter an even grid size N"<<endl;
        exit(1);                    
    }

    //filling out top and bottom 2 rows of matrix M to ensure 1D grip loops
    M[0][N-2]=M[0][N-1]=M[N-2][0]=M[N-1][1]=-1;
    M[0][0]=M[N-2][N-2]=G_a;
    M[0][1]=M[N-2][N-1]=1;//omega/gn_0;
    M[1][1]=M[N-1][N-1]=G_b;
    M[1][0]=M[N-1][N-2]=1;//omegaconj/gn_0;

    //filling out all other rows of M
    for (int i=2; i<N-1; i++){
        if (i%2==0){ //for GPa
            M[i][i-2]=-1;
            M[i][i]=G_a;
            M[i][i+1]=omega/(g*n_0);
            M[i][i+2]=-1;
        }else{ //for GPb
            M[i][i-2]=-1;
            M[i][i-1]=omega/(g*n_0);
            M[i][i]=G_b;
            M[i][i+2]=-1;
        }
    }


}

void gpsolver::temporalDiscretiser(matrix & M){

}





