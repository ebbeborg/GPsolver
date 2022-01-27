//Source code for the FDM and Runge functions which discretise the 2-component coupled GP eqns 

#include "gpsolver.h"
#include <stdlib.h>
#include <cmath>

//generates initial psi
void Init_psi_generator(dcomp psi[]){
    //placeholder
}

//solves eigenproblem (resulting from discretisation) using RK4 method to get psi(a0,b0,a1,b1,...,aN-1,bN-1) at +dt
void GPsolver::RK4(dcomp psi[]){ //remember to multiply Mk's by -i
    
    //declaring variables for RK4
    dcomp k[N];
    dcomp k_1[N], Mk_1[N];
    dcomp k_2[N], Mk_2[N];
    dcomp k_3[N], Mk_3[N];
    dcomp k_4[N], Mk_4[N];

    //1st RK4 iteration
    for (int i=0; i<N; i++){
        k_1[i]=psi[i];
    }
    spatialDiscretiser(k_1, Mk_1); //calculating Mk_1

    //2nd RK4 iteration
    for (int i=0; i<N; i++){
        k_2[i]=psi[i]+dt*I*Mk_1[i]/2.;
    }
    spatialDiscretiser(k_2, Mk_2); //calculating Mk_2
    
    //3rd RK4 iteration
    for (int i=0; i<N; i++){
        k_3[i]=psi[i]+dt*I*Mk_2[i]/2.;
    }
    spatialDiscretiser(k_3, Mk_3); //calculating Mk_3
    
    //4th RK4 iteration
    for (int i=0; i<N; i++){
        k_4[i]=psi[i]+dt*I*Mk_3[i];
    }
    spatialDiscretiser(k_4, Mk_4); //calculating Mk_4
    
    //new wavefunction after dt time increment
    for (int i=0; i<N; i++){
        psi[i]=dt/6.*(Mk_1[i]+2.*Mk_2[i]+2.*Mk_3[i]+Mk_4[i]);
    }    
}

//spatially discretises RHS of coupled GP eqn in 1D using FDM and calculates Mk(a0,b0,a1,b1,...,aN-1,bN-1) 
void GPsolver::spatialDiscretiser(dcomp k[], dcomp Mk[]){

    double C[N]; //constant introduced for convenience 
    Const_calc(k, C); //calculates constant for each component at each gridpoint

    Mk[0]=C[0]*k[0]+omega/g*n_0*k[1]-k[2]-k[N-2]; //1st row of Mk (condensate a)
    Mk[1]=omega/g*n_0*k[0]+C[1]*k[1]-k[3]-k[N-1]; //2nd row of Mk (condensate b)
    Mk[N-2]=-k[0]-k[N-4]+C[N-2]*k[N-2]+omega/g*n_0*k[N-1]; //2nd last row of Mk (condensate a)
    Mk[N-1]=-k[1]-k[N-3]+omega/g*n_0*k[N-2]+C[N-1]*k[N-1]; //last row of Mk (condensate b)

    //calculating Mk for each gridpoint 
    for (int i=2; i<N-2; i++){
        if (i%2==0){ //even entries are for condensate a
            Mk[i]=-k[i-2]+C[i]*k[i]+omega/g*n_0*k[i+1]-k[i+2];
        }else{ //odd entries for condensate b
            Mk[i]=-k[i-2]+omega/g*n_0*k[i-1]+C[i]*k[i]-k[i+2]; //omega needs to be complex conj
        }
    }
}

//Calculates convenient constant for spatial discretisation C(a0,b0,a1,b1,...,aN-1,bN-1)
void GPsolver::Const_calc(dcomp k[], double C[]){
    for (int i=0; i<N; i++){
        if (i%2==0){ //even entries are for condensate a
            2/pow(dx,2)+V_a/(g*n_0)+pow(abs(k[i]),2)/n_0+g_ab*pow(abs(k[i+1]),2)/(g*n_0);
        }else{ //odd entries for condensate b
            2/pow(dx,2)+V_a/(g*n_0)+pow(abs(k[i]),2)/n_0+g_ab*pow(abs(k[i-1]),2)/(g*n_0);
        }
    }
}





