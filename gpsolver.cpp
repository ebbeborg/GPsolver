//Source code for the FDM and Runge functions which discretise the 2-component coupled GP eqns 

#include "gpsolver.h"
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <fstream>

extern const dcomp I=dcomp(0.,-1.); //defining I=-i

//generates initial psi
void GPsolver::Init_psi_generator(dcomp psi[]){
    
    std::cout<<n_0<<std::endl;
    //opening up results file
    std::ofstream output;
    output.open("results.txt");

    //generating initial psi for each gridpoint and saving result
    for (int i=0; i<N; i++){
        psi[i]=sqrt(n_0); //since norm(psi)=n
        output<<psi[i]<<","; 
    }
    output<<"/n";
    output.close();
}

//solves eigenproblem (resulting from discretisation) using RK4 method to get psi(a0,b0,a1,b1,...,aN-1,bN-1) at +dt
void GPsolver::RK4(dcomp psi[]){ //remember to multiply Mk's by I=-i

    //declaring variables for RK4
    double dt; //timestep set to 1 right now
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
    
    //opening up results file
    std::ofstream output;
    output.open("results.txt", std::ios_base::app);

    //calculating new psi after dt time increment
    for (int i=0; i<N; i++){
        psi[i]=dt/6.*(Mk_1[i]+2.*Mk_2[i]+2.*Mk_3[i]+Mk_4[i]);
        output << psi[i]<<","; //saving results to file
    }
    output<<"\n"; //new line so that next iteration of psi can be appended correctly
    output.close();
}

//spatially discretises RHS of coupled GP eqn in 1D using FDM and calculates Mk(a0,b0,a1,b1,...,aN-1,bN-1) 
void GPsolver::spatialDiscretiser(dcomp k[], dcomp Mk[]){

    double C[N]; //constant introduced for convenience 
    Const_calc(k, C); //calculates constant for each component at each gridpoint

    //calculating Mk for each gridpoint, (N+i)%N to make grid loop 
    for (int i=0; i<N; i++){
        if (i%2==0){ //even entries are for condensate a
            Mk[i]=-k[(N+i-2)%N]+C[i]*k[i]+omega/g*n_0*k[i+1]-k[(N+i+2)%N];
        }else{ //odd entries for condensate b
            Mk[i]=-k[(N+i-2)%N]+omega/g*n_0*k[i-1]+C[i]*k[i]-k[(N+i+2)%N]; //omega needs to be complex conj
        }
    }
}

//Calculates convenient constant for spatial discretisation C(a0,b0,a1,b1,...,aN-1,bN-1)
void GPsolver::Const_calc(dcomp k[], double C[]){
    for (int i=0; i<N; i++){
        if (i%2==0){ //even entries are for condensate a
            2/pow(dx,2)+V_a/(g*n_0)+norm(k[i])/n_0+g_ab*norm(k[i+1])/(g*n_0);
        }else{ //odd entries for condensate b
            2/pow(dx,2)+V_a/(g*n_0)+norm(k[i])/n_0+g_ab*norm(k[i-1])/(g*n_0);
        }
    }
}





