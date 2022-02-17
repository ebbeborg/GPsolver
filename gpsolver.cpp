//Source code for the FDM and Runge functions which discretise the 2-component coupled GP eqns 

#include "gpsolver.h"
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <fstream>

dcomp I=dcomp(0.,1.); //defining i

//generates initial psi
void GPsolver::Init_psi_generator(dcomp psi[], bool excitation, double x[]){

    //creating spatial grid
    for (int i=1; i<gridsize; i++){
        x[i]=x[i-1]+dx;
    }

    //opening up results file
    std::ofstream output;
    output.open("results/results.txt");

    //generating initial psi for each gridpoint and saving result
    for (int i=0; i<N; i++){
        psi[i]=sqrt(n_0/2); //since norm(psi)=n
        
        if(excitation){ //add excitation at x_0
            if(i%2==0){ //condensate a
                psi[i]+=0.001*exp(I*(k_0*x[i/2])-pow((x[i/2]-x_0)/width,2)/2); 
            }else{ //condensate b
                psi[i]-=0.001*exp(I*(k_0*x[(i-1)/2])-pow((x[(i-1)/2]-x_0)/width,2)/2);
            }
        }

        output<<norm(psi[i]); 
        if(i<N-1){
            output<<",";
        }
    }

    output<<"\r\n";
    output.close();
}

//solves eigenproblem (resulting from discretisation) using RK4 method to get psi(a0,b0,a1,b1,...,aN-1,bN-1) at +dt
void GPsolver::RK4(dcomp psi[]){ //remember to multiply Mk's by I=-i

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
        k_2[i]=psi[i]+dt*Mk_1[i]/2.;
    }
    spatialDiscretiser(k_2, Mk_2); //calculating Mk_2

    //3rd RK4 iteration
    for (int i=0; i<N; i++){
        k_3[i]=psi[i]+dt*Mk_2[i]/2.;
    }
    spatialDiscretiser(k_3, Mk_3); //calculating Mk_3
    
    //4th RK4 iteration
    for (int i=0; i<N; i++){
        k_4[i]=psi[i]+dt*Mk_3[i];
    }
    spatialDiscretiser(k_4, Mk_4); //calculating Mk_4

    //calculating new psi after dt time increment
    for (int i=0; i<N; i++){
        psi[i]+=dt/6.*(Mk_1[i]+2.*Mk_2[i]+2.*Mk_3[i]+Mk_4[i]);
    }
}

//spatially discretises RHS of coupled GP eqn in 1D using FDM and calculates -iMk(a0,b0,a1,b1,...,aN-1,bN-1) 
void GPsolver::spatialDiscretiser(dcomp k[], dcomp Mk[]){

    dcomp C[N]; //constant introduced for convenience 
    Const_calc(k, C); //calculates constant for each component at each gridpoint

    //calculating -iMk for each gridpoint, (N+i)%N to make grid loop 
    for (int i=0; i<N; i++){
        if (i%2==0){ //even entries are for condensate a
            //Mk[i]=-I*(-k[(N+i-2)%N]/pow(dx,2)+C[i]*k[i]+omega/g*n_0*k[i+1]-k[(N+i+2)%N]/pow(dx,2));
            Mk[i]=-I*(-k[(N+i-2)%N]/pow(dx,2)+C[i]*k[i]+omega*k[i+1]-k[(N+i+2)%N]/pow(dx,2));
        }else{ //odd entries for condensate b
            //Mk[i]=-I*(-k[(N+i-2)%N]/pow(dx,2)+omega/g*n_0*k[i-1]+C[i]*k[i]-k[(N+i+2)%N]/pow(dx,2)); //omega needs to be complex conj
            Mk[i]=-I*(-k[(N+i-2)%N]/pow(dx,2)+omega*k[i-1]+C[i]*k[i]-k[(N+i+2)%N]/pow(dx,2));
        }
    }
}

//Calculates convenient constant for spatial discretisation C(a0,b0,a1,b1,...,aN-1,bN-1)
void GPsolver::Const_calc(dcomp k[], dcomp C[]){
    for (int i=0; i<N; i++){
        if (i%2==0){ //even entries are for condensate a
            //C[i]=2/pow(dx,2)+V_a/(g*n_0)+norm(k[i])/n_0+g_ab*norm(k[i+1])/(g*n_0);
            //C[i]=2/pow(dx,2)+g_ab*(norm(k[i+1])-norm(k[i]))/(g*n_0)+abs(omega)/(g*n_0);
            C[i]=2/pow(dx,2)+g_ab*(norm(k[i+1])-norm(k[i]))/(g*n_0)+abs(omega); 
        }else{ //odd entries for condensate b
            //C[i]=2/pow(dx,2)+V_a/(g*n_0)+norm(k[i])/n_0+g_ab*norm(k[i-1])/(g*n_0);
            //C[i]=2/pow(dx,2)+g_ab*(norm(k[i-1])-norm(k[i]))/(g*n_0)+abs(omega)/(g*n_0);
            C[i]=2/pow(dx,2)+g_ab*(norm(k[i-1])-norm(k[i]))/(g*n_0)+abs(omega);
        }
    }
}





