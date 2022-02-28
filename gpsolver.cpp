//Source code for the GPsolver class functions which discretise and solve the 2-component coupled GP eqns in time
#include "gpsolver.h"
#include <fstream>

dcomp I=dcomp(0.,1.); //defining complex i=sqrt(-1)

void GPsolver::Gridspace(double x[]){
    //generating spatial grid 
    for (int i=1; i<gridsize; i++){
        x[i]=x[i-1]+dx;
    }
}

//generating coherent coupling for modulation of system
void GPsolver::Modulator(double omega[]){
    for (int i=0; i<N; i++){
        if(i<N/2){
            omega[i]=omegaLHS;
        }else{
            omega[i]=omegaRHS;
        }
    }
}

//calculates Bogoliubov positive (u) and negative (v) mode amplitudes for given k_0
void GPsolver::Bogoliubov_mode_amplitudes(double &u, double &v){
    u=1.00007;
    v=-0.01219;
}

//generates initial psi and excitation
void GPsolver::Init_psi_generator(dcomp psi[], bool excitation, double x[]){

    //calculating bogoliubov excitation mode amplitudes u and v using initial wavevector k_0
    double u, v;
    Bogoliubov_mode_amplitudes(u, v);

    //opening up results file
    std::ofstream output;
    output.open("results/results.txt");

    //generating initial psi for each gridpoint and saving result
    for (int i=0; i<N; i++){
        psi[i]=sqrt(n_0/2); //since norm(psi)=n/2
        
        if(excitation){ //add excitation at x_0
            if(i%2==0){ //condensate a
                psi[i]=psi[i]*(1.+A*exp(-pow((x[i/2]-x_0)/width,2)/2)*(u*exp(I*k_0*x[i/2])+v*exp(-I*k_0*x[i/2]))); 
            }else{ //condensate b
                psi[i]=psi[i]*(1.-A*exp(-pow((x[(i-1)/2]-x_0)/width,2)/2)*(u*exp(I*k_0*x[(i-1)/2])+v*exp(-I*k_0*x[(i-1)/2]))); 
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
void GPsolver::RK4(dcomp psi[], double omega[]){ 

    //declaring variables for RK4
    dcomp psi_temp[N];
    //slopes at various time increments of psi
    dcomp k_1[N];
    dcomp k_2[N];
    dcomp k_3[N];
    dcomp k_4[N];

    //1st iteration, calculating slope k_1 at initial psi(t0)
    Spatial_discretiser(psi, k_1, omega);
    
    //2nd iteration, calculating slope k_2 at 1st psi_temp increment
    for (int i=0; i<N; i++){
        psi_temp[i]=psi[i]+dt*k_1[i]/2.;
    }
    Spatial_discretiser(psi_temp, k_2, omega);
    
    //3rd iteration, calculating slope k_3 at 2nd psi_temp increment
    for (int i=0; i<N; i++){
        psi_temp[i]=psi[i]+dt*k_2[i]/2.;
    }
    Spatial_discretiser(psi_temp, k_3, omega);
    
    //4th iteration, calculating slope k_4 at 3rd psi_temp increment
    for (int i=0; i<N; i++){
        psi_temp[i]=psi[i]+dt*k_3[i];
    }
    Spatial_discretiser(psi_temp, k_4, omega);
    
    //calculating new psi after dt time increment
    for (int i=0; i<N; i++){
        psi[i]=psi[i]+dt/6.*(k_1[i]+2.*k_2[i]+2.*k_3[i]+k_4[i]);
    }
}

//spatially discretises RHS of coupled GP eqn in 1D using FDM and calculates slope k=dpsi/dt=-iMpsi(a0,b0,a1,b1,...,aN-1,bN-1) 
void GPsolver::Spatial_discretiser(dcomp psi_temp[], dcomp k[], double omega[]){

    dcomp C[N]; //constant introduced for convenience 
    Const_calc(psi_temp, C); //calculates constant for each component at each gridpoint

    //calculating -idpsi/dt for each gridpoint (4th order scheme), (N+i)%N to make grid loop
    for (int i=0; i<N; i++){
        if (i%2==0){ //even entries are for condensate a
            k[i]=-I*(((psi_temp[(N+i-4)%N]/4.)-(4.*psi_temp[(N+i-2)%N])-(4.*psi_temp[(N+i+2)%N])+(psi_temp[(N+i+4)%N]/4.))/(3.*pow(dx,2))+C[i]*psi_temp[i]+omega[i+1]*psi_temp[i+1]);
        }else{ //odd entries for condensate b
            k[i]=-I*(((psi_temp[(N+i-4)%N]/4.)-(4.*psi_temp[(N+i-2)%N])-(4.*psi_temp[(N+i+2)%N])+(psi_temp[(N+i+4)%N]/4.))/(3.*pow(dx,2))+C[i]*psi_temp[i]+omega[i-1]*psi_temp[i-1]);
        }
    }
}

//Calculates convenient constant for RHS of discretised coupled GP eqns C(a0,b0,a1,b1,...,aN-1,bN-1)
void GPsolver::Const_calc(dcomp psi_temp[], dcomp C[]){
    
    double mu; //chemical potential in units of gn
    Chem_potential(mu); //calculates mu

    for (int i=0; i<N; i++){
        if (i%2==0){ //even entries are for condensate a
            C[i]=5/(2*pow(dx,2))+norm(psi_temp[i])/n_0+g_ab*norm(psi_temp[i+1])/(g*n_0)-mu;
        }else{ //odd entries for condensate b
            C[i]=5/(2*pow(dx,2))+norm(psi_temp[i])/n_0+g_ab*norm(psi_temp[i-1])/(g*n_0)-mu;
        }
    }
}

//calculates dimensionless chemical potential mu, to make GP eqns timeless
void GPsolver::Chem_potential(double mu){
    
    mu=(1+g_ab/g)/2;
}





