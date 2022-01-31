//Header file for FDM and Runge functions

//Spatially discretises the RHS of coupled GP eqns (for 2-component homogenous BEC) in 1D using finite differential method
//Returns RHS of coupled GP eqns in the form of a matrix M and coloumn vector x (the discretised order parameter) 

#ifndef GPSOLVER_H
#define GPSOLVER_H

#include <complex>
typedef std::complex<double> dcomp; //defining complex data type

//storing relevant condensate parameters
class BEC_parameters {
    public:
        //discretisation parameters
        int gridsize=100; //number of points/nodes on our condensate 1D grid
        int N=2*gridsize; //size of wavefunction variable (twice the gridsize to accomodate both components a&b)
        double dx=L/gridsize; //grid spacing
        double runtime=10; //total time
        double dt=0.01; //time stepsize
        //condensate parameters
        double L=20;
        double n_0=1/L; //total density of 2D condensate n_a=n_b=n_0/2, 10^6 /m
        double V_a=0, V_b=0; //external potential for homogenous system
        double g=1, g_ab=0.1; //interaction constants
        double omega=0; //coherent coupling
};

//declaring functions
class GPsolver: public BEC_parameters { 
        public:
            //generates initial psi
            void Init_psi_generator(dcomp psi[]);

            //solves eigenproblem (resulting from discretisation) using RK4 method to get psi(a0,b0,a1,b1,...,aN-1,bN-1) at +dt
            void RK4(dcomp psi[]);
            
            //spatially discretises RHS of coupled GP eqn in 1D using FDM and calculates Mk(a0,b0,a1,b1,...,aN-1,bN-1)  
            void spatialDiscretiser(dcomp k[], dcomp Mk[]);
                        
            //Calculates convenient constant for spatial discretisation C(a0,b0,a1,b1,...,aN-1,bN-1) 
            void Const_calc(dcomp k[], double C[]);
};

#endif