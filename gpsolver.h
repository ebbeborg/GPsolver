//Header file for FDM and Runge functions

//Spatially discretises the RHS of coupled GP eqns (for 2-component homogenous BEC) in 1D using finite differential method
//Returns RHS of coupled GP eqns in the form of a matrix M and coloumn vector x (the discretised order parameter) 

#ifndef GPSOLVER_H
#define GPSOLVER_H

#include <complex>
typedef std::complex<double> dcomp; //defining complex data type

//storing relevant condensate parameters (in healing lengths)
class BEC_parameters {
    public:
        //discretisation parameters
        int gridsize=800; //number of points/nodes on our condensate 1D grid
        int N=2*gridsize; //size of psi[N] (twice the gridsize to accomodate both components a&b)
        double runtime=20; //total time
        double dt=0.001; //time stepsize
        //ground state homogenous symmetric condensate parameters (GS1)
        double L=200; //length of system
        double dx=L/gridsize; //grid spacing
        double n_0=1/L; //total density of 2D condensate n_a=n_b=n_0/2
        double V_a=0, V_b=0; //external potential for homogenous system
        double g=1, g_ab=0.8; //interaction constants
        double omega=0; //LHS=0, omegaRHS=0; //coherent coupling on both sides of discontinuity
        //excitation wavepacket parameters
        double x_0=100; //initial position of packet
        double k_0=1.5; //wavevector of packet 
        double width=5;//packet spatial width
        //double A=sqrt(n_0/2)/10; //packet peak amplitude (1/10th of steady state)
        double u=1.0023, v=-0.0677; //amplitude of positive and negative bogoliubov exctitation modes
};

//declaring functions
class GPsolver: public BEC_parameters { 
        public:
            //generates spatial grid and coherent coupling distribution
            void System_generator(double x[]);

            //generates initial psi
            void Init_psi_generator(dcomp psi[], bool excitation, double x[]);

            //solves eigenproblem (resulting from discretisation) using RK4 method to get psi(a0,b0,a1,b1,...,aN-1,bN-1) at +dt
            void RK4(dcomp psi[]);
            
            //spatially discretises RHS of coupled GP eqn in 1D using FDM and calculates slope k=dpsi/dt=-iMpsi(a0,b0,a1,b1,...,aN-1,bN-1) 
            void Spatial_discretiser(dcomp k[], dcomp Mk[]);
                        
            //Calculates convenient constant for RHS of discretised coupled GP eqns C(a0,b0,a1,b1,...,aN-1,bN-1) 
            void Const_calc(dcomp k[], dcomp C[]);

            //Calculates dimensionless chemical potential mu, to make GP eqns timeless
            void Chem_potential(double mu);
};

#endif