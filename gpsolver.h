//Header file for GPsolver class functions and BEC parameters
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
        double omegaLHS=0, omegaRHS=0; //coherent coupling on both sides of discontinuity
        //excitation wavepacket parameters
        double x_0=100; //initial position of packet
        //double k_fundamental=2*M_PI/L; //fundamental wavevector of system
        double k_0=2;//80*k_fundamental; //wavevector of packet in terms of fundamental wavevector
        double width=5; //packet spatial width
        double A=0.01; //packet peak amplitude (1/100th of steady state)
};

//declaring functions
class GPsolver: public BEC_parameters { 
        public:
            //generates spatial grid and modulation of coherent coupling 
            void Gridspace(double x[]);

            //generating coherent coupling for modulation of system
            void Modulator(double omega[]);

            //calculates Bogoliubov positive (u) and negative (v) mode amplitudes for given k_0
            void Bogoliubov_mode_amplitudes(double &u, double &v);

            //generates initial psi
            void Init_psi_generator(dcomp psi[], bool excitation, double x[]);

            //solves eigenproblem (resulting from discretisation) using RK4 method to get psi(a0,b0,a1,b1,...,aN-1,bN-1) at +dt
            void RK4(dcomp psi[], double omega[]);
            
            //spatially discretises RHS of coupled GP eqn in 1D using FDM and calculates slope k=dpsi/dt=-iMpsi(a0,b0,a1,b1,...,aN-1,bN-1) 
            void Spatial_discretiser(dcomp psi_temp[], dcomp k[], double omega[]);
                        
            //Calculates convenient constant for RHS of discretised coupled GP eqns C(a0,b0,a1,b1,...,aN-1,bN-1) 
            void Const_calc(dcomp k[], dcomp C[]);

            //Calculates dimensionless chemical potential mu, to make GP eqns timeless
            void Chem_potential(double mu);
};

#endif