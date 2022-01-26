//Header file for FDM and Runge functions

//Spatially discretises the RHS of coupled GP eqns (for 2-component homogenous BEC) in 1D using finite differential method
//Returns RHS of coupled GP eqns in the form of a matrix M and coloumn vector x (the discretised order parameter) 

#ifndef GPSOLVER_H
#define GPSOLVER_H

//storing relevant condensate parameters
class BEC_parameters {
    public:
        //discretisation parameters
        int gridsize; //number of points/nodes on our condensate 1D grid
        int N=2*gridsize; //size of our wavefunction variable (twice the gridsize to accomodate both components a&b)
        double dx; //grid spacing
        //double dt=1; //timestep set to 1 right now
        double runtime; //total time
        //condensate parameters
        double n_0; //initial density of 2D condensate
        double V; //external potential
        double g, g_ab; //interaction constants
        double omega; //coherent coupling
};

class GPsolver: public BEC_parameters { 
        public:
            
            //Finds constant for spatially discretised wavefunction at each grid point 
            double C(double k[]);
            
            //spatially discretises RHS of coupled GP eqn  in 1D using FDM andreturns -iMk
            double spatialDiscretiser(double k[]);
            
            //solves eigenproblem (resulting from discretisation) using Runge-Kutter method
            void RK4(double psi[]);
};

#endif