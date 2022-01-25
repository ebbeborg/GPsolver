//Header file for FDM and Runge functions

//Spatially discretises the RHS of coupled GP eqns (for 2-component homogenous BEC) in 1D using finite differential method
//Returns RHS of coupled GP eqns in the form of a matrix M and coloumn vector x (the discretised order parameter) 

#ifndef GPSOLVER_H
#define GPSOLVER_H

//organising relevant condensate parameters
class BEC_parameters {
    public:
        //declaring parameters
        double gridsize; //number of points/nodes on our condensate 1D grid
        double N; //size of our wavefunction variable
        double dx; //grid spacing
        double n_0; //initial density of 2D condensate
        double V; //external potential
        double g, g_ab; //interaction constants
        double omega; //coherent coupling
};

class GPsolver {
        public:
            
            //Finds constant for spatially discretised wavefunction at each grid point 
            double C(double k[], BEC_parameters& parameters);
            
            //spatially discretises RHS of coupled GP eqn  in 1D using FDM andreturns -iMk
            double spatialDiscretiser(double k[], BEC_parameters& parameters);
            
            //solves eigenproblem (resulting from discretisation) using Runge-Kutter method
            void RK4(double psi_init[], BEC_parameters& parameters);

        private:
};

#endif