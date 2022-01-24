//Header file for FDM and Runge functions

//Spatially discretises the RHS of coupled GP eqns (for 2-component homogenous BEC) in 1D using finite differential method
//Returns RHS of coupled GP eqns in the form of a matrix M and coloumn vector x (the discretised order parameter) 

#ifndef GPSOLVER_H
#define GPSOLVER_H

#include <vector>

//typedef std::vector<std::vector<double>> matrix;

class gpsolver {
        public:
            
            //Finds constant for spatially discretised wavefunction at each grid point G
            double C(double psi[], double g, double g_ab, double n_0, double dx);
            
            //spatially discretises RHS of coupled GP eqn  in 1D using FDM
            void spatialDiscretiser(int N, double psi[], double g, double g_ab, double n_0, double dx , double omega);

            //solves eigenproblem (resulting from discretisation) using Runge-Kutter method
            void RK4(int N, double psi[], double g, double g_ab, double n_0, double dx, double omega);

        private:
};

#endif