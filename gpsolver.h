//Header file for FDM and Runge functions

//Spatially discretises the RHS of coupled GP eqns (for 2-component homogenous BEC) in 1D using finite differential method
//Returns RHS of coupled GP eqns in the form of a matrix M and coloumn vector x (the discretised order parameter) 

#ifndef GPSOLVER_H
#define GPSOLVER_H

#include <vector>;

typedef std::vector<std::vector<double>> matrix;

class gpsolver {
        public:
            
            //Finds constant for spatially discretised wavefunction at each grid point G
            double discretisedConstant(double, double, double, double);
            
            //spatially discretises RHS of coupled GP eqn  in 1D using FDM
            void spatialDiscretiser(int, matrix &, double, double, double, double, double);

            //solves eigenproblem (resulting from discretisation) using Runge-Kutter method
            void temporalDiscretiser(matrix &);

        private:
};

#endif