//Header file for FDM and Runge functions

//Spatially discretises the RHS of coupled GP eqns (for 2-component homogenous BEC) in 1D using finite differential method
//Returns RHS of coupled GP eqns in the form of a matrix M and coloumn vector x (the discretised order parameter) 

#ifndef GPSOLVER_H
#define GPSOLVER_H

class gpsolver {
        public:
            //spatially discretises RHS of coupled GP eqn  in 1D using FDM
            void spatialDiscretiser(int, int);

            //solves eigenproblem (resulting from discretisation) using Runge-Kutter method
            void temporalDiscretiser(int, int);

        private:
};


#endif