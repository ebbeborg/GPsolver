#include "gpsolver.h"


using namespace std;

int main()
{
    int parameter1, parameter2;
    //input spatial discretisation parameters
    
    //
    gpsolver GPa;
    gpsolver GPb;

    GPa.spatialDiscretiser(parameter1, parameter2);
    GPb.spatialDiscretiser(parameter1, parameter2);
    
    GPa.temporalDiscretiser(parameter1, parameter2);
    GPb.temporalDiscretiser(parameter1, parameter2);

    return 0;
}