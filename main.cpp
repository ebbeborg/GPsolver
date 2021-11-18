#include "gpsolver.h"


using namespace std;

int main()
{

    //input spatial discretisation parameters
    
    //
    gpsolver GPa;
    gpsolver GPb;

    GPa.spatialDiscretiser();
    GPb.spatialDiscretiser();
    
    GPa.temporalDiscretiser();
    GPb.temporalDiscretiser();

    return 0;
}