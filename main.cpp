//ihd/dt(x_a)=((-1/2m)(hbar^2)(Del^2) + g_a(x_a)^2 +g_ab(x_b)^2)x_a + omega(x_b)
//ihd/dt(x_b)=((-1/2m)(hbar^2)(Del^2) + g_b(x_b)^2 +g_ab(x_a)^2)x_b + omega^*(x_a)

#include "gpsolver.h"
#include <iostream>;

using namespace std;

int main()
{
    int N; //gridsize
    float delta_x; //grid spacing
    float n_0; //initial density of 2D condensate
    float V; //external potential
    float g, g_ab;
    float G_a, G_b; //constant introduced to write dimensionless discretised GP equation in more convenient way
    float coherentCoupling;
    gpsolver homoGP;

    //input BEC condensate parameters
    cout<<"Input initial density of condensate n_0:"; 
    cin>>n_0;
    cout<<"Input intraspecies (a on a, b on b) contact constant g:"; 
    cin>>g;
    cout<<"Input interspecies (a on b, b on a) contact constant g_ab:"; 
    cin>>g_ab;
    cout<<"Input coherent coupling:"; 
    cin>>coherentCoupling;
    //input discretisation parameters
    cout<<"Input gridsize:"; 
    cin>>N;
    cout<<"Input grid spacing:"; 
    cin>>delta_x;

    //declaring matrix where RHS of discretised GP eqn will live 
    float M[N][N]={};

    //dimensionalising constant G
    G_a=homoGP.discretisedConstant(g, g_ab, n_0, delta_x);
    G_b=homoGP.discretisedConstant(g, g_ab, n_0, delta_x);

    homoGP.spatialDiscretiser(N, M, g, n_0, G_a, G_b, coherentCoupling);
    homoGP.temporalDiscretiser(M);

    return 0;
}