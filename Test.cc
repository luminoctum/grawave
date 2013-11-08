#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <vector>
//#include "OdeSystem.hh"
#include "FiniteMethod.hh"
#include "Advection.hh"

using namespace std;
using namespace Eigen;

int main(){
    srand(0);
    ArrayXXf a = ArrayXXf::Random(10, 10);
    ArrayXXf u = ArrayXXf::Random(11, 10);
    ArrayXXf v = ArrayXXf::Random(10, 11);
    Halo bc;
    Difference diff(1);
    Interpolate interp(2);
    Zalesak advect(6);
    //bc.set_periodic();
    bc.set_row(Dirichlet, a.row(0));
    bc.set_col(Dirichlet, a.col(0));
    bc.update(a);

    cout << a << endl << endl;
    cout << u << endl << endl;
    cout << v << endl << endl;
    cout << advect(1, u, a, bc, 1) << endl << endl;
    cout << advect(1, v, a, bc, 2) << endl;
}
