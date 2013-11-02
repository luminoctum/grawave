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
    Grid one = Grid::Zero(3, 4) + 1.;
    Grid a(3,4);
    Grid uwind(4, 4), vwind(3,5);
    Grid b, c, d, flux;
    a << 1,2,3,7,
        5,2,7,3,
        3,1,9,4;
    uwind << 1,2,-2,-1,
         1,2,-2,-1,
         1,2,-2,-1,
         1,2,-2,-1;
    vwind << 1,-2,3,4,-5,
          1,2,3,-4,5,
          1,-2,3,4,5;
    Halo bc;
    bc.set_left_right(Dirichlet, a.row(0));
    bc.set_bottom_top(Dirichlet, a.col(0));
    //bc.set_periodic();
    bc.update(a);
    Zalesak zadjust(6);
    Difference diff(1); 
    Interpolate interp(2);
    DifferenceN difn(2);
    Integral integ;

    //cout << bc << endl;
    cout << a << endl << endl;
    cout << interp.upwind(uwind, a, bc, 1) << endl << endl;
    //cout << uwind << endl << endl;
    //cout << diff(a, bc, 2) << endl << endl;
    //cout << diff(a, 2) << endl << endl;
    //cout << interp.upwind(uwind, a, bc, 1) << endl << endl;
    flux = zadjust(10 * one, uwind, a, bc, 1);
    b = a - diff(uwind * interp.upwind(uwind, a, bc, 1), 1);
    c = a - diff(uwind * interp(a, bc, 1), 1);
    d = a - diff(flux, 1);
    cout << b << endl << endl;
    cout << c << endl << endl;
    cout << d << endl << endl;
    cout << flux << endl << endl;

    /*
    cout << a.col(0) << endl;
    cout << Grid::Zero(3, 1) << endl;
    //a.row(1) = a.row(1) - bc.left.value;
    //a.col(1) = a.col(1) - bc.top.value;
    cout << a << endl;
    cout << bc << endl;
    cout << diff(a,bc,2) << endl;
    cout << diff(a,2) << endl;
    cout << endl << interp(a, bc, 2) << endl;
    cout << endl << interp(a, 2) << endl;
    cout << endl << fluxdiv(diff, interp, uwind.transpose(), a, bc, 1) << endl;
    */
}
