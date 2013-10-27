#ifndef CARTESIAN2D
#define CARTESIAN2D
#include "OdeSystem.hh"
#include "FiniteMethod.hh"

class Cartesian2d : public: OdeSystem{
protected:
    CenteredDifference diff;
    Diffusion dissip;

public:
    Cartesian2d(std::string control_file):
        OdeSystem(control_file),
        diff(nrows, ncols),
        dissip(nrows, ncols, 4){}

    void operator() (const StateType &var, StateType &dvar, float){
        dvar[0] = diff(var[1], vattr[1].bc, 1);
        dvar[1] = - gp["N2"] * diff(var[2], var[2].bc, 1);
        halo_update(dvar);
    }

    void update(float t){
        var[2] = poisson_solver(var[0], vattr[0]);
    }

protected:
    void set_boundary_conditions(){
        Boundary zeta, b, psi;
        zeta.left = zeta.right << Periodic;
        zeta.bottom = zeta.top << Dirichlet | ZERO2(nrows, 1);
        b.left = b.right << Periodic;
        b.bottom = b.top << Dirichlet | ZERO2(nrows, 1);
        psi.left = psi.right << Periodic;
        psi.bottom = psi.top << Dirichlet | ZERO2(nrows, 1);

        vattr.emplace_back("zeta", 0, zeta);
        vattr.emplace_back("buoyancy", 0, b);
        vattr.emplace_back("psi", 0, psi);
    }
}
