#ifndef SHALLOWWATER
#define SHALLOWWATER
#include "OdeSystem.hh"
#include "FiniteMethod.hh"
//#include "Advection.hh"

class ShallowWater : public OdeSystem{
protected:
    Grid phix, phiy, reynolds, buffer;
    Difference<1> diff;
    Difference<2> diff2;
    DifferenceN<4> del2;
    Interpolate<2> half;

public:
    ShallowWater() : OdeSystem(){
        // initialize constants
        Grid fx(nrows + 1, ncols);
        Grid fy(nrows, ncols + 1);
        for (size_t i = 0; i < nrows + 1; i++)
            for (size_t j = 0; j < ncols; j++)
                fx(i, j) = sp["f0"] + sp["beta"] * (j * dy - ylen / 2.);
        for (size_t i = 0; i < nrows; i++)
            for (size_t j = 0; j < ncols + 1; j++)
                fy(i, j) = sp["f0"] + sp["beta"] * ((j - 0.5) * dy - ylen / 2.);
        gp["fx"] = fx;
        gp["fy"] = fy;
    }
    void set_boundary_conditions(){
        Halo phi, uwind, vwind, tracer;
        /* double periodic
        phi.set_all(Periodic);
        uwind.set_all(Periodic);
        vwind.set_all(Periodic);
        tracer.set_all(Periodic);
        */
        /* channel */
        phi.set_row(Periodic);
        phi.set_col(Neumann);
        uwind.set_row(Periodic);
        uwind.set_row_ghost();
        uwind.set_col(Neumann);
        vwind.set_row(Periodic);
        vwind.set_col(Dirichlet);
        vwind.set_col_ghost();
        tracer.set_col(Periodic);
        tracer.set_row(Periodic);
        /* box
        phi.set_all(Dirichlet);
        tracer.set_all(Dirichlet);
        uwind.set_row(Dirichlet);
        //uwind.set_row_ghost();
        uwind.set_col(Neumann);
        vwind.set_row(Neumann);
        vwind.set_col(Dirichlet);
        //vwind.set_col_ghost();
        */

        attr.emplace_back("phi", 0, phi);
        attr.emplace_back("uwind", 1, uwind);
        attr.emplace_back("vwind", 2, vwind);
        attr.emplace_back("tracer", 0, tracer);
    }
    void operator() (const State &var, State &dvar, float){
        phix = half.x(var[0], attr[0].hal);
        phiy = half.y(var[0], attr[0].hal);
        reynolds = half.y(var[1], attr[1].hal) * half.x(var[2], attr[2].hal) 
            / half.quad(var[0], attr[0].hal);
        dvar[0] = - diff.x(var[1]) / dx - diff.y(var[2]) / dy;
        dvar[1] = - 0.5 * diff.x(var[0] * var[0], attr[0].hal * attr[0].hal) / dx
            + gp["fx"] * half.y(half.x(var[2], attr[2].hal))
            - diff2.x(var[1] * var[1] / half.x(var[0], attr[0].hal)) / dx
            - diff.y(reynolds) / dy;
        dvar[2] = - 0.5 * diff.y(var[0] * var[0], attr[0].hal * attr[0].hal) / dy
            - gp["fy"] * half.x(half.y(var[1], attr[1].hal))
            - diff2.y(var[2] * var[2] / half.y(var[0], attr[0].hal)) / dy
            - diff.x(reynolds) / dx;
        dvar[3] = - diff.x(var[1] * half.x(var[3], attr[3].hal) / phix) / dx
            - diff.y(var[2] * half.y(var[3], attr[3].hal) / phiy) / dy;
        for (size_t i = 0; i < 4; i++){
            dvar[i] += 0.03 / dt * (del2.x(var[i], attr[i].hal) + del2.y(var[i], attr[i].hal));
        };
        clean_up(dvar);
    }
};
#endif
