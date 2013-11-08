#ifndef SHALLOWWATER
#define SHALLOWWATER
#include "OdeSystem.hh"
#include "FiniteMethod.hh"
#include "Advection.hh"

class ShallowWater : public OdeSystem{
protected:
    Grid phix, phiy, reynolds, buffer;
    Difference diff;
    DifferenceN dissip;
    Interpolate interp;
    Zalesak aflux;

public:
    ShallowWater(std::string control_file): 
        OdeSystem(control_file), 
        diff(1), interp(2), dissip(4), aflux(6){
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
    void operator() (const State &var, State &dvar, float){
        #define rotate90(var, i, j)  interp(interp(var[i], attr[i].hal, j), i) 
        phix = interp(var[0], attr[0].hal, 1);
        phiy = interp(var[0], attr[0].hal, 2);
        reynolds = interp(var[1], attr[1].hal, 2) * interp(var[2], attr[2].hal, 1) 
            / interp(var[0], attr[0].hal, 0);
        dvar[0] = - diff(var[1], 1) / dx - diff(var[2], 2) / dy;
        dvar[1] = - 0.5 * diff(var[0] * var[0], attr[0].hal * attr[0].hal, 1) / dx
            + gp["fx"] * rotate90(var, 2, 1)
            - diff.periodic(interp(var[1], 1) * interp(var[1], 1) / var[0], 1) / dx
            - diff(reynolds, 2) / dy;
        dvar[2] = - 0.5 * diff(var[0] * var[0], attr[0].hal * attr[0].hal, 2) / dy
            - gp["fy"] * rotate90(var, 1, 2)
            - diff.periodic(interp(var[2], 2) * interp(var[2], 2) / var[0], 2) / dy
            - diff(reynolds, 1) / dx;
        dvar[3] = - diff(var[1] * interp(var[3], attr[3].hal, 1) / phix, 1) / dx
            - diff(var[2] * interp(var[3], attr[3].hal, 2) / phiy, 2) / dy;
        //dvar[3] = - diff(aflux(dt / dx, var[1] / phix, var[3], attr[3].hal, 1), 1) / dx
        //    - diff(aflux(dt / dx, var[2] / phiy, var[3], attr[3].hal, 2), 2) / dy;
        for (size_t i = 0; i < 4; i++){
            dvar[i] += 0.03 / dt * (dissip(var[i], attr[i].hal, 1) + dissip(var[i], attr[i].hal, 2));
        };
        clean_up(dvar);
        #undef rotate90
    }
    void set_boundary_conditions(){
        Halo phi, uwind, vwind, tracer;
        /* double periodic
        phi.set_all(Periodic);
        uwind.set_all(Periodic);
        vwind.set_all(Periodic);
        tracer.set_all(Periodic);
        */
        /* channel
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
        */
        /* box */
        phi.set_all(Dirichlet);
        tracer.set_all(Dirichlet);
        uwind.set_row(Dirichlet);
        //uwind.set_row_ghost();
        uwind.set_col(Neumann);
        vwind.set_row(Neumann);
        vwind.set_col(Dirichlet);
        //vwind.set_col_ghost();

        attr.emplace_back("phi", 0, phi);
        attr.emplace_back("uwind", 1, uwind);
        attr.emplace_back("vwind", 2, vwind);
        attr.emplace_back("tracer", 0, tracer);
    }
protected:
    void ncwrite(float t){
        NcFile dataFile(ncfile.fname.c_str(),NcFile::Write);
        for (size_t i = 0; i < attr.size(); i++){
            buffer = (
                i == 0 ? var[0] : (
                attr[i].type == 0 ? var[i] / var[0] : (
                var[i] / interp(var[0], attr[0].hal, attr[i].type))));
            dataFile.get_var(attr[i].name.c_str())->put_rec(&buffer(0, 0), ncfile.current);
        }
        dataFile.get_var("time")->put_rec(&t, ncfile.current);
        ncfile.current++;
    }
};
#endif
