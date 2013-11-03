#ifndef SHALLOWWATER
#define SHALLOWWATER
#include "OdeSystem.hh"
#include "FiniteMethod.hh"
#include "Advection.hh"

class ShallowWater : public OdeSystem{
protected:
    Grid phix, phiy, buffer;
    Difference diff;
    DifferenceN dissip;
    Interpolate interp;

public:
    ShallowWater(std::string control_file): 
        OdeSystem(control_file), 
        diff(1), interp(2), dissip(4),
        phix(nrows + 1, ncols), phiy(nrows, ncols + 1){
        Grid f(nrows, ncols);
        for (size_t i = 0; i < nrows; i++)
            for (size_t j = 0; j < ncols; j++)
                f(i, j) = sp["f0"] + sp["beta"] * (j * dy - ylen / 2.);
        gp["f"] = f;
        set_boundary_conditions();
    }
    #define rotate90(var, i, j)  interp(interp(var[i], vattr[i].bc, j), i) 
    void operator() (const State &var, State &dvar, float){
        phix = interp(var[0], vattr[0].bc, 1);
        phiy = interp(var[0], vattr[0].bc, 2);
        dvar[0] = - diff(var[1], 1) / dx - diff(var[2], 2) / dy;
        dvar[1] = - 0.5 * diff(var[0] * var[0], vattr[0].bc, 1) / dx
            + interp(
                gp["f"] * interp(var[2], 2)
                - diff(var[1] * var[1] / phix, 1) / dx
                - diff(var[2] * rotate90(var, 1, 2) / phiy, 2) / dy, 
                vattr[1].bc, 1);
        dvar[2] = - 0.5 * diff(var[0] * var[0], vattr[0].bc, 2) / dy
            + interp(
                - gp["f"] * interp(var[1], 1)
                - diff(var[1] * rotate90(var, 2, 1) / phix, 1) / dx
                - diff(var[2] * var[2] / phiy, 2) / dy,
                vattr[2].bc, 2);
        dvar[3] = - diff(var[1] * interp(var[3], vattr[3].bc, 1) / phix, 1) / dx
            - diff(var[2] * interp(var[3], vattr[3].bc, 2) / phiy, 2) / dy;
        //dvar[3] = - diff(adv.upwind(var[1], var[3], 1) / phix, 1) / dx
        //    - diff(adv.upwind(var[2], var[3], 2) / phiy, 2) / dy;
        /*
        for (size_t i = 0; i < 4; i++){
            dvar[i] += 0.03 / dt * (dissip(var[i], 1) + dissip(var[i], 2));
        }*/
        halo_update(dvar);
    }
    #undef rotate90
protected:
    void set_boundary_conditions(){
        Halo phi, uwind, vwind, tracer;
        phi.set_periodic();
        uwins.set_periodic();
        vwins.set_periodic();
        tracer.set_periodic();
        /*
        phi.left = phi.right << Neumann | ZERO2(1, ncols);
        phi.bottom = phi.top << Neumann | ZERO2(nrows, 1);
        uwind.left = uwind.right << Dirichlet | ZERO2(1, ncols);
        uwind.bottom = uwind.top << Neumann | ZERO2(nrows + 1, 1);
        vwind.left = vwind.right << Neumann | ZERO2(1, ncols + 1);
        vwind.bottom = vwind.top << Dirichlet | ZERO2(nrows, 1);
        tracer.left = tracer.right << Neumann | ZERO2(1, ncols);
        tracer.bottom = tracer.top << Neumann | ZERO2(nrows, 1);
        */

        vattr.emplace_back("phi", 0, phi);
        vattr.emplace_back("uwind", 1, uwind);
        vattr.emplace_back("vwind", 2, vwind);
        vattr.emplace_back("tracer", 0, tracer);
    }
    void halo_update(State &dvar){
        for (size_t i = 0; i < dvar.size(); i++){
            if (dvar[i].size() == 0)
                dvar[i] = ZERO2(var[i].rows(), var[i].cols());
            if (dvar[i].rows() != var[i].rows() || dvar[i].cols() != var[i].cols()){
                // raise error
                std::cerr << "Error 1" << std::endl;
                exit(-1);
            }
            if (vattr[i].bc.left.type == Dirichlet)
                dvar[i].row(0) = ZERO2(1, var[i].cols());
            if (vattr[i].bc.right.type == Dirichlet)
                dvar[i].row(var[i].rows() - 1) = ZERO2(1, var[i].cols());
            if (vattr[i].bc.bottom.type == Dirichlet)
                dvar[i].col(0) = ZERO2(var[i].rows(), 1);
            if (vattr[i].bc.top.type == Dirichlet)
                dvar[i].col(var[i].cols() - 1) = ZERO2(var[i].rows(), 1);
        }
    }
    void ncwrite(float t){
        //std::cout << "Now writing..." << std::endl;
        NcFile dataFile(ncfile.fname.c_str(),NcFile::Write);
        for (size_t i = 0; i < vattr.size(); i++){
            buffer = (
                i == 0 ? var[0] : (
                vattr[i].type == 0 ? var[i] / var[0] : (
                var[i] / interp(var[0], vattr[0].bc, vattr[i].type))));
            dataFile.get_var(vattr[i].name.c_str())->put_rec(&buffer(0, 0), ncfile.current);
        }
        dataFile.get_var("time")->put_rec(&t, ncfile.current);
        ncfile.current++;
    }
};
#endif
