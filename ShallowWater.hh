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
        diff(1), interp(2), dissip(4){
        Grid f(nrows, ncols);
        for (size_t i = 0; i < nrows; i++)
            for (size_t j = 0; j < ncols; j++)
                f(i, j) = sp["f0"] + sp["beta"] * (j * dy - ylen / 2.);
        gp["f"] = f;
    }
    #define rotate90(var, i, j)  interp(interp(var[i], attr[i].hal, j), i) 
    void operator() (const State &var, State &dvar, float){
        phix = interp(var[0], attr[0].hal, 1);
        phiy = interp(var[0], attr[0].hal, 2);
        dvar[0] = - diff(var[1], 1) / dx - diff(var[2], 2) / dy;
        dvar[1] = - 0.5 * diff(var[0] * var[0], attr[0].hal * attr[0].hal, 1) / dx;
            /*
            + interp(
                gp["f"] * interp(var[2], 2)
                - diff(var[1] * var[1] / phix, 1) / dx
                - diff(var[2] * rotate90(var, 1, 2) / phiy, 2) / dy, 
                attr[1].hal, 1);
                */
        dvar[2] = - 0.5 * diff(var[0] * var[0], attr[0].hal * attr[0].hal, 2) / dy;
            /*
            + interp(
                - gp["f"] * interp(var[1], 1)
                - diff(var[1] * rotate90(var, 2, 1) / phix, 1) / dx
                - diff(var[2] * var[2] / phiy, 2) / dy,
                attr[2].hal, 2);
                */
        dvar[3] = - diff(var[1] * interp(var[3], attr[3].hal, 1) / phix, 1) / dx
            - diff(var[2] * interp(var[3], attr[3].hal, 2) / phiy, 2) / dy;
        //dvar[3] = - diff(adv.upwind(var[1], var[3], 1) / phix, 1) / dx
        //    - diff(adv.upwind(var[2], var[3], 2) / phiy, 2) / dy;
        /*
        for (size_t i = 0; i < 4; i++){
            dvar[i] += 0.03 / dt * (dissip(var[i], 1) + dissip(var[i], 2));
        }*/
        check_dimension(dvar);
    }
    #undef rotate90
    void set_boundary_conditions(){
        Halo phi, uwind, vwind, tracer;
        phi.set_periodic();
        uwind.set_periodic();
        vwind.set_periodic();
        tracer.set_periodic();

        attr.emplace_back("phi", 0, phi);
        attr.emplace_back("uwind", 1, uwind);
        attr.emplace_back("vwind", 2, vwind);
        attr.emplace_back("tracer", 0, tracer);
    }
protected:
    void ncwrite(float t){
        //std::cout << "Now writing..." << std::endl;
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
