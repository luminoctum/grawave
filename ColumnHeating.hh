#ifndef COLUMNHEATING
#define COLUMNHEATING
#include "ShallowWater.hh"
#include "SlabHeating.hh"
#include "MicroPhysics.hh"

class ColumnHeating: public SlabHeating{
protected:
    Grid vwind_over_r;
public:
    ColumnHeating(std::string control_file): 
        SlabHeating(control_file),
        vwind_over_r(nrows + 1, ncols){
        gp["mass0"] = gp["mass"];
        for (size_t i = 0; i < ncols; i++){
            gp["mass"].col(i) *= gp["west_east"].transpose();
            gp["massx"].col(i) *= gp["west_eastb"].transpose();
        }
        for (size_t i = 0; i < ncols + 1; i++)
            gp["massy"].col(i) *= gp["west_east"].transpose();
    }
    void operator() (const StateType &var, StateType &dvar, float){
    #define rotate90(var, i, j, k) ( interp(interp(var[i], vattr[i].bc, j), k) )
        for (size_t i = 0; i < dvar.size(); i++) 
            dvar[i] = ZERO2(var[i].rows(), var[i].cols());

        wwind = integ(diff(var[0], 1), vattr[5].bc, -2) * dy / dx;
        phi = integ(var[2] * gp["T_"] / T0 * (1 + var[3]) / (1 + eps * gp["eta_"]), 
                vattr[5].bc, 2) * dy
            + integ(grav * gp["T_"] / T0 
                * (1 - eps) / (1 + eps * var[3]) 
                * (var[3] - gp["eta_"]) / (1 + eps * gp["eta_"]), 
                vattr[6].bc, 2) * dy;
        for (size_t i = 0; i < ncols; i++)
            vwind_over_r.col(i) = (var[1] / gp["massx"]).col(i) / gp["west_eastb"].transpose();
        dvar[0] = (sp["f0"] + vwind_over_r) * var[1]
            - diff(gp["mass"] * interp(phi, 2), vattr[2].bc, 1) / dx
            + interp(gp["mass0"] * interp(phi, 2), vattr[2].bc, 1)
            + interp( 
                - diff(var[0] * var[0] / gp["massx"], 1) / dx 
                - diff(wwind * rotate90(var, 0, 2, 1) / gp["massy"], 2) / dy,
                vattr[0].bc, 1);
        dvar[1] = - (sp["f0"] + vwind_over_r) * var[0]
            + interp( 
                - diff(var[0] * var[1] / gp["massx"], 1) / dx 
                - diff(wwind * rotate90(var, 1, 2, 1) / gp["massy"], 2) / dy,
                vattr[1].bc, 1);
        dvar[2] = - 1 / gp["mass"] * interp(gp["N2"] * wwind, 2) + H2O.heat;
        for (size_t i = 2; i < 5; i++){
            dvar[i] += - 1 / gp["mass"] * (
                diff(zadjust(interp, dt / (dx * gp["mass"]), var[0], var[i], vattr[i].bc, 1), 1) / dx
                + diff(zadjust(interp, dt / (dx * gp["mass"]), wwind, var[i], vattr[i].bc, 2), 2) / dy
                );
        }
        for (size_t i = 0; i < 3; i++) 
            dvar[i] += 0.03 / dt * (dissip(var[i], 1) + dissip(var[i], 2));
        for (size_t i = 0; i < 2; i++) 
            dvar[i] += - gp["absorbx"] * var[i] / dt;
        dvar[2] += - gp["absorb"] * var[2] / dt;
        halo_update(dvar);
    #undef rotate90
    }
    void set_boundary_conditions(){
        Boundary uwind, wwind, buoyancy, etaH2O, etaNH3, phi, RH;

        uwind.left = uwind.right << Dirichlet | ZERO2(1, ncols);
        uwind.bottom = uwind.top << Neumann | ZERO2(nrows + 1, 1);
        wwind.left << Neumann | ZERO2(1, ncols + 1);
        wwind.right << Dirichlet | ZERO2(1, ncols + 1);
        wwind.bottom << Neumann | ZERO2(nrows, 1);
        wwind.top << Dirichlet | ZERO2(nrows, 1);
        buoyancy.left << Neumann | ZERO2(1, ncols);
        buoyancy.right << Dirichlet | ZERO2(1, ncols);
        buoyancy.bottom << Neumann | ZERO2(nrows, 1);
        buoyancy.top << Dirichlet | ZERO2(nrows, 1);
        etaH2O.left = etaH2O.right << Neumann | ZERO2(1, ncols);
        etaH2O.bottom << Dirichlet | gp["eta_H2O"].col(0).maxCoeff() + ZERO2(nrows, 1);
        etaH2O.top << Neumann | ZERO2(nrows, 1);
        etaNH3.left = etaNH3.right << Neumann | ZERO2(1, ncols);
        etaNH3.bottom = etaNH3.top << Neumann | ZERO2(nrows, 1);
        phi.left << Neumann | ZERO2(1, ncols + 1);
        phi.right << Dirichlet | ZERO2(1, ncols + 1);
        phi.bottom << Dirichlet | ZERO2(nrows, 1);
        phi.top << Neumann | ZERO2(nrows, 1);

        vattr.emplace_back("uwind",     1,  uwind);
        vattr.emplace_back("vwind",     1,  uwind);
        vattr.emplace_back("buoyancy",  0,  buoyancy);
        vattr.emplace_back("eta_H2O",   0,  etaH2O);
        vattr.emplace_back("eta_NH3",   0,  etaNH3);
        vattr.emplace_back("wwind",     2,  wwind);
        vattr.emplace_back("phi",       0,  phi);
        vattr.emplace_back("RH_H2O",    0,  RH);
        vattr.emplace_back("RH_NH3",    0,  RH);
    }
};

#endif
