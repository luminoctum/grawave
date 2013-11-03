#ifndef SLABHEATING
#define SLABHEATING
#include "ShallowWater.hh"
#include "MicroPhysics.hh"

class SlabHeating : public OdeSystem{
protected:
    Grid wwind, phi;
    float eps;
    Ammonia NH3;
    Water H2O;

    float PREF, grav, T0, Cp;
    Grid svp, ptol, Buffer;

    BiasedDifference diff;
    Interpolate interp;
    Integral integ;
    Diffusion dissip;
    Zalesak zadjust;

public:
    SlabHeating(std::string control_file): 
        OdeSystem(control_file),
        diff(nrows + 1, ncols + 1),
        interp(nrows + 1, ncols + 1, 6),
        dissip(nrows + 1, ncols + 1, 4),
        integ(nrows + 1, ncols + 1),
        zadjust(nrows, ncols){
        load_nc_file();
        eps = 8.7135;
        PREF = 1.E5;
        grav = 10.44;
        T0 = 134.8;
        Cp = 13947.26;
        ptol = PREF * gp["mass"];
        H2O.heat = Grid::Zero(nrows, ncols);
    }
    void operator() (const State &var, State &dvar, float){
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
        dvar[0] = sp["f0"] * var[1] - diff(gp["mass"] * interp(phi, 2), vattr[2].bc, 1) / dx
            + interp( 
                - diff(var[0] * var[0] / gp["massx"], 1) / dx 
                - diff(wwind * rotate90(var, 0, 2, 1) / gp["massy"], 2) / dy,
                vattr[0].bc, 1);
        dvar[1] = - sp["f0"] * var[0]
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
    void update(float t){
        var[5] = integ(diff(var[0], 1), vattr[5].bc, -2) * dy / dx;
        var[6] = integ(var[2] * gp["T_"] / T0 * (1 + var[3]) / (1 + eps * gp["eta_"]), 
                    vattr[5].bc, 2) * dy
               + integ(grav * gp["T_"] / T0 
                    * (1 - eps) / (1 + eps * var[3]) 
                    * (var[3] - gp["eta_"]) / (1 + eps * gp["eta_"]), 
                    vattr[6].bc, 2) * dy;
        // condense ammonia
        svp = NH3.sat_vapor_pressure(gp["T_"] * (1. + var[2] / grav));
        var[8] = var[4] / (1 + var[4]) * ptol / svp;
        var[8] = var[8].min(1.0).max(0.);
        var[4] = var[8] * svp / (ptol - var[8] * svp);
        // condese water
        svp = H2O.sat_vapor_pressure(gp["T_"] * (1. + var[2] / grav));
        H2O.heat = grav * (var[3] - svp / (ptol - svp)).max(0) * eps * H2O.Lf / (Cp * gp["T_"] * dt);
        var[7] = var[3] / (1 + var[3]) * ptol / svp;
        var[7] = var[7].min(1.0).max(0.);
        var[3] = var[7] * svp / (ptol - var[7] * svp);
    }
    void set_boundary_conditions(){
        Boundary uwind, wwind, buoyancy, etaH2O, etaNH3, phi, RH;

        uwind.left = uwind.right << Dirichlet | ZERO2(1, ncols);
        uwind.bottom = uwind.top << Neumann | ZERO2(nrows + 1, 1);
        wwind.left = wwind.right << Dirichlet | ZERO2(1, ncols + 1);
        wwind.bottom << Neumann | ZERO2(nrows, 1);
        wwind.top << Dirichlet | ZERO2(nrows, 1);
        buoyancy.left = buoyancy.right << Dirichlet | ZERO2(1, ncols);
        buoyancy.bottom << Neumann | ZERO2(nrows, 1);
        buoyancy.top << Dirichlet | ZERO2(nrows, 1);
        etaH2O.left = etaH2O.right << Neumann | ZERO2(1, ncols);
        etaH2O.bottom << Dirichlet | gp["eta_H2O"].col(0).maxCoeff() + ZERO2(nrows, 1);
        etaH2O.top << Neumann | ZERO2(nrows, 1);
        etaNH3.left = etaNH3.right << Neumann | ZERO2(1, ncols);
        etaNH3.bottom = etaNH3.top << Neumann | ZERO2(nrows, 1);
        phi.left = phi.right << Dirichlet | ZERO2(1, ncols + 1);
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
            if (i <= 1) Buffer = var[i] / gp["massx"];
            else if (i == 5) Buffer = var[i] / gp["massy"];
            else Buffer = var[i];
            dataFile.get_var(vattr[i].name.c_str())->put_rec(&Buffer(0, 0), ncfile.current);
        }
        dataFile.get_var("time")->put_rec(&t, ncfile.current);
        ncfile.current++;
    }
public:
    void set_NH3_mixr(float mixr){
        svp = NH3.sat_vapor_pressure(gp["T_"] * (1. + var[2] / grav));
        var[8] = 1. + Grid::Zero(nrows, ncols);
        var[4] = var[8] * svp / ptol;
        var[4] = var[4].min(mixr);
        var[8] = var[4] * ptol / svp;
        gp["etaNH3_"] = var[4];
    }
};

#endif
