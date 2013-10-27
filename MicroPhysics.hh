#ifndef MICROPHYSICS
#define MICROPHYSICS
#include <map>
#include <vector>
#include <cmath>
#include "Include.hh"
#include "Tools.hh"

class Species{
public:
    int nphase;
    std::vector<std::string> phase_name;
    std::map<std::string, Grid> phase;
    Grid heat, rh;
    float as, bs, al, bl, T3, p3, T00, Lf, Ls, mu, R;
    //std::string file_thermo_tbl;
public:
    Species(){};
    Species(std::string str){
        std::ifstream infile;
        std::string name;
        //file_thermo_tbl = "thermodynamics.tbl";
        infile.open("thermodynamics.tbl", std::ios::in);
        if (!infile) {ASSERT_FILE_NOT_FOUND("thermodynamics.tbl");}
        else {
            sclocate(infile, "Saturation vapor pressure");
            getline(infile, name);
            while(!infile.eof()){
                infile >> name >> as >> bs >> al >> bl
                    >> T3 >> p3 >> T00 >> Lf >> Ls;
                a_s[name] = as;
                b_s[name] = bs;
                a_l[name] = al;
                b_l[name] = bl;
                T_3[name] = T3;
                p_3[name] = p3;
                T_00[name] = T00;
                L_f[name] = Lf;
                L_s[name] = Ls;
            }
        }
        set_thermo_data(str);
    }
    void set_thermo_data(std::string name){
        as = a_s[name];
        bs = b_s[name];
        al = a_l[name];
        bl = b_l[name];
        T3 = T_3[name];
        p3 = p_3[name];
        T00 = T_00[name];
        Lf = L_f[name];
        Ls = L_s[name];
    }
    virtual float sat_vapor_pressure(float temp){}
    virtual Grid sat_vapor_pressure(const Grid &temp){}
private:
    std::map<std::string, float> a_s, b_s, a_l, b_l,
        T_3, p_3, T_00, L_f, L_s, M_u;
};

class Water : public Species{
public:
    Water(): Species("H2O"){
        mu = 18.0E-3;
        R = 8.31 / mu;}
    float sat_vapor_pressure(float temp){
        float t_c, omega, sat_vapor_p;
        if (temp <= T3) {
            sat_vapor_p = pow(10., as + bs / temp);
        }
        else if (temp <  T3 && temp >  T00) {
            omega = (temp - T3) / (T3 - T00);
            t_c = temp - 273.15;
            t_c = 17.67 * t_c / (t_c + 243.5);
            sat_vapor_p = 100 * omega * 6.112 * pow(10., t_c)
                + (1. - omega) * pow(10., al + bl / temp);
        }
        else {
            sat_vapor_p = pow(10., al + bl / temp);
        }
        return sat_vapor_p;
    }
    Grid sat_vapor_pressure(const Grid &temp){
        Grid Buffer = temp.max(85.).min(420.);
        return 611.2 * exp(17.67 * (Buffer - 273.15) / (Buffer - 29.65));
    }
};

class Ammonia : public Species{
public:
    Ammonia(): Species("NH3"){}
    float sat_vapor_pressure(float temp){
        float sat_vapor_p;
        if (temp < T00) {
            sat_vapor_p = pow(10., as + bs / temp);
        } else {
            sat_vapor_p = pow(10., al + bl / temp);
        }
        return sat_vapor_p;
    }
    Grid sat_vapor_pressure(const Grid &temp){
        Grid Buffer = temp.max(85.).min(420.);
        return (Buffer < T00).select(
                exp((as + bs / Buffer) * log(10.)),
                exp((al + bl / Buffer) * log(10.)));
    }
};

#endif
