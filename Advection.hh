#ifndef ADVECTION
#define ADVECTION
#include "FiniteMethod.hh"

class Zalesak{ 
protected:
    Grid var_td, var_min, var_max, flux,
         flux_1x, flux_1y, flux_hx, flux_hy,
         flux_ax, flux_ay, RAdjustx, RAdjusty,
         RApx, RApy, RAmx, RAmy;
    Interpolate interp2, interph;
    Difference diff;
    float eps;
    int order;
public:
    Zalesak(){}
    Zalesak(int _order) : eps(1.E-6), interp2(2), interph(_order), diff(1){}

    Grid operator() (const Grid &dtdx, const Grid &wind, const Grid &var, const Halo &hal, int axis){
        // Boundary condition needs to be checked
        int nrows = var.rows(), ncols = var.cols();
        Grid dvar;

        var_min.resize(nrows, ncols);
        var_max.resize(nrows, ncols);
        if (axis == 1){
            dvar.resize(nrows + 1, ncols);
            RAmx.resize(nrows + 2, ncols);
            RApx.resize(nrows + 2, ncols);

            flux_1x = wind * interp2.upwind(wind, var, hal, axis);
            var_td = var - dtdx * diff(flux_1x, axis);
            flux_hx = wind * interph(var, hal, axis);
            flux_ax = flux_hx - flux_1x;

            for (size_t i = 0; i < nrows; i++){
                if (i == 0){
                    var_max.row(i) = var.row(i).max(var.row(i + 1)).
                        max(var_td.row(i)).max(var_td.row(i + 1));
                    var_min.row(i) = var.row(i).min(var.row(i + 1)).
                        min(var_td.row(i)).min(var_td.row(i + 1));
                } else if (i == nrows - 1){
                    var_max.row(i) = var.row(i).max(var.row(i - 1)).
                        max(var_td.row(i)).max(var_td.row(i - 1));
                    var_min.row(i) = var.row(i).min(var.row(i - 1)).
                        min(var_td.row(i)).min(var_td.row(i - 1));
                } else {
                    var_max.row(i) = var.row(i).max(var.row(i - 1)).max(var.row(i + 1)).
                        max(var_td.row(i)).max(var_td.row(i - 1)).max(var_td.row(i + 1));
                    var_min.row(i) = var.row(i).min(var.row(i - 1)).min(var.row(i + 1)).
                        min(var_td.row(i)).min(var_td.row(i - 1)).min(var_td.row(i + 1));
                }
            }
            flux = flux_ax.block(0,0,nrows,ncols).max(ZERO2(nrows, ncols)) 
                - flux_ax.block(1,0,nrows,ncols).min(ZERO2(nrows, ncols));
            RApx.block(0,0,1,ncols) = ZERO2(1, ncols);
            RApx.block(1,0,nrows,ncols) = (flux.abs() < eps).select(0, 
                ((var_max - var_td) / (dtdx * flux)).min(ZERO2(nrows, ncols) + 1.));
            RApx.block(nrows+1,0,1,ncols) = ZERO2(1, ncols);
            flux = flux_ax.block(1,0,nrows,ncols).max(ZERO2(nrows, ncols)) 
                - flux_ax.block(0,0,nrows,ncols).min(ZERO2(nrows, ncols));
            RAmx.block(0,0,1,ncols) = ZERO2(1, ncols);
            RAmx.block(1,0,nrows,ncols) = (flux.abs() < eps).select(0,
                ((var_td - var_min) / (dtdx * flux)).min(ZERO2(nrows, ncols) + 1.));
            RAmx.block(nrows+1,0,1,ncols) = ZERO2(1, ncols);

            RAdjustx = (flux_ax > 0).select(
                    RAmx.block(0,0,nrows+1,ncols).min(RApx.block(1,0,nrows+1,ncols)),
                    RApx.block(0,0,nrows+1,ncols).min(RAmx.block(1,0,nrows+1,ncols))
                    );
            dvar = flux_1x + RAdjustx * flux_ax;
        } else if (axis == 2){
            dvar.resize(nrows, ncols + 1);
            RAmy.resize(nrows, ncols + 2);
            RApy.resize(nrows, ncols + 2);

            flux_1y = wind * interp2.upwind(wind, var, hal, axis);
            var_td = var - dtdx * diff(flux_1y, axis);
            flux_hy = wind * interph(var, hal, axis);
            flux_ay = flux_hy - flux_1y;

            for (size_t i = 0; i < ncols; i++){
                if (i == 0){
                    var_max.col(i) = var.col(i).max(var.col(i + 1)).
                        max(var_td.col(i)).max(var_td.col(i + 1));
                    var_min.col(i) = var.col(i).min(var.col(i + 1)).
                        min(var_td.col(i)).min(var_td.col(i + 1));
                } else if (i == ncols - 1){
                    var_max.col(i) = var.col(i).max(var.col(i - 1)).
                        max(var_td.col(i)).max(var_td.col(i - 1));
                    var_min.col(i) = var.col(i).min(var.col(i - 1)).
                        min(var_td.col(i)).min(var_td.col(i - 1));
                } else {
                    var_max.col(i) = var.col(i).max(var.col(i - 1)).max(var.col(i + 1)).
                        max(var_td.col(i)).max(var_td.col(i - 1)).max(var_td.col(i + 1));
                    var_min.col(i) = var.col(i).min(var.col(i - 1)).min(var.col(i + 1)).
                        min(var_td.col(i)).min(var_td.col(i - 1)).min(var_td.col(i + 1));
                }
            }
            flux = flux_ay.block(0,0,nrows,ncols).max(ZERO2(nrows, ncols)) 
                - flux_ay.block(0,1,nrows,ncols).min(ZERO2(nrows, ncols));
            RApy.block(0,0,nrows,1) = ZERO2(nrows, 1);
            RApy.block(0,1,nrows,ncols) = (flux.abs() < eps).select(0, 
                ((var_max - var_td) / (dtdx * flux)).min(ZERO2(nrows, ncols) + 1.));
            RApy.block(0,ncols+1,nrows,1) = ZERO2(nrows, 1);
            flux = flux_ay.block(0,1,nrows,ncols).max(ZERO2(nrows, ncols)) 
                - flux_ay.block(0,0,nrows,ncols).min(ZERO2(nrows, ncols));
            RAmy.block(0,0,nrows,1) = ZERO2(nrows, 1);
            RAmy.block(0,1,nrows,ncols) = (flux.abs() < eps).select(0,
                ((var_td - var_min) / (dtdx * flux)).min(ZERO2(nrows, ncols) + 1.));
            RAmy.block(0,ncols+1,nrows,1) = ZERO2(nrows, 1);

            RAdjusty = (flux_ay > 0).select(
                    RAmy.block(0,0,nrows,ncols+1).min(RApy.block(0,1,nrows,ncols+1)),
                    RApy.block(0,0,nrows,ncols+1).min(RAmy.block(0,1,nrows,ncols+1))
                    );
            dvar = flux_1y + RAdjusty * flux_ay;
        } else {ASSERT_VARIABLE_OUT_OF_RANGE("axis");}

        return dvar;
    }
    Grid operator() (float dtdx, const Grid &wind, const Grid &var, const Halo &hal, int axis){
        // Boundary condition needs to be checked
        int nrows = var.rows(), ncols = var.cols();
        Grid dvar;

        var_min.resize(nrows, ncols);
        var_max.resize(nrows, ncols);
        if (axis == 1){
            dvar.resize(nrows + 1, ncols);
            RAmx.resize(nrows + 2, ncols);
            RApx.resize(nrows + 2, ncols);

            flux_1x = wind * interp2.upwind(wind, var, hal, axis);
            var_td = var - dtdx * diff(flux_1x, axis);
            flux_hx = wind * interph(var, hal, axis);
            flux_ax = flux_hx - flux_1x;

            for (size_t i = 0; i < nrows; i++){
                if (i == 0){
                    var_max.row(i) = var.row(i).max(var.row(i + 1)).
                        max(var_td.row(i)).max(var_td.row(i + 1));
                    var_min.row(i) = var.row(i).min(var.row(i + 1)).
                        min(var_td.row(i)).min(var_td.row(i + 1));
                } else if (i == nrows - 1){
                    var_max.row(i) = var.row(i).max(var.row(i - 1)).
                        max(var_td.row(i)).max(var_td.row(i - 1));
                    var_min.row(i) = var.row(i).min(var.row(i - 1)).
                        min(var_td.row(i)).min(var_td.row(i - 1));
                } else {
                    var_max.row(i) = var.row(i).max(var.row(i - 1)).max(var.row(i + 1)).
                        max(var_td.row(i)).max(var_td.row(i - 1)).max(var_td.row(i + 1));
                    var_min.row(i) = var.row(i).min(var.row(i - 1)).min(var.row(i + 1)).
                        min(var_td.row(i)).min(var_td.row(i - 1)).min(var_td.row(i + 1));
                }
            }
            flux = flux_ax.block(0,0,nrows,ncols).max(ZERO2(nrows, ncols)) 
                - flux_ax.block(1,0,nrows,ncols).min(ZERO2(nrows, ncols));
            RApx.block(0,0,1,ncols) = ZERO2(1, ncols);
            RApx.block(1,0,nrows,ncols) = (flux.abs() < eps).select(0, 
                ((var_max - var_td) / (dtdx * flux)).min(ZERO2(nrows, ncols) + 1.));
            RApx.block(nrows+1,0,1,ncols) = ZERO2(1, ncols);
            flux = flux_ax.block(1,0,nrows,ncols).max(ZERO2(nrows, ncols)) 
                - flux_ax.block(0,0,nrows,ncols).min(ZERO2(nrows, ncols));
            RAmx.block(0,0,1,ncols) = ZERO2(1, ncols);
            RAmx.block(1,0,nrows,ncols) = (flux.abs() < eps).select(0,
                ((var_td - var_min) / (dtdx * flux)).min(ZERO2(nrows, ncols) + 1.));
            RAmx.block(nrows+1,0,1,ncols) = ZERO2(1, ncols);

            RAdjustx = (flux_ax > 0).select(
                    RAmx.block(0,0,nrows+1,ncols).min(RApx.block(1,0,nrows+1,ncols)),
                    RApx.block(0,0,nrows+1,ncols).min(RAmx.block(1,0,nrows+1,ncols))
                    );
            dvar = flux_1x + RAdjustx * flux_ax;
        } else if (axis == 2){
            dvar.resize(nrows, ncols + 1);
            RAmy.resize(nrows, ncols + 2);
            RApy.resize(nrows, ncols + 2);

            flux_1y = wind * interp2.upwind(wind, var, hal, axis);
            var_td = var - dtdx * diff(flux_1y, axis);
            flux_hy = wind * interph(var, hal, axis);
            flux_ay = flux_hy - flux_1y;

            for (size_t i = 0; i < ncols; i++){
                if (i == 0){
                    var_max.col(i) = var.col(i).max(var.col(i + 1)).
                        max(var_td.col(i)).max(var_td.col(i + 1));
                    var_min.col(i) = var.col(i).min(var.col(i + 1)).
                        min(var_td.col(i)).min(var_td.col(i + 1));
                } else if (i == ncols - 1){
                    var_max.col(i) = var.col(i).max(var.col(i - 1)).
                        max(var_td.col(i)).max(var_td.col(i - 1));
                    var_min.col(i) = var.col(i).min(var.col(i - 1)).
                        min(var_td.col(i)).min(var_td.col(i - 1));
                } else {
                    var_max.col(i) = var.col(i).max(var.col(i - 1)).max(var.col(i + 1)).
                        max(var_td.col(i)).max(var_td.col(i - 1)).max(var_td.col(i + 1));
                    var_min.col(i) = var.col(i).min(var.col(i - 1)).min(var.col(i + 1)).
                        min(var_td.col(i)).min(var_td.col(i - 1)).min(var_td.col(i + 1));
                }
            }
            flux = flux_ay.block(0,0,nrows,ncols).max(ZERO2(nrows, ncols)) 
                - flux_ay.block(0,1,nrows,ncols).min(ZERO2(nrows, ncols));
            RApy.block(0,0,nrows,1) = ZERO2(nrows, 1);
            RApy.block(0,1,nrows,ncols) = (flux.abs() < eps).select(0, 
                ((var_max - var_td) / (dtdx * flux)).min(ZERO2(nrows, ncols) + 1.));
            RApy.block(0,ncols+1,nrows,1) = ZERO2(nrows, 1);
            flux = flux_ay.block(0,1,nrows,ncols).max(ZERO2(nrows, ncols)) 
                - flux_ay.block(0,0,nrows,ncols).min(ZERO2(nrows, ncols));
            RAmy.block(0,0,nrows,1) = ZERO2(nrows, 1);
            RAmy.block(0,1,nrows,ncols) = (flux.abs() < eps).select(0,
                ((var_td - var_min) / (dtdx * flux)).min(ZERO2(nrows, ncols) + 1.));
            RAmy.block(0,ncols+1,nrows,1) = ZERO2(nrows, 1);

            RAdjusty = (flux_ay > 0).select(
                    RAmy.block(0,0,nrows,ncols+1).min(RApy.block(0,1,nrows,ncols+1)),
                    RApy.block(0,0,nrows,ncols+1).min(RAmy.block(0,1,nrows,ncols+1))
                    );
            dvar = flux_1y + RAdjusty * flux_ay;
        } else {ASSERT_VARIABLE_OUT_OF_RANGE("axis");}

        return dvar;
    }
};

class Arakawa{
protected:
    Grid Jpp, Jpx, Jxp, buffer;
    Difference diff;
public:
    Arakawa(){}
    Arakawa(int _nrows, int _ncols, int _order) : diff(1),
        Jpp(_nrows, _ncols), Jpx(_nrows, _ncols), Jxp(_nrows, _ncols){};

    Grid operator() (const Grid &psi, const Halo &psib, const Grid &var, const Halo &varb){
        Jpp = diff(psi, psib, 1) * diff(var, varb, 2)
            - diff(psi, psib, 2) * diff(var, varb, 1);
        Jpx = revolve(psi * diff(var, varb, 2), -1)
            - revolve(psi * diff(var, varb, 2), +1)
            - revolve(psi * diff(var, varb, 1), -2)
            + revolve(psi * diff(var, varb, 1), +2);
        Jxp = revolve(var * diff(psi, psib, 1), -2)
            - revolve(var * diff(psi, psib, 1), +2)
            - revolve(var * diff(psi, psib, 2), -1)
            + revolve(var * diff(psi, psib, 2), +1);
        return ((Jpp + Jxp + Jxp) / 3.);
    }
protected:
    Grid revolve(Grid var, int axis){
        buffer.resize(var.rows(), var.cols());
        if (axis == 1){
            for (size_t i = 0; i < var.rows(); i++)
                buffer.row(i) = var.row((i - 1 + var.rows()) % var.rows());
        } else if (axis == -1){
            for (size_t i = 0; i < var.rows(); i++)
                buffer.row(i) = var.row((i + 1) % var.rows());
        } else if (axis == 2){
            for (size_t i = 0; i < var.cols(); i++)
                buffer.col(i) = var.col((i - 1 + var.cols()) % var.cols());
        } else if (axis == -2){
            for (size_t i = 0; i < var.cols(); i++)
                buffer.col(i) = var.col((i + 1) % var.cols());
        } else{
            //raise error
        }
        return buffer;
    };
};

#endif
