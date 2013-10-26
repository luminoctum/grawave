#ifndef ADVECTION
#define ADVECTION
#include "FiniteMethod.hh"

class Advection{
protected:
    Interpolate interp2;
    BiasedDifference diff;
    int nrows, ncols;
public:
    Advection(){}
    Advection(int _nrows, int _ncols) : 
        nrows(_nrows), ncols(_ncols),
        interp2(_nrows + 1, _ncols + 1, 2),
        diff(_nrows + 1, _ncols + 1){}

};

class Zalesak : public Advection{
protected:
    Grid var_td, var_min, var_max,
         flux_1, flux_h, flux_A, flux,
         RAp, RAm, RAdjust;
    float eps;
public:
    Zalesak(){}
    Zalesak(int _nrows, int _ncols) : Advection(_nrows, _ncols), eps(1.E-6),
        var_td(_nrows, _ncols), var_min(_nrows, _ncols), var_max(_nrows, _ncols),
        flux_1(_nrows + 1, _ncols + 1), flux_h(_nrows + 1, _ncols + 1), 
        flux_A(_nrows + 1, _ncols + 1), flux(_nrows, _ncols),
        RAp(_nrows + 2, _ncols + 2), RAm(_nrows + 2, _ncols + 2), 
        RAdjust(_nrows + 1, _ncols + 1){}

    Grid operator() (FiniteMethod &interph, const Grid &scale, 
            const Grid &wind, const Grid &var, const Boundary &varb, int axis){
        if (axis == 1){
            // Boundary condition not considered here
            flux_1.block(0,0,nrows+1,ncols) = wind * interp2(wind, var, varb, axis);
            var_td = var - scale * diff(flux_1.block(0,0,nrows+1,ncols), axis);
            flux_h.block(0,0,nrows+1,ncols) = wind * interph(var, varb, axis);
            flux_A.block(0,0,nrows+1,ncols) = (flux_h - flux_1).block(0,0,nrows+1,ncols);

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
            flux = flux_A.block(0,0,nrows,ncols).max(Grid::Zero(nrows, ncols)) 
                - flux_A.block(1,0,nrows,ncols).min(Grid::Zero(nrows, ncols));
            RAp.block(0,0,1,ncols) = Grid::Zero(1, ncols);
            RAp.block(1,0,nrows,ncols) = (flux.abs() < eps).select(0, 
                ((var_max - var_td) / (scale * flux)).min(Grid::Zero(nrows, ncols) + 1.));
            RAp.block(nrows+1,0,1,ncols) = Grid::Zero(1, ncols);
            flux = flux_A.block(1,0,nrows,ncols).max(Grid::Zero(nrows, ncols)) 
                - flux_A.block(0,0,nrows,ncols).min(Grid::Zero(nrows, ncols));
            RAm.block(0,0,1,ncols) = Grid::Zero(1, ncols);
            RAm.block(1,0,nrows,ncols) = (flux.abs() < eps).select(0,
                ((var_td - var_min) / (scale * flux)).min(Grid::Zero(nrows, ncols) + 1.));
            RAm.block(nrows+1,0,1,ncols) = Grid::Zero(1, ncols);

            RAdjust.block(0,0,nrows+1,ncols) = (flux_A.block(0,0,nrows+1,ncols) > 0).select(
                    RAm.block(0,0,nrows+1,ncols).min(RAp.block(1,0,nrows+1,ncols)),
                    RAp.block(0,0,nrows+1,ncols).min(RAm.block(1,0,nrows+1,ncols))
                    );
            return (flux_1 + RAdjust * flux_A).block(0,0,nrows+1,ncols);
        } else if (axis == 2){
            flux_1.block(0,0,nrows,ncols+1) = wind * interp2(wind, var, varb, axis);
            var_td = var - scale * diff(flux_1.block(0,0,nrows,ncols+1), axis);
            flux_h.block(0,0,nrows,ncols+1) = wind * interph(var, varb, axis);
            flux_A.block(0,0,nrows,ncols+1) = (flux_h - flux_1).block(0,0,nrows,ncols+1);

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
            flux = flux_A.block(0,0,nrows,ncols).max(Grid::Zero(nrows, ncols)) 
                - flux_A.block(0,1,nrows,ncols).min(Grid::Zero(nrows, ncols));
            RAp.block(0,0,nrows,1) = Grid::Zero(nrows, 1);
            RAp.block(0,1,nrows,ncols) = (flux.abs() < eps).select(0, 
                ((var_max - var_td) / (scale * flux)).min(Grid::Zero(nrows, ncols) + 1.));
            RAp.block(0,ncols+1,nrows,1) = Grid::Zero(nrows, 1);
            flux = flux_A.block(0,1,nrows,ncols).max(Grid::Zero(nrows, ncols)) 
                - flux_A.block(0,0,nrows,ncols).min(Grid::Zero(nrows, ncols));
            RAm.block(0,0,nrows,1) = Grid::Zero(nrows, 1);
            RAm.block(0,1,nrows,ncols) = (flux.abs() < eps).select(0,
                ((var_td - var_min) / (scale * flux)).min(Grid::Zero(nrows, ncols) + 1.));
            RAm.block(0,ncols+1,nrows,1) = Grid::Zero(nrows, 1);

            RAdjust.block(0,0,nrows,ncols+1) = (flux_A.block(0,0,nrows,ncols+1) > 0).select(
                    RAm.block(0,0,nrows,ncols+1).min(RAp.block(0,1,nrows,ncols+1)),
                    RAp.block(0,0,nrows,ncols+1).min(RAm.block(0,1,nrows,ncols+1))
                    );
            return (flux_1 + RAdjust * flux_A).block(0,0,nrows,ncols+1);
        } else {ERROR_1;}
    }
};

class Arakawa : public Advection{
protected:
    Grid Jpp, Jpx, Jxp, Buffer;
    CenteredDifference Diff;
public:
    Arakawa(){}
    Arakawa(int _nrows, int _ncols) : Advection(_nrows, _ncols),
        Jpp(_nrows, _ncols), Jpx(_nrows, _ncols), Jxp(_nrows, _ncols),
        Buffer(_nrows, _ncols){};

    Grid operator() (const Grid &psi, const Boundary &psib, const Grid &var, const Boundary &varb){
        Jpp = Diff(psi, psib, 1) * Diff(var, varb, 2)
            - Diff(psi, psib, 2) * Diff(var, varb, 1);
        Jpx = revolve(psi * Diff(var, varb, 2), -1)
            - revolve(psi * Diff(var, varb, 2), +1)
            - revolve(psi * Diff(var, varb, 1), -2)
            + revolve(psi * Diff(var, varb, 1), +2);
        Jxp = revolve(var * Diff(psi, psib, 1), -2)
            - revolve(var * Diff(psi, psib, 1), +2)
            - revolve(var * Diff(psi, psib, 2), -1)
            + revolve(var * Diff(psi, psib, 2), +1);
        return ((Jpp + Jxp + Jxp) / 3.);
    }
protected:
    Grid revolve(Grid var, int axis){
        if (axis == 1){
            for (size_t i = 0; i < var.rows(); i++)
                Buffer.row(i) = var.row((i - 1 + var.rows()) % var.rows());
        } else if (axis == -1){
            for (size_t i = 0; i < var.rows(); i++)
                Buffer.row(i) = var.row((i + 1) % var.rows());
        } else if (axis == 2){
            for (size_t i = 0; i < var.cols(); i++)
                Buffer.col(i) = var.col((i - 1 + var.cols()) % var.cols());
        } else if (axis == -2){
            for (size_t i = 0; i < var.cols(); i++)
                Buffer.col(i) = var.col((i + 1) % var.cols());
        } else{
            //raise error
        }
        return Buffer;
    };
};

#endif
