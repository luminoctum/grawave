#ifndef FINITEMethod
#define FINITEMethod

#include "Include.hh"
#include "Halo.hh"
#define MIDDLE(order, var, i) ( \
        order == 1 ? - 0.5 * (var(i + 1) - var(i)) : ( \
        order == 2 ? \
            static_cast<Grid>(0.5 * (var(i + 1) + var(i))) : ( \
        order == 3 ? \
            static_cast<Grid>(1./12. * (var(i + 2) - var(i - 1)) \
                - 1./4. * (var(i + 1) - var(i))) : ( \
        order == 4 ? \
            static_cast<Grid>(7./12. * (var(i + 1) + var(i)) \
                - 1./12. * (var(i + 2) + var(i - 1))) : ( \
        order == 5 ? \
            static_cast<Grid>(-1./60. * (var(i + 3) - var(i - 2)) \
                + 1./12. * (var(i + 2) - var(i - 1)) \
                - 1./6. * (var(i + 1) - var(i))) : \
        /* order == 6 */ \
            static_cast<Grid>(37./60. * (var(i + 1) + var(i)) \
                - 2./15. * (var(i + 2) + var(i - 1)) \
                + 1./60. * (var(i + 3) + var(i - 2))))))) \
        )
#define DIFFN(order, var, i) ( \
        order == 2 ? var(i-1)-2*var(i)+var(i+1) : ( \
        order == 4 ? \
            static_cast<Grid>(-var(i-2)+4*var(i-1)-6*var(i)+4*var(i+1)-var(i+2)) : \
        /* order == 6 */ \
            static_cast<Grid>(var(i-3)-6*var(i-2)+15*var(i-1)-20*var(i)+15*var(i+1)-6*var(i+2)+var(i+3))) \
        )

/* axis = 2  : up
 * axis = 1  : right
 * axis = -2 : down
 * axis = -1 : left
 */

template <int order>
class Difference{
    /* Calculate finite difference */
public:
    Grid x(const Grid &var){
        int nrows = var.rows(), ncols = var.cols();
        Grid dvar;
        if (order == 2){
            dvar.resize(nrows, ncols);
            dvar.row(0) = var.row(1) - var.row(0);
            for (size_t i = 1; i < nrows - 1; i++)
                dvar.row(i) = 0.5 * (var.row(i + 1) - var.row(i - 1));
            dvar.row(nrows - 1) = var.row(nrows - 1) - var.row(nrows - 2);
        } else if (order == 1){
            dvar.resize(nrows - 1, ncols);
            for (size_t i = 0; i < nrows - 1; i++)
                dvar.row(i) = var.row(i + 1) - var.row(i);
        } 
        return dvar;
    }
    Grid y(const Grid &var){
        int nrows = var.rows(), ncols = var.cols();
        Grid dvar;
        if (order == 2){
            dvar.resize(nrows, ncols);
            dvar.col(0) = var.col(1) - var.col(0);
            for (size_t i = 1; i < var.cols() - 1; i++)
                dvar.col(i) = 0.5 * (var.col(i + 1) - var.col(i - 1));
            dvar.col(ncols - 1) = var.col(ncols - 1) - var.col(ncols - 2);
        } else if (order == 1){
            dvar.resize(nrows, ncols - 1);
            for (size_t i = 0; i < ncols - 1; i++)
                dvar.col(i) = var.col(i + 1) - var.col(i);
        } 
        return dvar;
    }
    Grid x(const Grid &var, const Halo &hal){
        int nrows = var.rows(), ncols = var.cols();
        Grid dvar;
        if (order == 2){
            dvar.resize(nrows, ncols);
            dvar = this->x(var);
            dvar.row(0) = (var.row(1) - hal.left) / 2.;
            dvar.row(nrows - 1) = (hal.right - var.row(nrows - 2)) / 2.;
        } else if (order == 1){
            dvar.resize(nrows + 1, ncols);
            dvar.block(1, 0, nrows - 1, ncols) = this->x(var);
            dvar.row(0) = var.row(0) - hal.left;
            dvar.row(nrows) = hal.right - var.row(nrows - 1);
        } 
        return dvar;
    }
    Grid y(const Grid &var, const Halo &hal){
        int nrows = var.rows(), ncols = var.cols();
        Grid dvar;
        if (order == 2){
            dvar.resize(nrows, ncols);
            dvar = this->y(var);
            dvar.col(0) = (var.col(1) - hal.bottom) / 2.;
            dvar.col(ncols - 1) = (hal.top - var.col(ncols - 2)) / 2.;
        } else if (order == 1){
            dvar.resize(nrows, ncols + 1);
            dvar.block(0, 1, nrows, ncols - 1) = this->y(var);
            dvar.col(0) = var.col(0) - hal.bottom;
            dvar.col(ncols) = hal.top - var.col(ncols - 1);
        } 
        return dvar;
    }
};

template <int order>
class DifferenceN{
    /* calculate higher order difference */
public:
    Grid x(const Grid &var, const Halo &hal){
        int nrows = var.rows(), ncols = var.cols();
        Grid dvar;
        dvar.resize(nrows, ncols);
        dvar.row(0) = var.row(1) - 2 * var.row(0) + hal.left;
        for (size_t i = 1; i < nrows - 1; i++)
            dvar.row(i) = DIFFN(2*MIN(i, nrows-1-i, order/2), var.row, i);
        dvar.row(nrows - 1) = var.row(nrows - 2) - 2 * var.row(nrows - 1) + hal.right;
        return dvar;
    }
    Grid y(const Grid &var, const Halo &hal){
        int nrows = var.rows(), ncols = var.cols();
        Grid dvar;
        dvar.resize(nrows, ncols);
        dvar.col(0) = var.col(1) - 2 * var.col(0) + hal.bottom;
        for (size_t i = 1; i < ncols - 1; i++)
            dvar.col(i) = DIFFN(2*MIN(i, ncols-1-i, order/2), var.col, i);
        dvar.col(ncols - 1) = var.col(ncols - 2) - 2 * var.col(ncols - 1) + hal.top;
        return dvar;
    }
};

class Integral{
    /* integrate over an axis, problem exist when you differential it on one axis
     * and integrate it over another axis */
protected:
    int order;
public:
    Integral(int _order = 0){}
    Grid operator() (const Grid &var, const Halo &hal, int axis){
        int nrows = var.rows(), ncols = var.cols();
        Grid dvar;
        if (axis == 2){
            dvar.resize(nrows, ncols - 1);
            dvar.col(0) = hal.bottom + var.col(0);
            for (size_t i = 1; i < ncols - 1; i++)
                dvar.col(i) = dvar.col(i - 1) + var.col(i);
        } else if (axis == -2){
            dvar.resize(nrows, ncols - 1);
            dvar.col(ncols - 2) = hal.top - var.col(ncols - 1);
            for (int i = ncols - 2; i > 0; i--)
                dvar.col(i - 1) = dvar.col(i) - var.col(i);
        } else if (axis == -1){
            dvar.resize(nrows - 1, ncols);
            dvar.row(nrows - 2) = hal.right - var.row(nrows - 1);
            for (int i = nrows - 2; i > 0; i--)
                dvar.row(i - 1) = dvar.row(i) - var.row(i);
        } else if (axis == 1){
            dvar.resize(nrows - 1, ncols);
            dvar.row(0) = hal.left + var.row(0);
            for (size_t i = 1; i < nrows - 1; i++)
                dvar.row(i) = dvar.row(i - 1) + var.row(i);
        } else {ASSERT_VARIABLE_OUT_OF_RANGE("axis");}   
        return dvar;
    }
};

class Interpolate{
    /* Make an interpolation to half grid */
protected:
    int order;
    Grid wsignx, wsigny, buffer;
public:
    Interpolate(){}
    Interpolate(int _order) : order(_order){
        if (_order != 2 && _order != 4 && _order != 6){
            ASSERT_VARIABLE_OUT_OF_RANGE("order");
        }
    }
    Grid operator() (const Grid &var, int axis){
        int nrows = var.rows(), ncols = var.cols();
        Grid dvar;
        if (axis == 1){
            dvar.resize(nrows - 1, ncols);
            for (size_t i = 0; i < nrows - 1; i++)
                dvar.row(i) = MIDDLE(2*MIN(i+1,nrows-1-i,order/2), var.row, i);
        } else if (axis == 2){
            dvar.resize(nrows, ncols - 1);
            for (size_t i = 0; i < ncols - 1; i++)
                dvar.col(i) = MIDDLE(2*MIN(i+1,ncols-1-i,order/2), var.col, i);
        } else {ASSERT_VARIABLE_OUT_OF_RANGE("axis");}
        return dvar;
    }
    Grid operator() (const Grid &var, const Halo &hal, int axis){
        int nrows = var.rows(), ncols = var.cols();
        Grid dvar;
        if (axis == 1){
            dvar.resize(nrows + 1, ncols);
            dvar.block(1, 0, nrows - 1, ncols) = (*this)(var, axis);
            dvar.row(0) = (var.row(0) + hal.left) / 2.;
            dvar.row(nrows) = (var.row(nrows - 1) + hal.right) / 2.;
        } else if (axis == 2){
            dvar.resize(nrows, ncols + 1);
            dvar.block(0, 1, nrows, ncols - 1) = (*this)(var, axis);
            dvar.col(0) = (var.col(0) + hal.bottom) / 2.;
            dvar.col(ncols) = (var.col(ncols - 1) + hal.top) / 2.;
        } else if (axis == 0){
            buffer.resize(nrows + 2, ncols + 2);
            buffer.block(1, 1, nrows, ncols) = var;
            buffer.block(0, 1, 1, ncols) = hal.left;
            buffer.block(nrows + 1, 1, 1, ncols) = hal.right;
            buffer.block(1, 0, nrows, 1) = hal.bottom;
            buffer.block(1, ncols + 1, nrows, 1) = hal.top;
            buffer(0, 0) = (buffer(0, 1) + buffer(1, 0)) / 2.;
            buffer(0, ncols + 1) = (buffer(0, ncols) + buffer(1, ncols + 1)) / 2.;
            buffer(nrows + 1, 0) = (buffer(nrows, 0) + buffer(nrows + 1, 1)) / 2.;
            buffer(nrows + 1, ncols + 1) = (buffer(nrows + 1, ncols) + buffer(nrows, ncols + 1)) / 2.;

            dvar = (*this)((*this)(buffer, 1), 2);
        }
        else {ASSERT_VARIABLE_OUT_OF_RANGE("axis");}
        return dvar;
    }
    Grid upwind(const Grid &wind, const Grid &var, int axis){
        int nrows = var.rows(), ncols = var.cols();
        Grid dvar;
        if (axis == 1){
            dvar.resize(nrows - 1, ncols);
            wsignx = (wind.block(1,0,nrows - 1,ncols) > 0).select(
                    ZERO2(nrows - 1, ncols) + 1., 
                    ZERO2(nrows - 1, ncols) - 1.);
            for (size_t i = 0; i < nrows - 1; i++)
                dvar.row(i) = MIDDLE(2*MIN(i+1, nrows-1-i, order/2), var.row, i) +
                    wsignx.row(i) * MIDDLE(2*MIN(i+1, nrows-1-i, order/2) - 1, var.row, i);
        } else if (axis == 2){
            dvar.resize(nrows, ncols - 1);
            wsigny = (wind.block(0,1,nrows,ncols - 1) > 0).select(
                    ZERO2(nrows, ncols - 1) + 1., 
                    ZERO2(nrows, ncols - 1) - 1.);
            for (size_t i = 0; i < ncols - 1; i++)
                dvar.col(i) = MIDDLE(2*MIN(i+1, ncols-1-i, order/2), var.col, i) +
                    wsigny.col(i) * MIDDLE(2*MIN(i+1, ncols-1-i, order/2) - 1, var.col, i);
        } else {ASSERT_VARIABLE_OUT_OF_RANGE("axis");}
        return dvar;
    }
    Grid upwind(const Grid &wind, const Grid &var, const Halo &hal, int axis){
        int nrows = var.rows(), ncols = var.cols();
        Grid dvar;
        if (axis == 1){
            dvar.resize(nrows + 1, ncols);
            dvar.block(1, 0, nrows - 1, ncols) = this->upwind(wind, var, axis);
            dvar.row(0) = (var.row(0) + hal.left) / 2. +
                (wind.row(0) > 0).select(ZERO2(1, ncols) + 1., ZERO2(1, ncols) - 1.)
                * 0.5 * (hal.left - var.row(0));
            dvar.row(nrows) = (var.row(nrows - 1) + hal.right) / 2. +
                (wind.row(nrows) > 0).select(ZERO2(1, ncols) + 1., ZERO2(1, ncols) - 1.)
                * 0.5 * (var.row(nrows - 1) - hal.right);
        } else if (axis == 2){
            dvar.resize(nrows, ncols + 1);
            dvar.block(0, 1, nrows, ncols - 1) = this->upwind(wind, var, axis);
            dvar.col(0) = (var.col(0) + hal.bottom) / 2. + 
                (wind.col(0) > 0).select(ZERO2(nrows, 1) + 1., ZERO2(nrows, 1) - 1.) 
                * 0.5 * (hal.bottom - var.col(0));
            dvar.col(ncols) = (var.col(ncols - 1) + hal.top) / 2. + 
                (wind.col(ncols) > 0).select(ZERO2(nrows, 1) + 1., ZERO2(nrows, 1) - 1.) 
                * 0.5 * (var.col(ncols - 1) - hal.top);
        } else {ASSERT_VARIABLE_OUT_OF_RANGE("axis");}
        return dvar;
    }
};

#endif
