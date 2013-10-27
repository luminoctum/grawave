#ifndef FINITEMethod
#define FINITEMethod

#include "Include.hh"
#define MAXDIM 500
#define Min(a, b, c) ( a < b ? (a < c ? a : c) : (b < c ? b : c) )
#define D_order_1(var, i) ( var(i + 1) - var(i) )
#define D_order_1_m(var, i) ( var(i) - var(i - 1) )
#define D_order_2(var, i) ( (var(i + 1) - var(i - 1)) / 2. )
#define D_bnd_low(var, varb) ( \
        varb.type == Neumann ? varb.value : \
            static_cast<Grid>(var - varb.value))
#define D_bnd_high(var, varb) ( \
        varb.type == Neumann ? varb.value :  \
            static_cast<Grid>(varb.value - var))
#define HalfPointAvg(order, var, i) ( \
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
#define HalfPointLow(var, varb) ( \
        varb.type == Neumann ? var - 0.5 * varb.value : \
            static_cast<Grid>((var + varb.value) / 2.))
#define HalfPointHigh(var, varb) ( \
        varb.type == Neumann ? var + 0.5 * varb.value : \
            static_cast<Grid>((var + varb.value) / 2.))
#define BndLow(var, varb) ( \
        varb.type == Neumann ? var - varb.value : varb.value )
#define BndHigh(var, varb) ( \
        varb.type == Neumann ? var + varb.value : varb.value )
#define Nabla(order, var, i) ( \
        order == 2 ? var(i-1)-2*var(i)+var(i+1) : ( \
        order == 4 ? \
            static_cast<Grid>(-var(i-2)+4*var(i-1)-6*var(i)+4*var(i+1)-var(i+2)) : \
        /* order == 6 */ \
            static_cast<Grid>(var(i-3)-6*var(i-2)+15*var(i-1)-20*var(i)+15*var(i+1)-6*var(i+2)+var(i+3))) \
        )

class FiniteMethod{
    /* axis = 2  : difference up
     * axis = 1  : difference right
     * axis = -2 : difference down
     * axis = -1 : difference left
     */
protected:
    Grid Buffer, Buffer2, Sum;
public:
    FiniteMethod(){}
    FiniteMethod(int _nrows, int _ncols) : 
        Buffer(_nrows, _ncols), Buffer2(_nrows, _ncols){}
    virtual Grid operator() (const Grid &, const Boundary &, int){};
    virtual Grid operator() (const Grid &, int){};
};

class CenteredDifference : public FiniteMethod{
#define BUFFER Buffer.block(0, 0, var.rows(), var.cols())
public:
    CenteredDifference(int _nrows = MAXDIM, int _ncols = MAXDIM) : FiniteMethod(_nrows, _ncols){}
    Grid operator() (const Grid &var, const Boundary &varb, int axis){
        if (axis == 1){
            BUFFER.row(0) = (D_bnd_low(var.row(0), varb.left) 
                    + D_order_1(var.row, 0)) / 2.;
            for (size_t i = 1; i < var.rows() - 1; i++)
                BUFFER.row(i) = D_order_2(var.row, i);
            BUFFER.row(var.rows() - 1) = (D_bnd_high(var.row(var.rows() - 1), varb.right) 
                    + D_order_1_m(var.row, var.rows() - 1)) / 2.;
        } else if (axis == 2){
            BUFFER.col(0) = (D_bnd_low(var.col(0), varb.bottom)
                    + D_order_1(var.col, 0)) / 2.;
            for (size_t i = 1; i < var.cols() - 1; i++)
                BUFFER.col(i) = D_order_2(var.col, i);
            BUFFER.col(var.cols() - 1) = (D_bnd_high(var.col(var.cols() - 1), varb.top)
                    + D_order_1_m(var.col, var.cols() - 1)) / 2.;
        } else {ASSERT_VARIABLE_OUT_OF_RANGE("axis");}
        return BUFFER;
    }
#undef BUFFER
};

class BiasedDifference : public FiniteMethod{
public:
    BiasedDifference(int _nrows = MAXDIM, int _ncols = MAXDIM) : FiniteMethod(_nrows, _ncols){}
    Grid operator() (const Grid &var, const Boundary &varb, int axis){
        int nrows = var.rows();
        int ncols = var.cols();
        if (axis == 1){
            Buffer.block(0,0,nrows+1,ncols).row(0) = D_bnd_low(var.row(0), varb.left);
            for (size_t i = 0; i < nrows - 1; i++)
                Buffer.block(0,0,nrows+1,ncols).row(i + 1) = D_order_1(var.row, i);
            Buffer.block(0,0,nrows+1,ncols).row(nrows) = D_bnd_high(var.row(nrows - 1), varb.right);
            return Buffer.block(0,0,nrows+1,ncols);
        } else if (axis == 2){
            Buffer.block(0,0,nrows,ncols+1).col(0) = D_bnd_low(var.col(0), varb.bottom);
            for (size_t i = 0; i < ncols - 1; i++)
                Buffer.block(0,0,nrows,ncols+1).col(i + 1) = D_order_1(var.col, i);
            Buffer.block(0,0,nrows,ncols+1).col(ncols) = D_bnd_high(var.col(ncols - 1), varb.top);
            return Buffer.block(0,0,nrows,ncols+1);
        } else{ASSERT_VARIABLE_OUT_OF_RANGE("axis");}
    }
    Grid operator() (const Grid &var, int axis){
        int nrows = var.rows();
        int ncols = var.cols();
        if (axis == 1){
            for (size_t i = 0; i < nrows - 1; i++)
                Buffer.block(0,0,nrows-1,ncols).row(i) = D_order_1(var.row, i);
            return Buffer.block(0,0,nrows-1,ncols);
        } else if (axis == 2){
            for (size_t i = 0; i < ncols - 1; i++)
                Buffer.block(0,0,nrows,ncols-1).col(i) = D_order_1(var.col, i);
            return Buffer.block(0,0,nrows,ncols-1);
        } else{ASSERT_VARIABLE_OUT_OF_RANGE("axis");}
    }
};

class Diffusion: public FiniteMethod{
#define BUFFER Buffer.block(0,0,nrows,ncols)
protected:
    int order; // has to be even
public:
    Diffusion(int _nrows = MAXDIM, int _ncols = MAXDIM, int _order = 4) 
        : FiniteMethod(_nrows, _ncols), order(_order){}
    Grid operator() (const Grid &var, int axis){
        int nrows = var.rows();
        int ncols = var.cols();
        if (axis == 1){
            BUFFER.row(0) = Grid::Zero(1, ncols);
            for (size_t i = 1; i < nrows - 1; i++)
                BUFFER.row(i) = Nabla(2*Min(i, nrows-1-i, order/2), var.row, i);
            BUFFER.row(nrows - 1) = Grid::Zero(1, ncols);
        } else if (axis == 2){
            BUFFER.col(0) = Grid::Zero(nrows, 1);
            for (size_t i = 1; i < ncols - 1; i++)
                BUFFER.col(i) = Nabla(2*Min(i, ncols-1-i, order/2), var.col, i);
            BUFFER.col(ncols - 1) = Grid::Zero(nrows, 1);
        } else{ASSERT_VARIABLE_OUT_OF_RANGE("axis");}
        return BUFFER;
    }
#undef BUFFER
};

class Integral : public FiniteMethod{
#define BUFFER Buffer.block(0,0,nrows,ncols+1)
public:
    Integral(int _nrows = MAXDIM, int _ncols = MAXDIM) : FiniteMethod(_nrows, _ncols){};
    Grid operator() (const Grid &var, const Boundary &varb, int axis){
        int nrows = var.rows();
        int ncols = var.cols();
        if (axis == 2){
            BUFFER.col(0) = varb.bottom.value;
            for (size_t i = 0; i < ncols; i++)
                BUFFER.col(i + 1) = BUFFER.col(i) + var.col(i);
        } else if (axis == -2){
            BUFFER.col(ncols) = varb.top.value;
            for (int i = ncols; i > 0; i--)
                BUFFER.col(i - 1) = BUFFER.col(i) + var.col(i - 1);
        } else{ASSERT_VARIABLE_OUT_OF_RANGE("axis");}   
        return BUFFER;
    }
#undef BUFFER
};

class Interpolate: public FiniteMethod{
protected:
    int order; //has to be even number
    Grid ones;
public:
    Interpolate(int _nrows = MAXDIM, int _ncols = MAXDIM, int _order = 2) 
        : FiniteMethod(_nrows, _ncols), order(_order){
            ones = Grid::Zero(_nrows, _ncols) + 1.;
    }
    Grid operator() (const Grid &var, int axis){
        int nrows = var.rows();
        int ncols = var.cols();
        if (axis == 1){
            for (size_t i = 0; i < nrows - 1; i++)
                Buffer.block(i,0,1,ncols) = 
                    HalfPointAvg(2*Min(i+1, nrows-1-i, order/2), var.row, i);
            return Buffer.block(0,0,nrows-1,ncols);
        } else if (axis == 2){
            for (size_t i = 0; i < ncols - 1; i++)
                Buffer.block(0,i,nrows,1) = 
                    HalfPointAvg(2*Min(i+1, ncols-1-i, order/2), var.col, i);
            return Buffer.block(0,0,nrows,ncols-1);
        } else{ASSERT_VARIABLE_OUT_OF_RANGE("axis");}
    }
    Grid operator() (const Grid &var, const Boundary &varb, int axis){
        int nrows = var.rows();
        int ncols = var.cols();
        if (axis == 1){
            Buffer.block(1,0,nrows-1,ncols) = (*this)(var, axis);
            Buffer.block(0,0,1,ncols) = HalfPointLow(var.row(0), varb.left);
            Buffer.block(nrows,0,1,ncols) = HalfPointHigh(var.row(nrows - 1), varb.right);
            return Buffer.block(0,0,nrows+1,ncols);
        } else if (axis == 2){
            Buffer.block(0,1,nrows,ncols-1) = (*this)(var, axis);
            Buffer.block(0,0,nrows,1) = HalfPointLow(var.col(0), varb.bottom);
            Buffer.block(0,ncols,nrows,1) = HalfPointHigh(var.col(ncols - 1), varb.top);
            return Buffer.block(0,0,nrows,ncols+1);
        } else{ASSERT_VARIABLE_OUT_OF_RANGE("axis");}
    }
    Grid operator() (const Grid &wind, const Grid &var, int axis){
        int nrows = var.rows();
        int ncols = var.cols();
        if (axis == 1){
            Buffer2 = wind.block(1,0,nrows-1,ncols);
            Buffer2 = (Buffer2 > 0).select(ones.block(0,0,nrows-1,ncols), -ones.block(0,0,nrows-1,ncols));
            for (size_t i = 0; i < nrows - 1; i++)
                Buffer.block(i,0,1,ncols) = 
                    HalfPointAvg(2*Min(i+1, nrows-1-i, order/2), var.row, i) +
                    Buffer2.row(i) * HalfPointAvg(2*Min(i+1, nrows-1-i, order/2) - 1, var.row, i);
            return Buffer.block(0,0,nrows-1,ncols);
        } else if (axis == 2){
            Buffer2 = wind.block(0,1,nrows,ncols-1);
            Buffer2 = (Buffer2 > 0).select(ones.block(0,0,nrows,ncols-1), -ones.block(0,0,nrows,ncols-1));
            for (size_t i = 0; i < ncols - 1; i++)
                Buffer.block(0,i,nrows,1) = 
                    HalfPointAvg(2*Min(i+1, ncols-1-i, order/2), var.col, i) +
                    Buffer2.col(i) * HalfPointAvg(2*Min(i+1, ncols-1-i, order/2) - 1, var.col, i);
            return Buffer.block(0,0,nrows,ncols-1);
        } else{ASSERT_VARIABLE_OUT_OF_RANGE("axis");}
    }
    Grid operator() (const Grid &wind, const Grid &var, const Boundary &varb, int axis){
        int nrows = var.rows();
        int ncols = var.cols();
        if (axis == 1){
            Buffer.block(1,0,nrows-1,ncols) = (*this)(wind, var, axis);
            Buffer.block(0,0,1,ncols) = 
                HalfPointLow(var.row(0), varb.left) +
                (wind.row(0) > 0).select(ones.block(0,0,1,ncols), - ones.block(0,0,1,ncols)) 
                * 0.5 * (BndLow(var.row(0), varb.left) - var.row(0));
            Buffer.block(nrows,0,1,ncols)  = 
                HalfPointHigh(var.row(nrows-1), varb.right) +
                (wind.row(nrows) > 0).select(ones.block(0,0,1,ncols), - ones.block(0,0,1,ncols)) 
                * 0.5 * (var.row(nrows-1) - BndHigh(var.row(nrows-1), varb.right));
            return Buffer.block(0,0,nrows+1,ncols);
        } else if (axis == 2){
            Buffer.block(0,1,nrows,ncols-1) = (*this)(wind, var, axis);
            Buffer.block(0,0,nrows,1) = 
                HalfPointLow(var.col(0), varb.bottom) +
                (wind.col(0) > 0).select(ones.block(0,0,nrows,1), - ones.block(0,0,nrows,1)) 
                * 0.5 * (BndLow(var.col(0), varb.bottom) - var.col(0));
            Buffer.block(0,ncols,nrows,1) = 
                HalfPointHigh(var.col(ncols-1), varb.top) +
                (wind.col(ncols) > 0).select(ones.block(0,0,nrows,1), - ones.block(0,0,nrows,1)) 
                * 0.5 * (var.col(ncols-1) - BndHigh(var.col(ncols-1), varb.top));
            return Buffer.block(0,0,nrows,ncols+1);
        } else{ASSERT_VARIABLE_OUT_OF_RANGE("axis");}
    }
};

Grid fluxdiv(FiniteMethod &diff, FiniteMethod &interp, const Grid &wind, 
        const Grid &var, const Boundary &varb, int axis){
    return diff(wind * interp(var, varb, axis), axis);
};
#endif
