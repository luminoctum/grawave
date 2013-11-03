#! /usr/bin/env python2.7
from OdeSystem import *
from pycli.ode.gridop import *

class ShallowWater(OdeSystem):
    id = '1'

    def set_variables(self):
        sigma   = 1.E6
        x0, y0 = self.xlen * 4/5., self.ylen/2.
        X, Y = meshgrid(self.xaxis, self.yaxis)
        phi = 10000 * (1 + 2 * \
                exp(-((X - x0)**2 + (Y - y0)**2)/(2. * sigma * sigma)))
        self.var['phi'] = phi.T;
        self.var['west_east'] = self.xaxis
        self.var['south_north'] = self.yaxis
        self.var['west_eastb'] = tohalf(self.xaxis, ext = 'both')
        self.var['south_northb'] = tohalf(self.yaxis, ext = 'both')
        self.var['tracer'] = X.T
    
    def write_parameters(self, f0 = 1.E-4, beta = 1.E-11):
        file = open(self.ctrl_file, 'a')
        file.write('// Parameters\n')
        file.write('%-20s%-8.3E\n' % ('"f0"', f0))
        file.write('%-20s%-8.3E\n' % ('"beta"', beta))
        file.close()

if __name__ == '__main__':
    name    =   'shallow water system'
    nx      =   400
    ny      =   100
    xlen    =   40.E6
    ylen    =   10.E6
    tstart  =   0.
    tend    =   180000.
    tstep   =   600.
    frame   =   1
    f0      =   1.E-4
    beta    =   1.E-11

    model = ShallowWater(name, nx, ny, xlen, ylen, tstart, tend, tstep, frame)
    model.write_parameters(f0, beta)
    model.set_variables()
    model.write_ncfile()
