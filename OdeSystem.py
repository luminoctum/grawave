#! /usr/bin/env python2.7
from pylab import *
from pycli.ncio.ncfile import *

class OdeSystem:
    ctrl_file   = 'control.in'
    var_file    = 'variables.lst'
    nc_file     = 'dynamics.nc'
    id          = '0'
    var         = {}
    def __init__(self, 
            name = '',
            nx = 0,
            ny = 0,
            xlen = 0,
            ylen = 0,
            tstart = 0.,
            tend = 0.,
            tstep = 0.,
            frame = 0):
        self.name = name
        self.nx = nx
        self.ny = ny
        self.xlen = xlen
        self.ylen = ylen
        self.tstart = tstart
        self.tend = tend
        self.tstep = tstep
        self.frame = frame
        self.xaxis = linspace(0, self.xlen, self.nx)
        self.yaxis = linspace(0, self.ylen, self.ny)

        self.write_domain()

    def set_variables(): pass

    def write_domain(self):
        file = open(self.ctrl_file, 'w')
        file.write('Experiment: ' + self.name + '\n\n')
        file.write('// Time and domain\n')
        file.write('%-20s%-8i\n' % ('num_of_grids_in_x', self.nx))
        file.write('%-20s%-8i\n' % ('num_of_grids_in_y', self.ny))
        file.write('%-20s%-8.3E\n' % ('x_length', self.xlen))
        file.write('%-20s%-8.3E\n' % ('y_length', self.ylen))
        file.write('%-20s%-8.3f\n' % ('time_start', self.tstart))
        file.write('%-20s%-8.3f\n' % ('time_end', self.tend))
        file.write('%-20s%-8.3f\n' % ('time_step', self.tstep))
        file.write('%-20s%-8i\n\n' % ('time_per_frame', self.frame))
        file.close()

    def write_ncfile(self):
        varlist = genfromtxt(self.var_file,
                dtype = 'string',
                delimiter = '"',
                skip_header = 1,
                usecols = [1, 3, 5, 7, 9])
        file = ncfile(self.nc_file, 'w')
        dims = {}
        for i in range(varlist.shape[0]):
            if self.id in varlist[i, 2]:
                if varlist[i, 1][0] == 'd':
                    if varlist[i, 0] == 'time':
                        file.add_dim(varlist[i, 0], [], varlist[i, 3], varlist[i, 4])
                    else:
                        file.add_dim(varlist[i, 0], self.var[varlist[i, 0]], varlist[i, 3], varlist[i, 4])
                    dims[varlist[i, 1][1]] = varlist[i, 0]
                else:
                    dim_name, dim_length = (), ()
                    for key in varlist[i, 1]: 
                        dim_name = dim_name + (dims[key],)
                        if dims[key] != 'time':
                            dim_length = dim_length + (len(self.var[dims[key]]),)
                    if varlist[i, 0] in self.var.keys():
                        file.add_var(varlist[i, 0], dim_name, self.var[varlist[i, 0]].T, varlist[i, 3], varlist[i, 4])
                    else:
                        file.add_var(varlist[i, 0], dim_name, zeros(dim_length), varlist[i, 3], varlist[i, 4])
        file.time[:] = 0
