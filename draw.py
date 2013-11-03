#! /usr/bin/env python2.7
from pylab import *
from pycli.ncio.ncfile import *
from pycli.other.tools import *

class struct: pass

infile  = open('control.in','r')
while not ('Time and domain' in infile.readline()): pass
nx      = int(infile.readline().split()[-1])
ny      = int(infile.readline().split()[-1])
xlen    = float(infile.readline().split()[-1])
ylen    = float(infile.readline().split()[-1])
tbegin  = float(infile.readline().split()[-1])
tend    = float(infile.readline().split()[-1])
step    = float(infile.readline().split()[-1])
frame   = int(infile.readline().split()[-1])
while not ('Parameters' in infile.readline()): pass
line = infile.readline().split()
params = {}
while len(line) > 1:
    id = line[0][1:-1]
    value = float(line[1])
    params[id] = value
    line = infile.readline().split()
infile.close()

data = Dataset('dynamics.nc', 'r')
xaxis = data.variables['west_east']
yaxis = data.variables['zplev']
times = data.variables['time']
X, Y = meshgrid(xaxis, yaxis)

phi = zeros((len(times), len(xaxis), len(yaxis)))
for i in range(len(times)):
    temp = data.variables['phi_'][:, :] + data.variables['phi'][i, :, :]
    phi[i, :, :] = temp.T
print phi.shape
figure()
ax = axes()
ylim([min(phi[-1,:,12]), max(phi[0,:,12])])
for i in range(0, len(times)/2, 80):
    ax.plot(xaxis, phi[i, :, 12])
show()
