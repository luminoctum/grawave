#! /usr/bin/env python2.7
from ShallowWater import *
from moist import *
from pycli.ode.gridop import *
from scipy.integrate import trapz
from pycli.data.SaturnData import *
from scipy.interpolate import UnivariateSpline as spline

def heating(xH2O = 0, np = 0, T1bar = 0, return_diag = 0):

    data = SaturnData()
    data.import_rss(pmin = 10.E3, pmax = 30.E3)
    atmos = Atmosphere(P0 = 1.E5, T0 = 134.8, grav = 10.44, xHe = 0.034)
    # First diag, used to generate mean temperature profile connects to the Voyger observation
    diag = ThermoDiagram(pmin = 30.E3, pmax = 40.E5, np = 200, atmos = atmos)
    H2O = Water()
    array1 = hstack([log(diag.ptol), log(data.rss['lindal85'].pres)])[::-1]
    array2 = hstack([diag._moist_adiabats(atmos.T0, H2O), data.rss['lindal85'].temp])[::-1]

    # Second diag, used to determine cloud level
    diag = ThermoDiagram(pmin = 1.E5, pmax = 40.E5, np = 200, atmos = atmos)
    diag._adiabats(T1bar, H2O, xH2O)
    pcloud = diag.cloudb_ptol

    # Third diag, used to determine the axes of the domain
    diag = ThermoDiagram(pmin = 10.E3, pmax = pcloud, np = np, atmos = atmos)

    temp_high = diag._adiabats(T1bar, H2O, xH2O)
    temp_bar = amap(lambda p: spline(array1, array2, k = 2, s = 0)(p), log(diag.ptol))
    theta_bar = amap(atmos._theta, temp_bar, diag.ptol)

    N2_bar = atmos.grav * fdiff(theta_bar, diag.zlev)
    N2_bar = hstack([N2_bar[0], N2_bar, N2_bar[-1]])
    N2_bar[N2_bar < 0.] = 0.
    for i in range(2, len(N2_bar) - 2):
        if N2_bar[i-1] == 0 and N2_bar[i+1] == 0: N2_bar[i] = 0.
        if N2_bar[i-2] == 0 and N2_bar[i+2] == 0: N2_bar[i] = 0.

    svp = H2O._sat_vapor_pressure(temp_bar)
    eta_bar = svp / (diag.ptol - svp)

    buoyancy = atmos.grav * (temp_high - temp_bar).clip(0, 1E6) / temp_bar
    if return_diag == 1: return diag 
    else: return temp_bar, N2_bar, eta_bar, buoyancy

class ColumnHeating(ShallowWater):
    id = '2'

    def set_variables(self):
        atmos = Atmosphere(P0 = 1.E5, T0 = 134.8, grav = 10.44, xHe = 0.034)
        H2O = Water()
        eps = H2O.mu / atmos.mu
        T_, N2_, eta_, b = heating(xH2O = self.xH2O, np = self.ny, T1bar = self.temph)

        b = array([b * exp(-x**2 / (2. * self.sigma**2)) for x in self.xaxis])
        T_ = array([T_ for x in self.xaxis])
        eta_ = array([eta_ for x in self.xaxis])
        phi_ = zeros((self.nx, self.ny + 1))
        for i in range(1, self.ny + 1):
            phi_[:, i] = phi_[:, i - 1] \
                    + (1 + eta_[:, i - 1]) / (1 + eps * eta_[:, i - 1]) * atmos.grav * T_[:, i - 1] / atmos.T0 \
                    * (self.yaxisb[i] - self.yaxisb[i - 1])
        svp = H2O._sat_vapor_pressure(T_ * (1 + b / atmos.grav))
        ptol = array([atmos.P0 * exp(- self.yaxis / atmos.H0) for x in self.xaxis])
        eta = svp / (ptol - svp)
        eta = eta.clip(0, self.xH2O)
        RH = eta / (1 + eta) * ptol / svp
        N2 = array([N2_ for x in self.xaxis])
        dphidz = (T_ / atmos.T0) * b * (1 + eta) / (1 + eps * eta) + \
                atmos.grav * T_ / atmos.T0 * (1 - eps) / (1 + eps * eta) * (eta - eta_) / ( 1 + eps * eta_)
        phi = zeros((self.nx, self.ny + 1))
        for i in range(1, self.ny + 1):
            phi[:, i] = phi[:, i - 1] + dphidz[:, i - 1] * (self.yaxisb[i] - self.yaxisb[i - 1])
        x0 = self.xlen * (1 - 0.4)
        absorb = zeros((self.nx, self.ny))
        for i, x in enumerate(self.xaxis):
            if x <= x0:
                absorb[i, :] = 0
            else:
                absorb[i, :] = self.gamma * (exp(((x - x0) / (self.xlen - x0))**4) - 1)
        absorbx = tohalf(absorb, axis = 0, ext = 'both')
        self.var['west_east'] = self.xaxis
        self.var['west_eastb'] = self.xaxisb
        self.var['zplev'] = self.yaxis
        self.var['zplevb'] = self.yaxisb
        self.var['buoyancy'] = b
        self.var['phi'] = phi
        self.var['mass'] = array([exp(- self.yaxis / atmos.H0) for x in self.xaxis])
        self.var['massx'] = array([exp(- self.yaxis / atmos.H0) for x in self.xaxisb])
        self.var['massy'] = array([exp(- self.yaxisb / atmos.H0) for x in self.xaxis])
        self.var['N2'] = N2
        self.var['T_'] = T_
        self.var['eta_'] = eta_
        self.var['phi_'] = phi_
        self.var['eta_H2O'] = eta
        #self.var['RH_H2O'] = 1. + zeros((self.nx, self.ny))
        self.var['RH_H2O'] = RH
        self.var['absorb'] = absorb
        self.var['absorbx'] = absorbx

    def write_parameters(self, f0 = 1.E04, kxx = 1.E-2, kzz = 1.E-2):
        ShallowWater.write_parameters(self, f0 = f0)
        file = open(self.ctrl_file, 'a')
        file.write('%-20s%-8.3E\n' % ('"kxx"', kxx))
        file.write('%-20s%-8.3E\n' % ('"kzz"', kzz))
        file.close()

if __name__ == '__main__':
    name    =   'column heating'
    nx      =   64
    ny      =   64
    xlen    =   5.E6
    tstart  =   0.
    tend    =   300000.
    tstep   =   5.0
    frame   =   20
    f0      =   2.128E-4
    beta    =   0.E-11
    kxx     =   1.E+0
    kzz     =   1.E+0
    sigma   =   2.E5
    gamma   =   0.04
    xH2O    =   0.012

    ax = draw_prepare()
    T2 = fsolve(lambda x: compute(ax, x, xH2O) - 134.8, 150)
    clf()
    #T2 = 144.86

    diag = heating(xH2O = xH2O, np = ny, T1bar = T2, return_diag = 1)
    xaxis = linspace(xlen, 0, nx, endpoint = False)
    xaxis = xaxis[::-1]
    yaxis = diag.zlev
    ylen  = yaxis[-1] - yaxis[0]
    model = ColumnHeating(name, nx, ny, xlen, ylen, tstart, tend, tstep, frame)
    model.write_parameters(f0 = f0, kxx = kxx, kzz = kzz)
    model.xaxis = xaxis
    model.xaxisb= tohalf(model.xaxis, ext = 'both')
    model.yaxis = yaxis
    model.yaxisb= tohalf(model.yaxis, ext = 'both')
    model.sigma = sigma
    model.gamma = gamma
    model.xH2O  = xH2O
    model.temph = T2
    model.set_variables()
    model.write_ncfile()
