#! /usr/bin/env python2.7
from pylab import *
from pycli.ode.gridop import *
from pycli.data.SaturnData import *
from pycli.other.mouse import *
from pycli.other.tools import *
from scipy.optimize import fsolve
from scipy.interpolate import UnivariateSpline as spline

class Species:
    def _sat_mixr(self, temp, ptol):
        svp = self._sat_vapor_pressure(temp)
        return svp / (ptol - svp)

class Water(Species):
    condensible = True
    mu = 18.E-3
    R = 461.67 # 8.31 / mu
    cp = 2077.5 # 4.5 * R
    cl = 4.18E3
    Lv = 2.6E6 # or 2.26E6
    def _sat_vapor_pressure(self, temp):
        return 611.2 * exp(17.67 * (temp - 273.15) / (temp - 29.65))
    def _temp_from_sat_mixr(self, smr, ptol):
        svp = ptol * smr / (1 + smr)
        return 243.5 / ((17.67 / log(svp / 611.2)) - 1) + 273.15

class Helium:
    condensible = False
    mu = 4.E-3
    R = 2077.5 # 8.31 / mu
    cp = 5193.75 # 2.5 * R

class Hydrogen:
    condensible = False
    mu = 2.E-3
    R = 4155 # 8.31 / mu
    cp = 14542.5 # 3.5 * R

class Atmosphere:
    def __init__(self, 
            P0 = 0, 
            T0 = 0, 
            grav = 0, 
            xHe = 0):
        H2, He, H2O = Hydrogen(), Helium(), Water()
        self.xHe = xHe
        self.mu = (H2.mu + He.mu * xHe) / (1. + xHe)
        self.R = 8.31 / self.mu
        # scale cp by 0.9 to match voyer observation
        self.cp = 0.9 * (H2.cp * H2.mu + H2.cp * He.mu * xHe) / (H2.mu + He.mu * xHe)
        self.P0, self.T0, self.grav = P0, T0, grav
        self.H0 = self.R * T0 / grav
    def _theta(self, temp, ptol):
        return temp * (self.P0 / ptol) ** (self.R / self.cp)
    def _sat_theta_e(self, temp, ptol, speci): # approximated formula
        smr = speci._sat_mixr(temp, ptol)
        smr = smr / (1 + smr)
        return self._theta(temp, ptol - speci._sat_vapor_pressure(temp)) \
                * exp(speci.Lv * speci.mu / self.mu * smr / (self.cp * temp))
    def _temp_v(self, temp, mixr, speci):
        return temp * (1 + mixr) / (1 + speci.mu / self.mu * mixr)
    def _sat_temp_v(self, temp, ptol, speci):
        smr = speci._sat_mixr(temp, ptol)
        return temp * (1 + smr) / ( 1 + speci.mu / self.mu * smr)
    def _temp_from_theta(self, theta, ptol):
        guess = theta * (ptol / self.P0) ** (self.R / self.cp)
        def zero(p, t): return lambda x: self._theta(x, p) - t
        return fsolve(zero(ptol, theta), guess, xtol = 1.E-5)[0]
    def _temp_from_theta_e(self, theta_e, ptol, speci):
        guess = min(400, theta_e * (ptol / self.P0) ** (self.R / self.cp))
        def zero(p, t): return lambda x: self._sat_theta_e(x, p, speci) - t
        return fsolve(zero(ptol, theta_e), guess, xtol = 1.E-5)[0]

class ThermoDiagram:
    def __init__(self, 
            pmin = 0, pmax = 0, np = 0,
            tmin = 0, tmax = 0, nt = 0, 
            atmos = 0):
        self.pmin, self.pmax, self.np = pmin, pmax, np
        self.tmin, self.tmax, self.nt = tmin, tmax, nt
        self.ptol = logspace(log10(self.pmax), log10(self.pmin), self.np)
        self.temp = linspace(self.tmin, self.tmax, self.nt)
        self.zlev = - atmos.H0 * log (self.ptol / atmos.P0)
        self.zlevh = tohalf(self.zlev)
        self.ptolh = atmos.P0 * exp ( - self.zlevh / atmos.H0)
        self.atmos = atmos
    def _adiabats(self, temp0, speci, xspeci):
        ltemp, ptol = [], self.ptol[::-1]
        for i, p in enumerate(ptol):
            temp = self.atmos._temp_from_theta_e(temp0, p, speci)
            svp = speci._sat_vapor_pressure(temp)
            mixr = svp / (p - svp)
            if mixr < xspeci: 
                ltemp.append(temp)
            else: 
                ltemp.append(temp)
                theta_cloud = self.atmos._theta(temp, p)
                self.cloudb_ptol = p
                self.cloudb_temp = temp
                self.cloudb_temp_v = self.atmos._temp_v(temp, mixr, speci)
                break
        if i < self.np - 1:
            for j, p in enumerate(ptol[i+1:]):
                temp = self.atmos._temp_from_theta(theta_cloud, p)
                ltemp.append(temp)
        return array(ltemp[::-1])
    def _moist_adiabats(self, temp0, speci):
        return amap(lambda p: self.atmos._temp_from_theta_e(temp0, p, speci), self.ptol)
    def _dry_adiabats(self, temp0):
        return amap(lambda p: self.atmos._temp_from_theta(temp0, p), self.ptol)
    def _sat_mixr(self, temp0, speci):
        temp = self._moist_adiabats(temp0, speci)
        svp = speci._sat_vapor_pressure(temp)
        smr = svp / (self.ptol - svp)
        return smr

def compute(ax = 0, T1bar = 0, xH2O = 0):
    if ax == 0:
        figure(100)
        ax = axes()
        flag = 1
    else: flag = 0
    data = SaturnData()
    data.import_rss()
    atmos = Atmosphere(P0 = 1.E5, T0 = T1bar, grav = 10.44, xHe = 0.034)
    diag = ThermoDiagram(pmin = 10E5, pmax = 40.E5, np = 200,
            tmin = 260., tmax = 360., nt = 100, atmos = atmos)
    H2O = Water()

    # moist & dry adiabats
    ax.plot(diag._adiabats(atmos.T0, H2O, xH2O), diag.ptol, 'b-', linewidth = 2)
    # virtual temperature
    ax.plot(amap(lambda t, x: atmos._temp_v(t, x, H2O), \
            diag._adiabats(atmos.T0, H2O, xH2O), \
            diag._sat_mixr(atmos.T0, H2O)), 
            diag.ptol, 'r-', linewidth = 2)
    # contours of virtual temperature 
    X, Y = meshgrid(diag.temp, diag.ptol)
    h = ax.contour(X, Y, 
            map(lambda x,y: atmos._sat_temp_v(x,y,H2O), X, Y), 
            [diag.cloudb_temp_v],
            colors = 'c', linewidths = 2)
    cvt = h.collections[0].get_paths()[0].vertices # vertices of const virtual temperature
    fitted = spline(cvt[:, 0], cvt[:, 1], k = 3, s= 0)
    # voyager observation
    ax.plot(data.rss['lindal85'].temp, data.rss['lindal85'].pres, 'k', linewidth = 4)
    # cloud level
    ax.plot([diag.tmin, diag.tmax], [diag.cloudb_ptol, diag.cloudb_ptol], 'g-', linewidth = 3)
    # const mixing ratio
    ax.plot(amap(lambda p: H2O._temp_from_sat_mixr(xH2O, p), diag.ptol), 
            diag.ptol, 'm-', linewidth = 2)

    cloudb_ptol = diag.cloudb_ptol
    cloubb_temp = diag.cloudb_temp
    cloudb_temp2 = fsolve(lambda x: fitted(x) - diag.cloudb_ptol, 
            diag.cloudb_temp - 20, xtol = 1.E-5)[0]
    atmos.T1 = fsolve(lambda x: 
            atmos._temp_from_theta_e(x, diag.cloudb_ptol, H2O) - cloudb_temp2, 
            atmos.T0 - 10, xtol = 1.E-5)[0]
    ax.plot(diag._adiabats(atmos.T1, H2O, xH2O)[diag.ptol < cloudb_ptol], 
            diag.ptol[diag.ptol < cloudb_ptol], 'b--', linewidth = 2)
    ax.plot(amap(lambda x, y: atmos._temp_v(x, y, H2O), \
            diag._adiabats(atmos.T1, H2O, xH2O)[diag.ptol < cloudb_ptol], \
            diag._sat_mixr(atmos.T1, H2O)[diag.ptol < cloudb_ptol]), 
            diag.ptol[diag.ptol < cloudb_ptol], 'r--', linewidth = 2)
    pval = array(range(1, 10) + range(10, 20, 2) + range(20, 50, 5))
    pval = hstack([arange(0.1, 1, 0.2), pval])
    tval = arange(100, 400, 10)
    ax.set_yscale('log')
    ax.set_yticks(pval * 1E5)
    ax.set_yticklabels(pval)
    ax.set_xticks(tval)
    ax.set_ylabel('Pressure (bar)', fontsize = 20, weight = 'bold')
    ax.set_xlabel('Temperature (K)', fontsize = 20, weight = 'bold')
    ax.set_xlim([diag.tmin, diag.tmax])
    ax.set_ylim([diag.pmax, diag.pmin])
    if flag: close()
    else:
        grid(True)
        ax2 = ax.twinx()
        ax2.set_ylim(min(diag.zlev / 1E3), max(diag.zlev / 1E3))
        ax2.set_ylabel('$\mathbf{Z^*}$ (km)', fontsize = 20, weight = 'bold')
    print 'T1 = %7.2f, T2 = %7.2f' % (atmos.T1, atmos.T0)
    return atmos.T1

def cape_cin_rad(ax = 0, xH2O = 0):
    atmos = Atmosphere(P0 = 1.E5, T0 = 134.8, grav = 10.44, xHe = 0.034)
    H2O = Water()
    T1 = 134.8
    if xH2O >= 0.012:
        T2 = fsolve(lambda x: compute(T1bar = x, xH2O = xH2O) - T1, 160, xtol = 1.E-5)[0]
    else:
        T2 = fsolve(lambda x: compute(T1bar = x, xH2O = xH2O) - T1, 138.5, xtol = 1.E-5)[0]
    #T2 = 144.8

    # connect to voyager data
    data = SaturnData()
    data.import_rss(pmin = 10.E3, pmax = 30.E3)
    # First diag, used to generate mean temperature profile connects to the Voyger observation
    diag = ThermoDiagram(pmin = 30.E3, pmax = 40.E5, np = 200, atmos = atmos)
    H2O = Water()
    array1 = hstack([log(diag.ptol), log(data.rss['lindal85'].pres)])[::-1]
    array2 = hstack([diag._moist_adiabats(atmos.T0, H2O), data.rss['lindal85'].temp])[::-1]

    # Second diag
    diag = ThermoDiagram(pmin = 10.E3, pmax = 40.E5, np = 200, atmos = atmos)
    temp_bar = amap(lambda p: spline(array1, array2, k = 2, s = 0)(p), log(diag.ptol))
    #temp_v1 = atmos._sat_temp_v(diag._moist_adiabats(T1, H2O), diag.ptol, H2O)
    temp_v1 = atmos._sat_temp_v(temp_bar, diag.ptol, H2O)
    temp_v2 = atmos._sat_temp_v(diag._moist_adiabats(T2, H2O), diag.ptol, H2O)
    icloud = find(temp_v2 - temp_v1 >= 0)[0]
    itop = find(temp_v2 - temp_v1 >= 0)[-1]
    pcloud = diag.ptol[icloud]
    ptop = diag.ptol[itop]
    buoy = atmos.grav * (temp_v2 - temp_v1) / temp_v1
    CAPE = trapz((buoy * temp_bar / atmos.T0)[icloud:itop], diag.zlev[icloud:itop])
    sigma = 5.67E-8
    coolt = atmos.cp * (pcloud - ptop) / (3 * sigma * atmos.grav) * (T1**-3 - T2**-3) \
            / (365 * 24 * 3600.)

    if (ax != 0): ax.plot(buoy, diag.ptol, 'b', linewidth = 2)
    return T2, pcloud, CAPE, coolt

def draw_thermo_diagram():
    ax      =   draw_prepare()
    T1bar   =   134.8
    xH2O    =   0.010
    compute(ax, T1bar, xH2O)
    show()

def draw_buoyancy():
    ax      =   draw_prepare()
    xH2O    =   [0.010, 0.012, 0.014, 0.016]
    for x in xH2O:
        print x
        cape_cin_rad(ax = ax, xH2O = x)
    pval = array(range(1, 10) + range(10, 20, 2) + range(20, 50, 5))
    pval = hstack([arange(0.1, 1, 0.2), pval])
    ax.set_yscale('log')
    ax.set_yticks(pval * 1E5)
    ax.set_yticklabels(pval)
    ax.set_ylabel('Pressure (bar)', fontsize = 20, weight = 'bold')
    ax.set_xlabel('Buoyancy ($\mathbf{m / s^2}$)', fontsize = 20, weight = 'bold')
    ax.set_xlim(0, 2)
    ax.set_ylim(30.E5, 0.1E5)
    ax2 = ax.twinx()
    ax2.set_ylim(-176.7, 119.6)
    ax2.set_ylabel('$\mathbf{Z^*}$ (km)', fontsize = 20, weight = 'bold')
    #ax.legend(['x = 0.010', 'x = 0.012', 'x = 0.014', 'x = 0.016'], loc = 0, frameon = False)
    show()

def draw_scaling():
    xH2O    =   linspace(0.010, 0.016, 7)
    #xH2O    =   [0.010, 0.012, 0.014, 0.016]
    data    =   []
    TEMP    =   10.
    CAPE    =   2.E5
    COOL    =   30.

    for x in xH2O:
        print x
        data.append(cape_cin_rad(xH2O = x))
    data = array(data)
    print data
    ax = draw_prepare()
    ax.plot(xH2O, data[:, 0] - 134.8, 'ro-', linewidth = 2)
    ax2 = ax.twinx()
    ax2.plot(xH2O, data[:, 3], 'bo-', linewidth = 2)
    ax.set_xlabel('Water Mixing ratio', fontsize = 20, weight = 'bold')
    ax.set_ylabel('Temperature increase at 1 bar (K)', fontsize = 20, weight = 'bold')
    ax2.set_ylabel('Cycling time (year)', fontsize = 20, weight = 'bold')
    ax.set_ylim([0, 25])
    show()

if __name__ == '__main__':
    #draw_thermo_diagram()
    #draw_buoyancy()
    #draw_scaling()
    atmos = Atmosphere(P0 = 1.E5, T0 = 134.8, grav = 10.44, xHe = 0.034)
    print atmos.R, atmos.H0, atmos.cp
