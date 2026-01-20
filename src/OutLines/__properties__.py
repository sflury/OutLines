from numpy import trapezoid,array,interp,log,log10,pi,logspace,argmax
from scipy.optimize import brentq,fminbound
from scipy.integrate import fixed_quad
from OutLines.__funcs__ import *
# hard-wired constants
global c,mH,mu,Msun,Lsun,yr,pi4,kpc
c    = 2.99792458e5 # km/s
mH   = 1.6733e-24 # g
mu   = 0.6 # mean molecular weight
Msun = 1.99e33 # g
Lsun = 3.9e33 # erg/s
yr   = 365.25*86400 # s/yr
kpc  = 3.086e21 # cm
pi4  = 4*pi

class Properties(object):
    '''
    Name:
        Properties

    Purpose:
        Inherit an OutLines Profile object and calculate outflow properties
        based on the characteristic parameters of the line profile. Results are
        stored as attributes and can be viewed using the \'print_props\' method.

    Arguments:
        :model (*OutLine object*): `OutLines` model profile for emission or
                        absorption

    Attributes:
        :props (*dict*): dictionary of the outflow properties
        :units (*dict*): dictionary of the outflow properties' units

    Methods:
        :update_params: updates parameters stored in the Properties object
        :calc_props: calculates outflow properties from parameters; called by
                        __init__ and update_params
        :print_props: prints outflow properties to the command line
    '''
    def __init__(self,Profile):
        # inherit the wind
        vars(self).update(vars(Profile))
        self.set_pars()
        self.calc_props()
    # print documentation
    def docs(self):
        print(self.__doc__)
    # sets parameters from the provided Profile model
    def set_pars(self):
        self.vinf = self.params['TerminalVelocity']*c
        self.beta = self.params['VelocityIndex']
        if self.settings['Geometry'] == 'Spherical':
            self.t0 = pi
            self.t1 = 0
        elif self.settings['Geometry'] == 'FilledCones':
            self.inc = self.params['Inclination']
            self.t0  = self.params['OpeningAngle']
            self.t1  = 0.
        elif self.settings['Geometry'] == 'HollowConesFixedCavity':
            self.inc = self.params['Inclination']
            self.t0  = self.params['OpeningAngle']
            self.t1  = self.params['OpeningAngle']-0.174533
        else:
            self.inc = self.params['Inclination']
            self.t0 = self.params['OpeningAngle']
            self.t1 = self.params['CavityAngle']
        if self.settings['FromRest']:
            self.vini = 0
        else:
            self.vini = self.params['LaunchVelocity']
        self.denpars = array([self.params[param] for param in self.npname])
    # options for updating the parameters
    def update_params(self,par_name,par_val):
        if hasattr(par_val,'__len__'):
            for name,val in zip(par_name,par_val):
                self.update_params(name,val)
        else:
            # if velocity > 1, convert to c units
            if 'Terminal' in par_name or 'Doppler' in par_name:
                if par_val > 1 :
                    par_val = par_val/2.99792458e5
            if 'Inclination' in par_name or 'Angle' in par_name:
                if par_val > 2*pi :
                    par_val = par_val*pi/180
            self.params[par_name] = par_val
            if self.settings['StaticComponent'] :
                if 'Static' in par_name or 'Doppler' in par_name:
                    self.statps[par_name] = par_val
                else:
                    self.outfps[par_name] = par_val
        self.set_pars()
        self.calc_props()
    # density profile
    def den(self,w1):
        return n[self.settings['DensityProfile']](w1,self.beta,self.vini,\
                        self.settings['VelocityField'],*self.denpars)
    # characteristic outflow velocity via Flury+ 2023
    def vel(self,x1):
        return v[self.settings['VelocityField']](x1,self.beta,self.vini)
    # outflow momentum density via Flury+ 2023
    def mom(self,x1):
        return self.den(self.vel(x1))*self.vel(x1)
    # calculate properties of the outflow
    def calc_props(self):
        # characteristic outflow radius
        xarr = logspace(0,2,2001)
        xout_guess = xarr[argmax(self.mom(xarr))]
        xout = fminbound(lambda x1: -self.mom(x1),1,2*xout_guess)
        # characteristic outflow velocity
        vout = self.vinf*self.vel(xout)
        # integrated density profile
        self.Rcal = lambda x1: self.den(self.vel(x1))
        Rcal = fixed_quad(self.Rcal,1,xout)[0]
        # mass outflow rate in km*Msun/yr (need to multiply by R0^2 n0)
        Mdot = pi4*(cos(self.t1)-cos(self.t0))*mu*mH*vout*Rcal*yr/Msun*kpc**2*1e5
        # relative momentum injection rate in dyne
        pdot = Mdot*vout*10**(30.7997-34)
        pdot_rel = ((cos(self.t1)-cos(self.t0))*Rcal*(vout/100)**2)
        #log10(pi4*(cos(self.t1)-cos(self.t0))*mu*mH*vout**2*Rcal * dyne)
        # relative energy injection rate in erg/s
        Edot =  (0.5*Mdot*vout**2) * 10**(35.7997-42)
        Edot_rel = ((cos(self.t1)-cos(self.t0))*Rcal*(vout/100)**3)
        self.props = {'x.out':xout,'v.out':vout,'Mdot':Mdot,'pdot':pdot,\
                      'Edot':Edot,'v.esc':vout/100,'pdot.esc':pdot_rel,\
                      'Edot.esc':Edot_rel}
        self.units = {'x.out':'','v.out':'km s^-1','Mdot':'Msun yr^-1 ',\
                    'pdot':'10^34 dyne','Edot':'10^42 erg s^-1',
                    'v.esc':'','pdot.esc':'','Edot.esc':''}
    # print outflow properties
    def print_props(self):
        print(3*' '+43*'-'+'\n  |'+13*' '+'MODEL PROPERTIES'+14*' '+'|\n'+3*' '+43*'-')
        for param,value in iter(self.props.items()):
            print(f'  | {param: >12s} : {value: >8.3f}  {self.units[param]: <16s} |')
        print(3*' '+43*'-')
        print('  | Mdot, pdot, Edot / R0^2 n0 [kpc^2 cm^-3]  |')
        print('  |  pdot.est, Edot.esc / v0 [100 km s^-1]    |\n'+3*' '+43*'-')
