from OutLines import __nebular__,__absorption__,__resonfluor__
from OutLines.__funcs__ import *
from OutLines.__params__ import *
# numpy functions for some basics
from numpy.random import randn
import sys
# wrapper to apply meta Profile object class to specific line profile cases
# with generalized methods and attributes
def profile_constructor(ProfileSubClass):
    '''
    Name:
        profile_constructor

    Purpose:
        Construct the Profile class given a particular line profile object

    Arguments:
        :ProfileSubClass (*class*): object class of `Nebular`, `Resonant`,
                `Fluorescent` or `Absorption`

    Returns:
        :Profile (*object*): An instance of the Profile object class inhereting
                the Nebular, Resonant, Fluorescent, or Absorption class with
                properties specific to the related spectral line

    '''
    class Profile(ProfileSubClass):
        '''
        Name:
            Profile

        Purpose:
            Generate a model line profile for the given outflow conditions
            at the given wavelength (and oscilator strength if applicable)

        Keyword Arguments:
            :VelocityField (*str*): radial velocity field indicating the assumed
                    treatment of the underlying physics and acceleration.
                    Options include
                         \'BetaCAK\' (1-1/x)^b from CAK theory - Castor et al. 1975
                         \'AccPlaw\' sqrt(1-x^(1-b)) from Steidel et al. 2010
                         \'VelPlaw\' 0.5(x-1)^b from various sources
            :Geometry (*str*): azimuthal geometry capturing the (an)isotropy
                    Options include
                        \'Sphere\' or \'Spherical\' isotropic case
                        \'FilledCones\' bidirectional cones
                        \'HollowCones\' or \'OpenCones\' cones with a cavity
                        \'HollowConesFixedCavity\' cones with fixed cavity size
            :DensityProfile (*str*): radial profile to use for the gas density.
                    Options include integration of many pulses or episodes
                         \'PowerLaw\'
                         \'Exponential\'
                         \'PowerLaw2\'
                     and individual pulse-like episodes or outbursts
                         \'LogNormal\'
                         \'Normal\'
                         \'Shell\'
                         \'FRED\'
                    or a combination of the two
                         \'Pulses\'
                         \'DampedPulses\'
            :Disk     (*bool*): boolean indicating whether a disk is present and
                    obstructs the cone posterior to the observer -- requires cone
                    geometry in order to take effect, default is False
            :FromRest (*bool*): boolean indicating whether the wind/shell/bubble
                    is launched from rest, default is True
            :Aperture (*bool*): boolean indicating whether the observed profile is
                    truncated due to the effects of a spherical aperture, default
                    is False

        Attributes:
            :Profile (*function*): function returning a line profile based on the chosen
                    velocity field, geometry, and density profile; takes wavelength and
                    paramaters (see the `params` attribute)
            :OutLinesModel (*function*): function returning the normalized
                    OutLines profile for a wind, bubble, or outflow without
                    fluxes, column densities, or other radiative transfer terms
            :npar (*int*): number of assigned parameters
            :params (*dict*): dictionary of parameters, initially assigned or later
                    updated by the user
            :bounds (*dict*): dictionary of parameter bounds, initially assigned or
                    later updated by the user

        Bound methods for viewing, updating, and computing the line profile:
            :get_profile: returns a model line profile given a wavelength array
            :get_params: returns the current parameters as a list
            :get_param_names: returns parameter names as a list
            :get_bounds: returns the current boundaries set for the parameters
            :update_params: updates parameter(s) given the name(s) and value(s)
            :update_bounds: updates bound(s) given the name(s) and value(s)
            :print_params: prints current parameter values to the command line
            :print_bounds: prints current parameter bounds to the command line
            :print_settings: prints profile model settings to the command line
            :set_params: resets parameter values to the OutLines defaults
            :set_bounds: resets paramater bounds to the OutLines defaults

        Bound methods for obtaining velocity quantiles of the line profile:
            :velocity_quantiles: returns velocity quantiles
            :W80: returns the W80 (0.9-0.1) velocity quantile
            :FWHM: returns the 0.75-0.25 velocity quantile

        Bound methods for modeling an observed line profile with OutLines:
            :init_params: samples 2*npar+1 values of the current parameters
                    assuming some intrinsic dispersion (default of 2^-6);
                    designed to generate valid initial walkers for MCMC fits
            :log_prior: returns a uniform bounded prior on the parameters with
                    parameters bounds defined in the Profile attribute
            :log_likelihood: returns the chi squared log likelihood of the line
                    profile model generated by OutLines given wavelength,
                    flux, and flux uncertainties
            :log_probability: returns the log probability of the parameters
                    given the data; sums the log prior and the log likelihood;
                    can be passed to MCMC codes like \`emcee\`
        '''
        def __init__(self,*args,VelocityField='BetaCAK',DensityProfile='PowerLaw',\
                    Geometry='Spherical',AddStatic=False,Disk=False,Aperture=False,FromRest=True):

            # keyword arguments specifying the model
            self.settings = {'Profile':ProfileSubClass.__name__,\
                             'VelocityField':VelocityField,\
                             'DensityProfile':DensityProfile,\
                             'Geometry':Geometry.replace('Open','Hollow'),\
                             'StaticComponent':AddStatic,\
                             'FromRest':FromRest,\
                             'Aperture':Aperture,\
                             'Disk':Disk
                             }

            # inheret the wind
            super(Profile,self).__init__(*args)

            # set central wavelength(s)
            if hasattr(args[0],'__len__'):
                self.w0 = array(args[0])
                # if absorption is present, set the oscilator strength if given
                if ProfileSubClass.__name__ != 'Nebular':
                    if len(args[0]) == len(args[1]):
                        self.fosc = array(args[1])
                    # otherwise, inform the user and call it quits
                    else:
                        print('OutLines requires the same number of oscilator'+\
                            '\nstrengths as central wavelengths for '+\
                            '\nAbsorption, Resonant, and Fluorescent line profiles')
                        sys.exit()
                    if ProfileSubClass.__name__ != 'Absorption':
                        if ProfileSubClass.__name__ != 'Fluorescent' :
                            if len(args[1]) == len(args[2]):
                                self.pline = array(args[2])
                            else:
                                print('OutLines requires the same number of oscilator'+\
                                    '\nstrengths and channel emission fractions for '+\
                                    '\nResonant and P Cygni line profiles')
                                sys.exit()
                        if ProfileSubClass.__name__ == 'Fluorescent' :
                            if len(args[1]) == len(args[2]) and len(args[1]) == len(args[3]):
                                self.fres  = array(args[2])
                                self.pline = array(args[3])
                            else:
                                print('OutLines requires the same number of oscilator'+\
                                    '\nstrengths and channel emission fractions for '+\
                                    '\nFluorescent line profiles')
                                sys.exit()
            else:
                self.w0 = array([args[0]])
                # if absorption is present, set the oscilator strength if given
                if ProfileSubClass.__name__ != 'Nebular':
                    self.fosc = array([args[1]])
                    if ProfileSubClass.__name__ != 'Absorption':
                        if ProfileSubClass.__name__ != 'Fluorescent' :
                            self.pline = array([args[2]])
                        else:
                            self.fres  = array([args[2]])
                            self.pline = array([args[3]])
            # set the number of lines
            self.nLines = len(self.w0)

            # set parameters, bounds, and profile model
            self.set_params()
            self.set_profile()
        # print documentation
        def docs(self):
            print(self.__doc__)
        # setters for the profile function, parameters, and boundaries
        def set_profile(self):
            if self.settings['StaticComponent']:
                self.Profile = self.__WithStatic__
                self.get_outflow    = self.__get__outflow_profile__
                self.get_static     = self.__get__static_profile__
            else:
                self.Profile = self.__WithoutStatic__
        def set_params(self):
            # set parameters, names, and labels
            # starting with velocities and fluxes
            if self.settings['StaticComponent']:
                par_init = [3e-5,1e-3,*BetaPars[self.settings['VelocityField']]]
                par_init += [ \
                    *ProfPars['Static'][self.settings['Profile']],\
                    *ProfPars['Outflow'][self.settings['Profile']] ]
                par_name  = ['DopplerWidth','TerminalVelocity','VelocityIndex']
                par_name += [ \
                    *ProfName['Static'][self.settings['Profile']],\
                    *ProfName['Outflow'][self.settings['Profile']]]
                par_labs  = [r'$\sigma_v$ [km s$^{-1}$]',r'$v_\infty$ [km s$^{-1}$]',r'$\beta$']
                par_labs += [ \
                    *ProfLabs['Static'][self.settings['Profile']],\
                    *ProfLabs['Outflow'][self.settings['Profile']] ]
            else:
                par_init = [1e-3,*BetaPars[self.settings['VelocityField']]]
                par_init += [ \
                    *ProfPars['Outflow'][self.settings['Profile']] ]
                par_name  = ['TerminalVelocity','VelocityIndex']
                par_name += [ \
                    *ProfName['Outflow'][self.settings['Profile']] ]
                par_labs  = [r'$v_\infty$ [km s$^{-1}$]',r'$\beta$']
                par_labs += [ \
                    *ProfLabs['Outflow'][self.settings['Profile']] ]
            # add in geometry
            par_init += [*GeomPars[self.settings['Geometry']]]
            par_name += [*GeomName[self.settings['Geometry']]]
            par_labs += [*GeomLabs[self.settings['Geometry']]]
            # optional disk
            if self.settings['Disk'] :
                par_init += [2]
                par_name += ['DiskRadius']
                par_labs += [r'$x_{disk}$']
            # optional launch velocity
            if not self.settings['FromRest'] :
                par_init += [1e-6]
                par_name += ['LaunchVelocity']
                par_labs += [r'$v_{0}$']
            # optional aperture limit
            if self.settings['Aperture'] :
                par_init += [1e-3]
                par_name += ['ApertureVelocity']
                par_labs += [r'$v_{apert}$']
            par_init += [*DensPars[self.settings['DensityProfile']]]
            par_name += [*DensName[self.settings['DensityProfile']]]
            par_labs += [*DensLabs[self.settings['DensityProfile']]]
            # sort and store parameters in dictionaries
            self.params = {}
            self.labels = {}
            if self.settings['StaticComponent'] :
                self.statps = {}
                self.outfps = {}
            for name,val,lab in zip(par_name,par_init,par_labs):
                if 'Flux' in name or 'Column' in name :
                    for i in range(self.nLines):
                        self.params[f'{name}{i+1}'] = val
                        self.labels[f'{name}{i+1}'] = lab+f'$_{i+1}$'
                else:
                    self.params[name] = val
                    self.labels[name] = lab
                if self.settings['StaticComponent'] :
                    if 'Static' in name or 'Doppler' in name:
                        if 'Flux' in name or 'Column' in name:
                            for i in range(self.nLines):
                                self.statps[f'{name}{i+1}'] = val
                        else:
                            self.statps[name] = val
                    else:
                        if 'Flux' in name or 'Column' in name :
                            for i in range(self.nLines):
                                self.outfps[f'{name}{i+1}'] = val
                        else:
                            self.outfps[name] = val
            self.npar = len(self.params)
            self.gpname = GeomName[self.settings['Geometry']]
            self.npname = DensName[self.settings['DensityProfile']]
            # set parameter bounds
            bounds = [[],[]]
            for l in [0,1]:
                if self.settings['StaticComponent']:
                    bounds[l] += [float(l),float(l),\
                        *BetaBounds[self.settings['VelocityField']][l],\
                        *ProfBounds['Static'][self.settings['Profile']][l],\
                        *ProfBounds['Outflow'][self.settings['Profile']][l] ]
                else:
                    bounds[l] += [float(l),\
                        *BetaBounds[self.settings['VelocityField']][l],\
                        *ProfBounds['Outflow'][self.settings['Profile']][l] ]
                bounds[l] += [*GeomBounds[self.settings['Geometry']][l]]
                if self.settings['Disk'] :
                    bounds[l] += [*DiskBounds[l]]
                if not self.settings['FromRest'] :
                    bounds[l] += [float(l)]
                if self.settings['Aperture'] :
                    bounds[l] += [float(l)]
                bounds[l] += [*DensBounds[self.settings['DensityProfile']][l]]
            # sort and store parameter bounds in dictionary
            self.bounds = {}
            for name,lower,upper in zip(par_name,bounds[0],bounds[1]):
                if 'Flux' in name or 'Column' in name :
                    for i in range(self.nLines):
                        name.strip(f'{i}')
                        self.bounds[f'{name}{i+1}'] = [lower,upper]
                else:
                    self.bounds[name] = [lower,upper]
        # getters for the names, labels, paramaters, bounds, and profile
        def get_params(self,comp='all'):
            if comp.lower() == 'outflow' and self.settings['StaticComponent']:
                return array(list(self.outfps.values()))
            elif comp.lower() == 'static' and self.settings['StaticComponent']:
                return array(list(self.statps.values()))
            else:
                return array(list(self.params.values()))
        def get_param_names(self):
            return array(list(self.params.keys()))
        def get_param_labels(self):
            return array(list(self.labels.values()))
        def get_bounds(self,astype=tuple):
            return astype(array(list(self.bounds.values())).T)
        def get_profile(self,wavelengths):
            return self.Profile(wavelengths,*self.get_params())
        def __get__static_profile__(self,wavelengths):
            return self.__OnlyStatic__(wavelengths,*self.get_params(comp='static'))
        def __get__outflow_profile__(self,wavelengths):
            return self.__WithoutStatic__(wavelengths,*self.get_params(comp='outflow'))
        # options for updating the parameters and bounds
        def update_params(self,par_name,par_val):
            if hasattr(par_val,'__len__'):
                for name,val in zip(par_name,par_val):
                    self.update_params(name,val)
            else:
                # if velocity > 1, convert to c units
                if 'Terminal' in par_name or 'Launch' in par_name or 'Doppler' in par_name or 'Aperture' in par_name :
                    if par_val > 1 :
                        par_val = par_val/2.99792458e5
                if 'Inclination' in par_name or 'Angle' in par_name:
                    if par_val > pi/2 :
                        if par_val > 90 :
                            print('Angles must be between 0 and 90 degrees')
                        par_val = par_val*pi/180
                self.params[par_name] = par_val
                if self.settings['StaticComponent'] :
                    if 'Static' in par_name or 'Doppler' in par_name:
                        self.statps[par_name] = par_val
                    else:
                        self.outfps[par_name] = par_val
        # options for updating the parameters and bounds
        def update_bounds(self,par_name,bound_val):
            if hasattr(bound_val[0],'__len__'):
                for name,val in zip(par_name,bound_val):
                    self.update_bounds(name,val)
            else:
                # if velocity > 1, convert to c units
                if 'Terminal' in par_name or 'Launch' in par_name or 'Doppler' in par_name or 'Aperture' in par_name :
                    if bound_val[1] > 1 :
                        bound_val[0] = bound_val[0]/2.99792458e5
                        bound_val[1] = bound_val[1]/2.99792458e5
                if 'Inclination' in par_name or 'Angle' in par_name:
                    if bound_val[1] > pi/2 :
                        if bound_val[1] > 90 :
                            print('Angles must be between 0 and 90 degrees')
                            bound_val[0] = bound_val[0]%90
                            bound_val[1] = bound_val[1]%90
                        bound_val[0] = bound_val[0]*pi/180
                        bound_val[1] = bound_val[1]*pi/180

                self.bounds[par_name] = bound_val
        # convenience function for printing settings to terminal
        def print_settings(self):
            print(3*' '+36*'-'+'\n  |'+11*' '+'MODEL SETTINGS'+11*' '+'|\n'+3*' '+36*'-')
            for i,w0 in enumerate(self.w0):
                print(f'  |     Line  {i+1: >2d}     : {w0: <15.3f} |')
            for option,setting in iter(self.settings.items()):
                if type(setting) == bool :
                    if setting:
                        print(f'  | {option: >16s} : {"Yes": <15s} |')
                    else:
                        print(f'  | {option: >16s} : {"No": <15s} |')
                elif 'Fixed' in setting :
                    s1 = setting.split('Fixed')[0]
                    s2 = 'Fixed'+setting.split('Fixed')[1]
                    print(f'  | {option: >16s} : {s1: <15s} |')
                    print(f'  | {"": >16s}   {s2: <15s} |')
                else:
                    print(f'  | {option: >16s} : {setting: <15s} |')
            print(3*' '+36*'-')
        # convenience function for printing parameters to terminal
        def print_params(self):
            print(3*' '+36*'-'+'\n  |'+10*' '+'MODEL PARAMETERS'+10*' '+'|\n'+3*' '+36*'-')
            for param,value in iter(self.params.items()):
                if 'Terminal' in param or 'Doppler' in param or 'Aperture' in param:
                    print(f'  | {param: >18s} : {value*2.99792458e5: >8.3f} km/s |')
                elif 'Inclination' in param or 'Angle' in param:
                    print(f'  | {param: >18s} : {value*180/pi: >8.3f}°     |')
                else:
                    print(f'  | {param: >18s} : {value: >8.3f}      |')
            print(3*' '+36*'-')
        # convenience function for printing parameter bounds to terminal
        def print_bounds(self):
            delim = [',','|']
            print(3*' '+61*'-'+'\n  |'+20*' '+'PARAMETER  BOUNDARIES'+20*' '+'|\n'+3*' '+61*'-')
            for param,bound in iter(self.bounds.items()):
                line = f'  | {param: >18s} : '
                for l in [0,1]:
                    if 'Terminal' in param or 'Doppler' in param or 'Aperture' in param:
                        line += f'{bound[l]*2.99792458e5: >13.3f} km/s {delim[l]}'
                    elif 'Inclination' in param or 'Angle' in param:
                        line += f'{bound[l]*180/pi: >13.3f}°     {delim[l]}'
                    else:
                        line += f'{bound[l]: >13.3f}      {delim[l]}'
                print(line)
            print(3*' '+61*'-')
        # compute velocity quantiles
        def velocity_quantiles(self,quantiles=array([0.1,0.9]),verbose=True):
            uscl = 101
            uarr = array(list(map(lambda i: float(i)/uscl,range(-uscl,uscl+1,1))))
            wave = self.w0[0] * (1+uarr*self.params['TerminalVelocity'])
            vels = uarr*2.99792458e5*self.params['TerminalVelocity']
            if self.settings['Profile'] != 'Absorption' :
                cdf  = cumulative_trapezoid(self.get_profile(wave),x=wave)
            elif self.settings['Profile'] == 'Absorption' :
                cdf  = cumulative_trapezoid(1-self.get_profile(wave),x=wave)
            else:
                print('Velocity Quantiles Not Supported for P Cygni.')
                print('Try again using Nebular or Absorption profiles')
                print('with the same parameters.')
                return nan
            cdf /= cdf[-1]
            vel_quant = interp(quantiles,cdf,vels[1:])
            if verbose:
                print(3*' '+34*'-'+'\n  |'+8*' '+'VELOCITY QUANTILES'+8*' '+'|\n'+3*' '+34*'-')
                for qi,vq in zip(quantiles,vel_quant):
                    print('  |   '+f'{qi: >8.6f}  :  {vq: >8.3f} km s^-1  |')
                print(3*' '+34*'-')
            return vel_quant
        # convenience functions for different common velocity quantiles
        def get_W80(self):
            vq = self.velocity_quantiles(quantiles=[0.1,0.9],verbose=False)
            return vq[1]-vq[0]
        def get_W50(self):
            vq = self.velocity_quantiles(quantiles=[0.25,0.75],verbose=False)
            return vq[1]-vq[0]
        def get_FWHF(self):
            return self.get_W50()
        def get_1sigma(self):
            vq = self.velocity_quantiles(quantiles=[0.1587,0.8413],verbose=False)
            return vq[1]-vq[0]
        # initiate parameters
        def init_params(self,disp=2**-6):
            nwalk = 2*self.npar+1
            init = self.get_params()*(1 + disp*randn(nwalk,self.npar))
            for i in range(nwalk):
                for j in range(self.npar):
                    while init[i,j] < self.get_bounds()[0][j] or \
                                            init[i,j] > self.get_bounds()[1][j]:
                        init[i,j] = self.get_params()[j]*(1 + disp*randn())
                if self.settings['Geometry'] == 'HollowCones' :
                    j = where(self.get_param_names()=='OpeningAngle')[0]
                    if init[i,j+1] >= init[i,j] - 0.0872664626 :
                        init[i,j+1] = init[i,j] - 0.0872664626+disp*randn()
            return init,nwalk
        # log likelihood for chi squared
        def log_likelihood(self,theta,x,y,yerr):
            ybar = self.Profile(x,*theta)
            ln_prb = -0.5*((y-ybar)/(yerr))**2-2*log(yerr)
            return nansum(ln_prb)
        # boundary condition check
        @staticmethod
        def __fun_bound__(t,tl,tu):
            return (tl<t)&(t<tu)&(isfinite(t))
        # log prior -- uniform bounded
        def log_prior(self,theta):
            tlower,tupper = self.get_bounds()
            tbound = list(map(self.__fun_bound__,theta,tlower,tupper))
            if all(tbound):
                return 0
            else:
                return -inf
        # log probability -- prior plus likelihood
        def log_probability(self,theta,x,y,var):
            lp = self.log_prior(theta)
            if not isfinite(lp):
                return -inf
            else:
                ll = self.log_likelihood(theta,x,y,var)
                if not isfinite(ll):
                    return -inf
                else:
                    return lp + ll

    return Profile
# various "parent" line profile classes
@profile_constructor
class Nebular():
    '''
    Name:
        Nebular

    Purpose:
        Compute a model nebular emission line profile for outflowing gas for the
        given atomic data and assumed velocity field, geometry, and density profile.

    Arguments:
        :w0 (*float*): rest cental wavelength of the line

    Keyword arguments:
            :VelocityField (*str*): radial velocity field indicating the assumed
                    treatment of the underlying physics and acceleration.
                    Options include
                         \'BetaCAK\' (1-1/x)^b from CAK theory
                         \'AccPlaw\' sqrt(1-x^(1-b)) from Steidel et al. 2010
                         \'VelPlaw\' 0.5(x-1)^b from various sources
            :Geometry (*str*): azimuthal geometry capturing the (an)isotropy
                    Options include
                        \'Sphere\' or \'Spherical\' isotropic case
                        \'FilledCones\' bidirectional cones
                        \'HollowCones\' or \'OpenCones\' cones with a cavity
                        \'HollowConesFixedCavity\' cones with fixed cavity size
            :DensityProfile (*str*): radial profile to use for the gas density.
                    Options include integration of many pulses or episodes
                         \'PowerLaw\'
                         \'Exponential\'
                         \'PowerLaw2\'
                     and individual pulse-like episodes or outbursts
                         \'LogNormal\'
                         \'Normal\'
                         \'Shell\'
                         \'FRED\'
                    or a combination of the two
                         \'Pulses\'
                         \'DampedPulses\'
            :Disk     (*bool*): boolean indicating whether a disk is present and
                    obstructs the cone posterior to the observer -- requires cone
                    geometry in order to take effect
            :Aperture (*bool*): boolean indicating whether the observed profile is
                    truncated due to the effects of a spherical aperture

    Returns:
        :Profile (*object*): An instance of the Profile object class inhereting
                the Nebular class with properties specific to an emission line

    Attributes:
        :Profile (*function*): function returning a line profile based on the chosen
                velocity field, geometry, and density profile; takes wavelength and
                paramaters (see the `params` attribute)
        :params (*dict*): dictionary of parameters, initially assigned or later
                updated by the user
        :npar (*int*): number of assigned parameters
        :bounds (*dict*): dictionary of parameter bounds, initially assigned or
                later updated by the user

    Bound methods:
        :get_params: returns the current parameters
        :update_params: updates parameter(s) given the name(s) and value(s)
        :get_bounds: returns the current boundaries set for the parameters
        :update_bounds: updates bound(s) given the name(s) and value(s)
        :print_params: prints current parameter values to the command line

    '''
    def __init__(self,w0):
        kwargs = {k:self.settings[k] for k in \
            ['VelocityField','DensityProfile','Geometry','Disk','FromRest','Aperture']}
        self.OutLinesModel = __nebular__.build_profile_model(**kwargs)
        pass
    # define profile with and without a static ISM component
    # order of parameters:
    # static flux(es), Doppler width, outflow flux(es), terminal velocity
    # velocity field index, geometry terms, density profile terms
    def __WithStatic__(self,w,sigv,vinf,beta,*pars):
        Fs = pars[:self.nLines]
        Fo = pars[self.nLines:2*self.nLines]
        po = pars[2*self.nLines:]
        OutLine = array(list(map(lambda w0: \
          self.OutLinesModel(w,w0,vinf,beta,*po), self.w0 )))
        static  = array(list(map(lambda w0: gauss(w,w0,sigv), self.w0 )))
        return Fs@static + Fo@OutLine
    def __WithoutStatic__(self,w,vinf,beta,*pars):
        Fo = pars[:self.nLines]
        po = pars[self.nLines:]
        OutLine = array(list(map(lambda w0: \
          self.OutLinesModel(w,w0,vinf,beta,*po), self.w0 )))
        return Fo@OutLine
    def __OnlyStatic__(self,w,sigv,*Fs):
        static  = array(list(map(lambda w0: gauss(w,w0,sigv), self.w0 )))
        return Fs@static
@profile_constructor
class Absorption():
    '''
    Name:
        Absorption

    Purpose:
        Compute a model spectral absorption line profile for outflowing gas for the
        given atomic data and assumed velocity field, geometry, and density profile.

    Arguments:
        :w0 (*float*): rest cental wavelength of the line
        :fosc (*float*): oscilator strength of the line

    Keyword arguments:
            :VelocityField (*str*): radial velocity field indicating the assumed
                    treatment of the underlying physics and acceleration.
                    Options include
                         \'BetaCAK\' (1-1/x)^b from CAK theory
                         \'AccPlaw\' sqrt(1-x^(1-b)) from Steidel et al. 2010
                         \'VelPlaw\' 0.5(x-1)^b from various sources
            :Geometry (*str*): azimuthal geometry capturing the (an)isotropy
                    Options include
                        \'Sphere\' or \'Spherical\' isotropic case
                        \'FilledCones\' bidirectional cones
                        \'HollowCones\' or \'OpenCones\' cones with a cavity
                        \'HollowConesFixedCavity\' cones with fixed cavity size
            :DensityProfile (*str*): radial profile to use for the gas density.
                    Options include integration of many pulses or episodes
                         \'PowerLaw\'
                         \'Exponential\'
                         \'PowerLaw2\'
                     and individual pulse-like episodes or outbursts
                         \'LogNormal\'
                         \'Normal\'
                         \'Shell\'
                         \'FRED\'
                    or a combination of the two
                         \'Pulses\'
                         \'DampedPulses\'
            :Disk     (*bool*): boolean indicating whether a disk is present and
                    obstructs the cone posterior to the observer -- requires cone
                    geometry in order to take effect
            :Aperture (*bool*): boolean indicating whether the observed profile is
                    truncated due to the effects of a spherical aperture

    Returns:
        :Profile (*object*): An instance of the Profile object class inhereting
                the Absorption class with properties specific to an absorption line

    Attributes:
        :Profile (*function*): function returning a line profile based on the chosen
                velocity field, geometry, and density profile; takes wavelength and
                paramaters (see the `params` attribute)
        :params (*dict*): dictionary of parameters, initially assigned or later
                updated by the user
        :npar (*int*): number of assigned parameters
        :bounds (*dict*): dictionary of parameter bounds, initially assigned or
                later updated by the user

    Bound methods:
        :get_params: returns the current parameters
        :update_params: updates parameter(s) given the name(s) and value(s)
        :get_bounds: returns the current boundaries set for the parameters
        :update_bounds: updates bound(s) given the name(s) and value(s)
        :print_params: prints current parameter values to the command line

    '''
    def __init__(self,w0,fosc):
        kwargs = {k:self.settings[k] for k in \
            ['VelocityField','DensityProfile','Geometry','Disk','FromRest','Aperture']}
        self.OutLinesModel = __absorption__.build_profile_model(**kwargs)
        pass
    # define profile with and without a static ISM component
    def __WithStatic__(self,w,sigv,vinf,beta,*pars):
        OutLine = array(list(map(lambda w0: \
          self.OutLinesModel(w,w0,vinf,beta,*pars[2*self.nLines:]),self.w0)))
        static  = array(list(map(lambda w0: voigt(w,w0,sigv),self.w0)))
        # classical absorption cross-section
        # log(sigma) = -14.8247 in cm^2, b in km/s, w in Ang
        # speed of light
        # c = 2.99792458e5 km/s
        logNs = array(pars[:self.nLines])
        logNo = array(pars[self.nLines:2*self.nLines])
        taus = -self.w0*self.fosc*10**(logNs-14.8247)/(sigv*2.99792458e5)
        tauo = -self.w0*self.fosc*10**(logNo-14.8247)/(vinf*2.99792458e5)
        return exp( taus@static + tauo@OutLine )
    def __WithoutStatic__(self,w,vinf,beta,*pars):
        OutLine = array(list(map(lambda w0: \
          self.OutLinesModel(w,w0,vinf,beta,*pars[self.nLines:]),self.w0)))
        # classical absorption cross-section
        # log(sigma) = -14.8247 in cm^2, b in km/s, w in Ang
        # speed of light
        # c = 2.99792458e5 km/s
        logNo = array(pars[:self.nLines])
        tauo = -self.w0*self.fosc*10**(logNo-14.8247)/(vinf*2.99792458e5)
        return exp( tauo@OutLine )
    def __OnlyStatic__(self,w,sigv,*pars):
        static  = array(list(map(lambda w0: voigt(w,w0,sigv),self.w0)))
        # classical absorption cross-section
        # log(sigma) = -14.8247 in cm^2, b in km/s, w in Ang
        # speed of light
        # c = 2.99792458e5 km/s
        logNs = array(pars[:self.nLines])
        taus = -self.w0*self.fosc*10**(logNs-14.8247)/(sigv*2.99792458e5)
        return exp( taus@static )
@profile_constructor
class Resonant():
    '''
    Name:
        Resonant

    Purpose:
        Compute a model spectral resonant line profile for
        outflowing gas for the given atomic data and assumed velocity field,
        geometry, and density profile.

    Arguments:
        :w0 (*float*): rest cental wavelength of the line
        :fosc (*float*): oscilator strength of the line
        :pline (*float*): relative spontaneous transition probability of the
                    line through the resonant channel:
                        p_r = A_r / (A_f + A_r)

    Keyword arguments:
            :VelocityField (*str*): radial velocity field indicating the assumed
                    treatment of the underlying physics and acceleration.
                    Options include
                         \'BetaCAK\' (1-1/x)^b from CAK theory
                         \'AccPlaw\' sqrt(1-x^(1-b)) from Steidel et al. 2010
                         \'VelPlaw\' 0.5(x-1)^b from various sources
            :Geometry (*str*): azimuthal geometry capturing the (an)isotropy
                    Options include
                        \'Sphere\' or \'Spherical\' isotropic case
                        \'FilledCones\' bidirectional cones
                        \'HollowCones\' or \'OpenCones\' cones with a cavity
                        \'HollowConesFixedCavity\' cones with fixed cavity size
            :DensityProfile (*str*): radial profile to use for the gas density.
                    Options include integration of many pulses or episodes
                         \'PowerLaw\'
                         \'Exponential\'
                         \'PowerLaw2\'
                     and individual pulse-like episodes or outbursts
                         \'LogNormal\'
                         \'Normal\'
                         \'Shell\'
                         \'FRED\'
                    or a combination of the two
                         \'Pulses\'
                         \'DampedPulses\'
            :Disk     (*bool*): boolean indicating whether a disk is present and
                    obstructs the cone posterior to the observer -- requires cone
                    geometry in order to take effect
            :Aperture (*bool*): boolean indicating whether the observed profile is
                    truncated due to the effects of a spherical aperture

    Returns:
        :Profile (*object*): An instance of the Profile object class inhereting
                the Absorption class with properties specific to an absorption line

    Attributes:
        :Profile (*function*): function returning a line profile based on the chosen
                velocity field, geometry, and density profile; takes wavelength and
                paramaters (see the `params` attribute)
        :params (*dict*): dictionary of parameters, initially assigned or later
                updated by the user
        :npar (*int*): number of assigned parameters
        :bounds (*dict*): dictionary of parameter bounds, initially assigned or
                later updated by the user

    Bound methods:
        :get_params: returns the current parameters
        :update_params: updates parameter(s) given the name(s) and value(s)
        :get_bounds: returns the current boundaries set for the parameters
        :update_bounds: updates bound(s) given the name(s) and value(s)
        :print_params: prints current parameter values to the command line

    '''
    def __init__(self,w0,fosc,pline):
        kwargs = {k:self.settings[k] for k in \
            ['VelocityField','DensityProfile','Geometry','Disk','FromRest','Aperture']}
        self.OutLinesModel = __resonfluor__.build_profile_model(**kwargs)
        pass
    # define profile with and without a static ISM component
    def __WithStatic__(self,w,sigv,vinf,beta,*pars):
        OutLine = array(list(map(lambda w0: \
          self.OutLinesModel(w,w0,vinf,beta,*pars[2*self.nLines:]),self.w0)))
        staticv  = array(list(map(lambda w0: voigt(w,w0,sigv),self.w0)))
        #staticg  = array(list(map(lambda w0: gauss(w,w0,sigv),self.w0)))
        # classical absorption cross-section
        # log(sigma) = -14.8247 in cm^2, b in km/s, w in Ang
        # speed of light
        # c = 2.99792458e5 km/s
        logNs = array(pars[:self.nLines])
        logNo = array(pars[self.nLines:2*self.nLines])
        # static
        tausr = self.w0*self.fosc*10**(logNs-14.8247)/(vinf*2.99792458e5)
        fescs = array(list(map(self.__EscapeFraction__,tausr,self.pline,staticv)))
        srcfs = array(list(map(self.__SourceFunction__,tausr,staticv)))
        # OutLine
        tauor = self.w0*self.fosc*10**(logNo-14.8247)/(vinf*2.99792458e5)
        fesco = array(list(map(self.__EscapeFraction__,tauor,self.pline,OutLine)))
        srcfo = array(list(map(self.__SourceFunction__,tauor,OutLine)))
        return (fesco*srcfo).sum(axis=0) + (fescs*srcfs).sum(axis=0)
    def __WithoutStatic__(self,w,vinf,beta,*pars):
        OutLine = array(list(map(lambda w0: \
          self.OutLinesModel(w,w0,vinf,beta,*pars[self.nLines:]),self.w0)))
        # classical absorption cross-section
        # log(sigma) = -14.8247 in cm^2, b in km/s, w in Ang
        # speed of light c = 2.99792458e5 km/s
        logNo = array(pars[:self.nLines])
        tauor = self.w0*self.fosc*10**(logNo-14.8247)/(vinf*2.99792458e5)
        fesco = array(list(map(self.__EscapeFraction__,tauor,self.pline,OutLine)))
        srcfo = array(list(map(self.__SourceFunction__,tauor,OutLine)))
        return (fesco*srcfo).sum(axis=0)
    def __OnlyStatic__(self,w,sigv,*pars):
        staticv  = array(list(map(lambda w0: voigt(w,w0,sigv),self.w0)))
        #staticg  = array(list(map(lambda w0: gauss(w,w0,sigv),self.w0)))
        # classical absorption cross-section
        # log(sigma) = -14.8247 in cm^2, b in km/s, w in Ang
        # speed of light
        # c = 2.99792458e5 km/s
        logNs = array(pars[:self.nLines])
        tausr = self.w0*self.fosc*10**(logNs-14.8247)/(vinf*2.99792458e5)
        fescs = array(list(map(self.__EscapeFraction__,tausr,self.pline,staticv)))
        srcfs = array(list(map(self.__SourceFunction__,tausr,staticv)))
        return (fescs*srcfs).sum(axis=0)
    @staticmethod
    def __EscapeFraction__(tau,pline,profile):
        pesc = (1-exp(-tau*profile))/(tau*profile)
        #pesc[isnan(pesc)] = 1
        fesc = pline*pesc/(1-pline*(1-pesc))
        fesc[profile<2**-24] = pline
        return fesc
    @staticmethod
    def __SourceFunction__(tau,profile):
        return 1-exp(-tau*profile)
@profile_constructor
class Fluorescent():
    '''
    Name:
        Fluorescent

    Purpose:
        Compute a model spectral fluorescent line profile for
        outflowing gas for the given atomic data and assumed velocity field,
        geometry, and density profile.

    Arguments:
        :w0 (*float*): rest cental wavelength of the line
        :fosc (*float*): oscilator strength of the line
        :pline (*float*): relative spontaneous transition probability of the
                    line through the fluorescent channel:
                        p_f = A_f / (A_f + A_r)

    Keyword arguments:
            :VelocityField (*str*): radial velocity field indicating the assumed
                    treatment of the underlying physics and acceleration.
                    Options include
                         \'BetaCAK\' (1-1/x)^b from CAK theory
                         \'AccPlaw\' sqrt(1-x^(1-b)) from Steidel et al. 2010
                         \'VelPlaw\' 0.5(x-1)^b from various sources
            :Geometry (*str*): azimuthal geometry capturing the (an)isotropy
                    Options include
                        \'Sphere\' or \'Spherical\' isotropic case
                        \'FilledCones\' bidirectional cones
                        \'HollowCones\' or \'OpenCones\' cones with a cavity
                        \'HollowConesFixedCavity\' cones with fixed cavity size
            :DensityProfile (*str*): radial profile to use for the gas density.
                    Options include integration of many pulses or episodes
                         \'PowerLaw\'
                         \'Exponential\'
                         \'PowerLaw2\'
                     and individual pulse-like episodes or outbursts
                         \'LogNormal\'
                         \'Normal\'
                         \'Shell\'
                         \'FRED\'
                    or a combination of the two
                         \'Pulses\'
                         \'DampedPulses\'
            :Disk     (*bool*): boolean indicating whether a disk is present and
                    obstructs the cone posterior to the observer -- requires cone
                    geometry in order to take effect
            :Aperture (*bool*): boolean indicating whether the observed profile is
                    truncated due to the effects of a spherical aperture

    Returns:
        :Profile (*object*): An instance of the Profile object class inhereting
                the Absorption class with properties specific to an absorption line

    Attributes:
        :Profile (*function*): function returning a line profile based on the chosen
                velocity field, geometry, and density profile; takes wavelength and
                paramaters (see the `params` attribute)
        :params (*dict*): dictionary of parameters, initially assigned or later
                updated by the user
        :npar (*int*): number of assigned parameters
        :bounds (*dict*): dictionary of parameter bounds, initially assigned or
                later updated by the user

    Bound methods:
        :get_params: returns the current parameters
        :update_params: updates parameter(s) given the name(s) and value(s)
        :get_bounds: returns the current boundaries set for the parameters
        :update_bounds: updates bound(s) given the name(s) and value(s)
        :print_params: prints current parameter values to the command line

    '''
    def __init__(self,w0,fosc_flu,fosc_res,pline):
        kwargs = {k:self.settings[k] for k in \
            ['VelocityField','DensityProfile','Geometry','Disk','FromRest','Aperture']}
        self.OutLinesModel = __resonfluor__.build_profile_model(**kwargs)
        pass
    # define profile with and without a static ISM component
    def __WithStatic__(self,w,sigv,vinf,beta,*pars):
        OutLine = array(list(map(lambda w0: \
          self.OutLinesModel(w,w0,vinf,beta,*pars[2*self.nLines:]),self.w0)))
        staticv  = array(list(map(lambda w0: voigt(w,w0,sigv),self.w0)))
        #staticg  = array(list(map(lambda w0: gauss(w,w0,sigv),self.w0)))
        # classical absorption cross-section
        # log(sigma) = -14.8247 in cm^2, b in km/s, w in Ang
        # speed of light
        # c = 2.99792458e5 km/s
        logNs = array(pars[:self.nLines])
        logNo = array(pars[self.nLines:2*self.nLines])
        # static
        tausr = self.w0*self.fres*10**(logNs-14.8247)/(vinf*2.99792458e5)
        tausf = self.w0*self.fosc*10**(logNs-14.8247)/(vinf*2.99792458e5)
        fescs = array(list(map(self.__EscapeFraction__,tausr,self.pline,staticv)))
        srcfs = array(list(map(self.__SourceFunction__,tausf,staticv)))
        # OutLine
        tauor = self.w0*self.fres*10**(logNo-14.8247)/(vinf*2.99792458e5)
        tauof = self.w0*self.fosc*10**(logNo-14.8247)/(vinf*2.99792458e5)
        fesco = array(list(map(self.__EscapeFraction__,tauor,self.pline,OutLine)))
        srcfo = array(list(map(self.__SourceFunction__,tauof,OutLine)))
        return (fesco*srcfo).sum(axis=0) + (fescs*srcfs).sum(axis=0)
    def __WithoutStatic__(self,w,vinf,beta,*pars):
        OutLine = array(list(map(lambda w0: \
          self.OutLinesModel(w,w0,vinf,beta,*pars[self.nLines:]),self.w0)))
        # classical absorption cross-section
        # log(sigma) = -14.8247 in cm^2, b in km/s, w in Ang
        # speed of light c = 2.99792458e5 km/s
        logNo = array(pars[:self.nLines])
        tauor = self.w0*self.fres*10**(logNo-14.8247)/(vinf*2.99792458e5)
        tauof = self.w0*self.fosc*10**(logNo-14.8247)/(vinf*2.99792458e5)
        fesco = array(list(map(self.__EscapeFraction__,tauor,self.pline,OutLine)))
        srcfo = array(list(map(self.__SourceFunction__,tauof,OutLine)))
        return (fesco*srcfo).sum(axis=0)
    def __OnlyStatic__(self,w,sigv,*pars):
        staticv  = array(list(map(lambda w0: voigt(w,w0,sigv),self.w0)))
        #staticg  = array(list(map(lambda w0: gauss(w,w0,sigv),self.w0)))
        # classical absorption cross-section
        # log(sigma) = -14.8247 in cm^2, b in km/s, w in Ang
        # speed of light
        # c = 2.99792458e5 km/s
        logNs = array(pars[:self.nLines])
        tausr = self.w0*self.fres*10**(logNs-14.8247)/(vinf*2.99792458e5)
        tausf = self.w0*self.fosc*10**(logNs-14.8247)/(vinf*2.99792458e5)
        fescs = array(list(map(self.__EscapeFraction__,tausr,self.pline,staticv)))
        srcfs = array(list(map(self.__SourceFunction__,tausf,staticv)))
        return (fescs*srcfs).sum(axis=0)
    @staticmethod
    def __EscapeFraction__(tau,pline,profile):
        pesc = (1-exp(-tau*profile))/(tau*profile)
        fesc = pline/(1-(1-pline)*(1-pesc))
        fesc[profile<2**-24] = pline
        return fesc
    @staticmethod
    def __SourceFunction__(tau,profile):
        return 1-exp(-tau*profile)
@profile_constructor
class PCygni():
    '''
    Name:
        PCygni

    Purpose:
        Compute a model spectral resonant P Cygni line profile for
        outflowing gas for the given atomic data and assumed velocity field,
        geometry, and density profile.

    Arguments:
        :w0 (*float*): rest cental wavelength of the line
        :fosc (*float*): oscilator strength of the line

    Keyword arguments:
            :VelocityField (*str*): radial velocity field indicating the assumed
                    treatment of the underlying physics and acceleration.
                    Options include
                         \'BetaCAK\' (1-1/x)^b from CAK theory
                         \'AccPlaw\' sqrt(1-x^(1-b)) from Steidel et al. 2010
                         \'VelPlaw\' 0.5(x-1)^b from various sources
            :Geometry (*str*): azimuthal geometry capturing the (an)isotropy
                    Options include
                        \'Sphere\' or \'Spherical\' isotropic case
                        \'FilledCones\' bidirectional cones
                        \'HollowCones\' or \'OpenCones\' cones with a cavity
                        \'HollowConesFixedCavity\' cones with fixed cavity size
            :DensityProfile (*str*): radial profile to use for the gas density.
                    Options include integration of many pulses or episodes
                         \'PowerLaw\'
                         \'Exponential\'
                         \'PowerLaw2\'
                     and individual pulse-like episodes or outbursts
                         \'LogNormal\'
                         \'Normal\'
                         \'Shell\'
                         \'FRED\'
                    or a combination of the two
                         \'Pulses\'
                         \'DampedPulses\'
            :Disk     (*bool*): boolean indicating whether a disk is present and
                    obstructs the cone posterior to the observer -- requires cone
                    geometry in order to take effect
            :Aperture (*bool*): boolean indicating whether the observed profile is
                    truncated due to the effects of a spherical aperture

    Returns:
        :Profile (*object*): An instance of the Profile object class inhereting
                the Absorption class with properties specific to an absorption line

    Attributes:
        :Profile (*function*): function returning a line profile based on the chosen
                velocity field, geometry, and density profile; takes wavelength and
                paramaters (see the `params` attribute)
        :params (*dict*): dictionary of parameters, initially assigned or later
                updated by the user
        :npar (*int*): number of assigned parameters
        :bounds (*dict*): dictionary of parameter bounds, initially assigned or
                later updated by the user

    Bound methods:
        :get_params: returns the current parameters
        :update_params: updates parameter(s) given the name(s) and value(s)
        :get_bounds: returns the current boundaries set for the parameters
        :update_bounds: updates bound(s) given the name(s) and value(s)
        :print_params: prints current parameter values to the command line

    '''
    def __init__(self,w0,fosc,pline):
        kwargs = {k:self.settings[k] for k in \
            ['VelocityField','DensityProfile','Geometry','Disk','FromRest','Aperture']}
        self.OutLinesModelA = __absorption__.build_profile_model(**kwargs)
        self.OutLinesModelR = __resonfluor__.build_profile_model(**kwargs)
        pass
    # define profile with and without a static ISM component
    def __WithStatic__(self,w,sigv,vinf,beta,*pars):
        OutLineA = array(list(map(lambda w0: \
          self.OutLinesModelA(w,w0,vinf,beta,*pars[2*self.nLines:]),self.w0)))
        OutLineR = array(list(map(lambda w0: \
          self.OutLinesModelR(w,w0,vinf,beta,*pars[2*self.nLines:]),self.w0)))
        staticv  = array(list(map(lambda w0: voigt(w,w0,sigv),self.w0)))
        #staticg  = array(list(map(lambda w0: gauss(w,w0,sigv),self.w0)))
        # classical absorption cross-section
        # log(sigma) = -14.8247 in cm^2, b in km/s, w in Ang
        # speed of light
        # c = 2.99792458e5 km/s
        logNs = array(pars[:self.nLines])
        logNo = array(pars[self.nLines:2*self.nLines])
        # static
        tausr = self.w0*self.fosc*10**(logNs-14.8247)/(vinf*2.99792458e5)
        fescs = array(list(map(self.__EscapeFraction__,tausr,self.pline,staticv)))
        srcfs = array(list(map(self.__SourceFunction__,tausr,staticv)))
        # OutLine
        tauor = self.w0*self.fosc*10**(logNo-14.8247)/(vinf*2.99792458e5)
        fesco = array(list(map(self.__EscapeFraction__,tauor,self.pline,OutLineR)))
        srcfo = array(list(map(self.__SourceFunction__,tauor,OutLineR)))
        return (fesco*srcfo).sum(axis=0) + (fescs*srcfs).sum(axis=0) + \
                exp( - tausr@staticv - tauor@OutLineA )
    def __WithoutStatic__(self,w,vinf,beta,*pars):
        OutLineA = array(list(map(lambda w0: \
          self.OutLinesModelA(w,w0,vinf,beta,*pars[2*self.nLines:]),self.w0)))
        OutLineR = array(list(map(lambda w0: \
          self.OutLinesModelR(w,w0,vinf,beta,*pars[2*self.nLines:]),self.w0)))
        # classical absorption cross-section
        # log(sigma) = -14.8247 in cm^2, b in km/s, w in Ang
        # speed of light c = 2.99792458e5 km/s
        logNo = array(pars[:self.nLines])
        tauor = self.w0*self.fosc*10**(logNo-14.8247)/(vinf*2.99792458e5)
        fesco = array(list(map(self.__EscapeFraction__,tauor,self.pline,OutLineR)))
        srcfo = array(list(map(self.__SourceFunction__,tauor,OutLineR)))
        return (fesco*srcfo).sum(axis=0) + exp( - tauor@OutLineA )
    def __OnlyStatic__(self,w,sigv,*pars):
        staticv  = array(list(map(lambda w0: voigt(w,w0,sigv),self.w0)))
        #staticg  = array(list(map(lambda w0: gauss(w,w0,sigv),self.w0)))
        # classical absorption cross-section
        # log(sigma) = -14.8247 in cm^2, b in km/s, w in Ang
        # speed of light
        # c = 2.99792458e5 km/s
        logNs = array(pars[:self.nLines])
        tausr = self.w0*self.fosc*10**(logNs-14.8247)/(vinf*2.99792458e5)
        fescs = array(list(map(self.__EscapeFraction__,tausr,self.pline,staticv)))
        srcfs = array(list(map(self.__SourceFunction__,tausr,staticv)))
        return (fescs*srcfs).sum(axis=0) + exp( - tausr@staticv)
    @staticmethod
    def __EscapeFraction__(tau,pline,profile):
        pesc = (1-exp(-tau*profile))/(tau*profile)
        fesc = pline*pesc/(1-pline*(1-pesc))
        fesc[profile<2**-24] = pline
        return fesc
    @staticmethod
    def __SourceFunction__(tau,profile):
        return 1-exp(-tau*profile)
