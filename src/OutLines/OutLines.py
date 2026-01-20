from numpy import trapezoid,log,log10,exp,sqrt,array,append,ones,zeros,absolute
from scipy.integrate import romberg
from scipy.optimize import brentq
'''
Name:
    tau_out

Purpose:
    Calculate the absorption line profile for a galactic outflow
    assuming a power-law density field, a beta-law velocity field,
    and the Sobolev approximation.

Arguments:
    :wave  (*np.ndarray*) : array of observed wavelengths
    :wave0 (*float*) : rest-frame central wavelength
    :vinf  (*float*) : terminal velocity in c units
    :alpha (*float*) : power-law index of density field
    :beta  (*float*) : power-law index of velocity field

Returns:
    :tau (*np.ndarray*) : unitless absorption line profile
'''
def tau_out(wave,wave0,vinf,alpha,beta):
    # observed velocity from wavelengths
    u = absolute(wave-wave0)/wave0 / vinf
    # normalized radius from velocity
    x = lambda w: ( 1-w**(1/beta) )**-1 if 1-w**(1/beta) != 0 else 1e20
    # normalized density profile to scale optical depth
    n = lambda w: ( 1-w**(1/beta) )**alpha
    # inverse of the radial velocity gradient
    dxdw = lambda w: w**(1/beta-1) / ( beta * (1-w**(1/beta))**2 )
    # max possible velocity of gas contributing to profile
    # constrained by arbitrary radius of the background continuum source
    wmin = array(list( map( lambda ui: \
        brentq(lambda wl: ui - wl/x(wl)*sqrt(x(wl)**2-1),ui,1), \
        u[(u<1)&(wave<wave0)]) ))
    # absorption profile given by integrating over density / velocity gradient
    tau1 = zeros(len(u))
    tau1[(u<1)&(wave<wave0)] = \
        array(list(map(lambda wi,wm: romberg(lambda wj: n(wj)*dxdw(wj),wi,wm), \
        u[(u<1)&(wave<wave0)],wmin)))
    return tau1
'''
Name:
    phi_out

Purpose:
    Calculate the nebular emission line profile for a galactic outflow
    assuming a power-law density field, a beta-law velocity field,
    and the Sobolev approximation.

Arguments:
    :wave  (*np.ndarray*) : array of observed wavelengths
    :wave0 (*float*) : rest-frame central wavelength
    :vinf  (*float*) : terminal velocity in c units
    :alpha (*float*) : power-law index of density field
    :beta  (*float*) : power-law index of velocity field

Returns:
    :phi (*np.ndarray*) : emission line profile in units of wavelength^(-1)
'''
def phi_out(w,w0,vinf,alpha,beta):
    # observed velocity from wavelengths
    y = absolute(w-w0)/w0 /vinf
    # normalized radius from velocity
    x = lambda u: ( 1-u**(1/beta) )**-1
    # normalized density profile to scale emissivity
    n = lambda u: ( 1-u**(1/beta) )**alpha
    # u = v/vinf (not y since y here is used for vobs/v)
    # with no occultation, the minimum velocity is just the observed velocity
    umin = array(list( map(lambda yi: max([yi,0.0]), y[y<1]) ))
    # when occultation by source occurs, need to solve
    # for minimum contributing velocity to observed y
    umin[w[y<1]>w0] = array(list( map(lambda yi:
        brentq(lambda u: yi - u/x(u)*sqrt(x(u)**2-1),yi,1-1e-9,maxiter=100),
        y[(y<1)&(w>w0)]) ))
    # by definition, velocity cannot exceed terminal velocity
    # since u = v/vinf, max value is simply 1
    umax = 1.
    # profile integrated at each velocity
    phi = zeros(len(y))
    phi[y<1] = array( list( map(
        lambda u1: romberg(lambda u: n(u)**2,u1,umax),umin) ) )
    # normalize
    norm = trapezoid(phi,x=w)
    return phi/norm
