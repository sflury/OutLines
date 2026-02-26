from numpy import array,append,ones,zeros,absolute,where,max,min,sum,cumsum,nansum,\
                    log,log10,exp,sqrt,cos,sin,arccos,arcsin,pi,inf,nan,isnan,\
                    isfinite,sign,arange,linspace,logspace,argmin,nanmin,\
                    argmax,nanmax,interp,diff,quantile,nanquantile,square,genfromtxt
from numpy.random import seed,default_rng
global uniform_sampler,normal_sampler
uniform_sampler = default_rng().random
normal_sampler  = default_rng().standard_normal
from scipy.optimize import brentq,newton,fminbound
from scipy.special import wofz
from functools import partial
###
### Density Profiles
###
##
## WINDS
##
# power law
def dens_plaw(u,beta,vini,VF,alpha):
    return x[VF](u,beta,vini)**-alpha
# exponential
def dens_exp(u,beta,vini,VF,r):
    return exp(r-r*x[VF](u,beta,vini))
# double power law
def dens_dplaw(u,beta,vini,VF,a1,a2,x1):
    xuv = x[VF](u,beta,vini)
    if hasattr(u,'__len__'):
        den = zeros(len(u))
        den[xuv < x1]  = xuv[xuv < x1]**-a1
        den[xuv >= x1] = x1**(a2-a1) * xuv[xuv >= x1]**-a2
        return den
    else:
        if xuv < x1 :
            return xuv**-a1
        else:
            return x1**(a2-a1) * xuv**-a2
##
## BUBBLES
##
# derivative of logistic column density
def dens_dlogic(u,beta,vini,VF,x0,k):
    expon = exp(-k*(x[VF](u,beta,vini)-x0))
    return k*expon/(1+expon)**2
# uniform shell
def dens_shell(u,beta,vini,VF,xi,sigma):
    if hasattr(u,'__len__'):
        den = zeros(len(u))
        den[ absolute( x[VF](u,beta,vini) - xi ) < sigma ] = 1
        return den
    else:
        if x[VF](u,beta,vini) >  xi - sigma and x[VF](u,beta,vini) < xi + sigma :
            return 1
        else:
            return 0
# fast rise, exponential decay
def dens_fred(u,beta,vini,VF,r1,r2,x1):
    dx = x[VF](u,beta,vini) - (x1-sqrt(r1*r2))
    if hasattr(u,'__len__'):
        den = zeros(len(u))
        den[dx>0] = exp(2*sqrt(r1/r2)) * exp(-r1/dx[dx>0]-dx[dx>0]/r2)
        return den
    else:
        if dx > 0:
            return exp(2*sqrt(r1/r2)) * exp(-r1/dx-dx/r2)
        else:
            return 0
# log normal
def dens_lognorm(u,beta,vini,VF,x1,sigma):
    return exp( -0.5*( log10(x[VF](u,beta,vini)) - x1 )**2/sigma**2 )
# normal
def dens_norm(u,beta,vini,VF,xi,sigma):
    return exp( -0.5*( x[VF](u,beta,vini) - xi )**2/sigma**2 )
##
## PULSE/SHELL ENSEMBLES
##
# ensemble with linear spacing
def dens_pulseslin(u,beta,vini,VF,sigma,dx,x0,pulse='Normal'):
    dens_pulse = n[pulse]
    return sum([dens_pulse(u,beta,vini,VF,x,sigma) for x in arange(x0,x[VF](1-2**-16,beta,vini)+dx,dx)],axis=0)
# ensemble with log spacing
def dens_pulseslog(u,beta,vini,VF,sigma,dx,x0,pulse='Normal'):
    dens_pulse = n[pulse]
    return sum([dens_pulse(u,beta,vini,VF,x,sigma) for x in 10**arange(log10(x0),log10(x[VF](1-2**-16,beta,vini))+dx,dx)],axis=0)
# damped ensemble with linear spacing
def dens_pulsedamplin(u,beta,vini,VF,r,sigma,dx,x0,pulse='Normal'):
    return dens_pulseslin(u,beta,vini,VF,sigma,dx,x0,pulse=pulse) * dens_exp(u,beta,vini,VF,r)
# damped ensemble with log spacing
def dens_pulsedamplog(u,beta,vini,VF,r,sigma,dx,x0,pulse='Normal'):
    return dens_pulseslog(u,beta,vini,VF,sigma,dx,x0,pulse=pulse) * dens_exp(u,beta,vini,VF,r)
# packet ensemble with linear spacing
def dens_pulsegrouplin(u,beta,vini,VF,s1,x1,sigma,dx,x0,pulse='Normal'):
    return dens_pulseslin(u,beta,vini,VF,sigma,dx,x0,pulse=pulse) * dens_norm(u,beta,vini,VF,x1,s1)
# packet ensemble with log spacing
def dens_pulsegrouplog(u,beta,vini,VF,s1,x1,sigma,dx,x0,pulse='Normal'):
    return dens_pulseslog(u,beta,vini,VF,sigma,dx,x0,pulse=pulse) * dens_norm(u,beta,vini,VF,x1,s1)
###
### Velocity Fields
###
##
## expressed as x = f(w)
##
# Castor-Lamers 1979 / Pauldrach 1986 beta-law approximation to CAK 1975
def x_cak(u,beta,vini):
    return ( 1-((u-vini)/(1-vini))**(1/beta) )**-1
#  velocity power law
def x_vplaw(u,beta,vini,A=0.5):
    return (((u-vini)/(1-vini))/A)**(1/beta)+1
# Steidel 2010 acceleration power law
def x_aplaw(u,beta,vini):
    return ( 1  -  ((u-vini)/(1-vini))**2 )**(1/(1-beta))
# Murray 2005 optically thick radiation pressure
def x_M2005(u,beta,vini):
    return exp( ((u-vini)/(1-vini))**2 )
# my own exponential law
def x_expon(u,beta,vini):
    return 1-log(1-(u-vini)/(1-vini))/beta
##
## expressed as w = f(x)
##
# Castor-Lamers 1979 / Pauldrach 1986 beta-law approximation to CAK 1975
def w_cak(xv,beta,vini):
    return (1-vini)*(1-1/xv)**beta + vini
#  velocity power law
def w_vplaw(xv,beta,vini,A=0.5):
    return A*(1-vini)*(xv-1)**beta + vini
# Steidel 2010 acceleration power law
def w_aplaw(xv,beta,vini):
    return (1-vini)*(1-xv**(1-beta))**0.5 + vini
# Murray 2005 optically thick radiation pressure
def w_M2005(xv,beta,vini):
    return (1-vini)*sqrt(log(x)) + vini
# my own exponential law
def w_expon(xv,beta,vini):
    return (1-vini)*(1-exp(-beta*(xv-1))) + vini
##
## related differentials -- velocity gradients
##
# generalized callable
def dxdw(w,beta,vini,VF):
    return dxdv[VF](w,beta,vini)
# inverse of the radial velocity gradient
# via central finite difference method
def dxdw_cfd(w,beta,VF,h=2**-20):
    return (x[VF](w+h,beta)-x[VF](w,beta))/h
# explicit inverse of the radial velocity gradient
def dxdw_cak(u,beta,vini):
    #return x_cak(u,beta)**2 * u**(1/beta-1)/beta
    coef = ((u-vini)/(1-vini))**(1/beta)/(beta*(u-vini))
    return coef * x_cak(u,beta,vini)**2
def dxdw_aplaw(u,beta,vini):
    #return 2*u/(beta-1) * x_aplaw(u,beta)**beta
    coef = 2*(u-vini)/((beta-1)*(1-vini)**2)
    return coef * x_aplaw(u,beta,vini)**beta
def dxdw_vplaw(u,beta,vini):
    #return (x_vplaw(u,beta)-1)/(beta*u)
    return (x_vplaw(u,beta,vini)-1)/(beta*(u-vini))
def dxdw_M2005(xv,beta,vini):
    coef = 2*(u-vini)/((1-vini)**2)
    return x_M2005(xv,beta,vini)*coeff
def dxdw_expon(u,beta,vini):
    return (beta*(1-u))**-1
###
### Dictionaries of Possible Models
###
# normalized density profiles
global n,x,v
n = {'PowerLaw':        dens_plaw,\
     'PowerLaw2':       dens_dplaw,\
     'Exponential':     dens_exp,\
     'DLogistic':       dens_dlogic,\
     'LogNormal':       dens_lognorm,\
     'Normal':          dens_norm,\
     'Shell':           dens_shell,\
     'FRED':            dens_fred,\
     'Pulses':          dens_pulseslin,\
     'DampedPulses':    dens_pulsedamplin,\
     'PacketPulses':    dens_pulsegrouplin,\
     'PulsesLog':       dens_pulseslog,\
     'DampedPulsesLog': dens_pulsedamplog,\
     'PacketPulsesLog': dens_pulsegrouplog,\
     }
# normalized radial profiles
x = {'VelPlaw':   x_vplaw,\
     'AccPlaw':   x_aplaw,\
     'BetaCAK':   x_cak,\
     'Expontl':   x_expon,\
     'M2005TK':   x_M2005}
# normalized velocity fields
v = {'VelPlaw':   w_vplaw,\
     'AccPlaw':   w_aplaw,\
     'BetaCAK':   w_cak,\
     'Expontl':   w_expon,\
     'M2005TK':   w_M2005}
# normalized velocity gradients
dxdv = {'VelPlaw':   dxdw_vplaw,\
        'AccPlaw':   dxdw_aplaw,\
        'BetaCAK':   dxdw_cak,\
        'Expontl':   dxdw_expon,
        'M2005TK':   dxdw_M2005}
###
### convenience functions for static gas
###
# Gaussian profile
def gauss(w,w0,sigv):
    return exp(-0.5*((w-w0)/(w0*sigv))**2)/sqrt(2*pi*(w0*sigv)**2)
# Voigt profile
def voigt(w,w0,sigv):
    return wofz((w-w0)/(w0*sigv)).real
# derivative of a Gaussian
def dog(w,w0,var=1e-7):
    return -w/var * gauss(w,w0,var)
###
### Geometry Root-Finders
###
# difference between observed and predicted velocity for source
def diff_u0(u0,beta,vini,yi,VF):
    return yi - u0/x[VF](u0,beta,vini)*sqrt(x[VF](u0,beta,vini)**2-1)
# solve for velocity at source edge
def solve_u0(beta,vini,VF,yi):
    # check if a sign change occurs
    if diff_u0(yi,beta,vini,yi,VF) *  diff_u0(1-2**-24,beta,vini,yi,VF) < 0 :
        # if so, root-find to determine the minimum velocity
        return brentq(diff_u0,yi,1-2**-24,maxiter=100,args=(beta,vini,yi,VF))
    # if not, field is close to terminal velocity and profile is ~ zero-valued
    else:
        return 1-2**-24
# difference between observed and predicted velocity for aperture
def diff_u1(u1,beta,vini,yi,ya,VF):
    return yi - u1/x[VF](u1,beta,vini)*sqrt(x[VF](u1,beta,vini)**2-x[VF](ya,beta,vini)**2)
# solve for velocity at aperture edge
def solve_u1(ya,beta,vini,VF,yi):
    # check if a sign change occurs
    if diff_u1(ya,beta,vini,yi,ya,VF) *  diff_u1(1-2**-24,beta,vini,yi,ya,VF) < 0 :
        # if so, root-find to determine the minimum velocity
        return brentq(diff_u1,ya,1-2**-24,maxiter=100,args=(beta,vini,yi,ya,VF))
    # if not, field is close to aperture so integration stops at aperture
    else:
        return ya
###
### Numerical Integrators
###
# monte carlo integrator, default n_samp=10^4 gives 1% error since E~sqrt(N)
def int_mc(func,w0,w1,u,vinf,beta,vini,geo,cone,par,VF,DP,n_samp=1e4):
    intgrd = lambda w: func(w,u,vinf,beta,vini,geo,cone,par,VF,DP)
    fun_mx = brentq(lambda w: -intgrd(w),w0,w1,maxiter=100)
    w_mc   = uniform_sampler(int(n_samp))*(w1-w0)+w0
    fun_mc = uniform_sampler(int(n_samp))*(fun_mx*(1+2**-8))
    weight = (w1-w0)*fun_mx
    n_incl = count_nonzero(fun_mc<=intgrd(w_mc))
    return weight*n_incl/n_samp
# gauss-legendre quadrature
def fixed_quad(func,a,b,args=()):
    lgndr_y = (b-a)*lgndr_x + a
    return (b-a)/2. * nansum(lgndr_w*func(lgndr_y,*args))
# trapezoid integration
def trapezoid(y,x):
    return sum( (x[1:]-x[:-1])*(y[:-1]+y[1:])/2 )
# cumulative trapezoid integration
def cumulative_trapezoid(y,x):
    cumtrapz = zeros(len(x))
    cumtrapz[1:] = cumsum( (x[1:]-x[:-1])*(y[:-1]+y[1:])/2 )
    return cumtrapz
###
### Non-spherical Geometries
###
# pre-calculate geometry terms which do not depend on velocity
# to speed up the integral
def precalc_geometry(incl,tO,tC,vdisk):
    geo = []
    for theta in [tO,tC]:
        # spherical cap projected as ellipse located at S
        # with horizontal axis g and vertical axis h
        h = 1.
        if theta < 2**-20 :        g = 0.
        elif theta > pi/2-2**-20 : g = 1.
        else:
            g = sin(theta)
            h *= sin(theta)
        # no inclination
        if incl < 2**-20 : S = 0.
        # with inclination
        else:
            S  = sin(incl)*cos(theta)
            h *= cos(incl)
        # coefficients for deprojected velocity -- ellipse intersection
        A = (g/h)**2-1
        B = 2*S*(g/h)**2
        C = ((S/h)**2-1)*g**2
        # store
        geo += [S,g,h,A,B,C,incl,theta]
        # if a disk term included, calculate radius and projection
        # in velocity space for comparison with deprojected velocity
        if vdisk > 1-2**-20 : geo += [1]
        elif vdisk < 2**-20 : geo += [0]
        else:                 geo += [vdisk]
        if incl < 2**-20 :    geo += [0,1]
        else:                 geo += [vdisk*cos(incl),cos(incl)**-2-1]
    return array(geo)
# compute angular length of deprojected velocity inscribed by hemisphere
# rdisk = radius of disk
# rdisk_sini = disk with inclination projection
# adisk = sin(i) projection coefficient for intersection
def hemi_inscr(u,vinf,i,t,rdisk,rdisk_sini,adisk,cone,w):
    # if w cannot project onto u
    # then return no circle
    if u > w : return 0.
    # if not inclined -- face-on
    if i <= 2**-20 :
        if cone == 'post' : return 0
        if cone == 'ante' : return 1
    # if fully inclined -- edge-on
    if i >= pi/2 - 2**-20 :
        return 0.5
    # radius of deprojected velocity from cos varrho = u/w
    rdprv = sqrt(1-(u/w)**2)
    # if disk spans all of the deprojected velocity
    if rdisk_sini > rdprv :
        #print(u,w,rdprv,rdisk_sini)
        # posterior to disk fully obstructed
        if cone == 'post' : return 0.
        # anterior to disk fully inscribed
        if cone == 'ante' : return 1.
    # otherwise, disk obstructs some of the deprojected velocity
    # hemisphere disk radius always 1
    # intersection points (pd,qd) of deprojected velocity with projected disk
    qd = sqrt((1-rdprv**2)/adisk)
    if cone == 'post' : return arccos(qd/rdprv)/pi
    else : return 1-arccos(qd/rdprv)/pi
# compute angular length of deprojected velocity inscribed by cone
# S = displacement of ellipse
# g = axis of ellipse in x coordinate
# h = axis of ellipse in y coordinate
# A,B,C = coefficients for intersection of deprojected velocity with bicone ellipse
# rdisk, adisk = radius of disk
# adisk = sin(i) projection coefficient for intersection
def cone_inscr(u,vinf,S,g,h,A,B,C,i,t,rdisk,rdisk_sini,adisk,cone,w):
    # if w cannot project onto u
    # then return no circle
    if u > w : return 0.
    # radius sin varrho of deprojected velocity
    # from cos varrho = u/w and sin^2 + cos^2 = 1
    rdprv = sqrt(1-(u/w)**2)
    # intersection points (p,q) of deprojected velocity with projected spherical cap
    C += rdprv**2
    qc = (B+array([-1,1])*sqrt(B**2-4*A*C))/(2*A)
    # pc = sqrt(rbnd**2-qc**2) # --> no need to calculate
    # unit arc length
    dl = 0.
    # if no intersection occurs, qc[0] (always smaller than qc[1]) > rdprv
    # or b^2 < 4ac, giving a NaN intersection
    # if u/w close to 1, then arc length -> 0 due to integral limit
    if qc[0] > rdprv or B**2 < 4*A*C or u/w > 1-2**-20 : return 0.
    # if deprojected velocity is fully inscribed by the ellipse
    # deprojected velocity is fully seen
    elif rdprv < h-S : return 1.
    # otherwise, full cone treatment
    else:
        # if first intersection below x-axis,
        # include full deprojected velocity minus the angle subtended by (pc,qc)
        if qc[0] < 0 and qc[0] > -rdprv : dl += 1-arccos(qc[0]/rdprv)/pi
        # if first intersection above x-axis,
        # include angle subtended by (pc,qc)
        elif qc[0] < rdprv : dl += arccos(qc[0]/rdprv)/pi
        # if second intersection occurs,
        # this is the elliptical lune projected by the other cone
        # and is the arc length of the angle subtended by (pc,qc)
        if qc[1] <= rdprv and i+t > pi/2 : dl += arccos(qc[1]/rdprv)/pi
    return dl
# compute angular length of deprojected velocity inscribed by cone + disk
# S = displacement of ellipse
# g = axis of ellipse in x coordinate
# h = axis of ellipse in y coordinate
# A,B,C = coefficients for intersection of deprojected velocity with bicone ellipse
# rdisk, adisk = radius of disk
# adisk = sin(i) projection coefficient for intersection
def cndk_inscr(u,vinf,S,g,h,A,B,C,i,t,rdisk,rdisk_sini,adisk,cone,w):
    # if w cannot project onto u
    # then return no circle
    if u > w : return 0.
    # radius of deprojected velocity from cos varrho = u/w
    rdprv = sqrt(1-(u/w)**2)
    # intersection points (p,q) of deprojected velocity with projected spherical cap
    C += rdprv**2
    qc = (B+array([-1,1])*sqrt(B**2-4*A*C))/(2*A)
    # intersection points (pd,qd) of deprojected velocity with projected disk
    qd = sqrt((rdprv**2-rdisk**2)/adisk)
    # pc = sqrt(rbnd**2-qc**2) # --> no need to calculate
    # unit arc length
    dl = 0.
    # if no intersection occurs, qc[0] (always smaller than qc[1]) > rdprv
    # or b^2 < 4ac, giving a NaN intersection
    # if u/w close to 1, then arc length -> 0 due to integral limit
    if qc[0] > rdprv or B**2 < 4*A*C or u/w > 1-2**-20 : return 0.
    # if disk present and obstructs all of the deprojected velocity
    elif rdisk_sini > rdprv and cone == 'post' : return 0.
    # if deprojected velocity is fully inscribed by the ellipse
    elif rdprv < h-S :
        # if dis is too small to obstruct deprojected velocity
        # or if cone is anterior to the disk, deprojected velocity is fully seen
        if rdisk < rdprv or cone == 'ante': return 1.
        # if disk present and large enough to obstruct some of the deprojected velocity,
        # and deprojected velocity is fully enclosed by the cone
        elif rdisk > rdprv and rdisk_sini < rdprv : return 2*arccos(qd/rdprv)/pi
    # if disk present and large enough to obstruct some of the deprojected velocity,
    # consider cones minus the arc length within cone
    # which is obstructed by disk using either
    # the first intersection (posterior -- cone or cone + lune obstructed)
    # or the second intersection (anterior -- lune is obstructed)
    elif rdisk > rdprv and rdisk > 2**-20 :
        # the cone is posterior to disk
        if cone == 'post' :
            # disk is below projected cone -- no obstruction
            if qc[0] > 0 and qc[0] > qd : dl += arccos(qc[0]/rdprv)/pi
            # disk obstructs lower part of projected cone
            elif (qc[0] > 0 and qc[0] < qd) or (qc[0] < 0 and qc[0] > -qd) : dl += arccos(qd/rdprv)/pi
            # disk obstructs middle part of projected cone
            elif qc[0] < 0 and qc[0] < -qd : dl += (2*arccos(qd/rdprv)-arccos(qc[0]/rdprv))/pi
            # projected lune from anterior cone -- unaffected by disk
            if qc[1] <= rdprv and i+t > pi/2 : dl += arccos(qc[1]/rdprv)/pi
        # if cone is anterior to disk
        elif cone == 'ante' :
            # if first intersection below x-axis,
            # include full deprojected velocity minus the angle subtended by (pc,qc)
            if qc[0] < 0 and qc[0] > -rdprv : dl += 1-arccos(qc[0]/rdprv)/pi
            # if first intersection above x-axis,
            # include angle subtended by (pc,qc)
            else : dl += arccos(qc[0]/rdprv)/pi
            # if second intersection occurs and inclination + theta > pi/2,
            # there is an elliptical lune projected by the posterior cone
            # and is the arc length of the angle subtended by (pc,qc)
            # but can be obstructed by disk
            if qc[1] < rdprv and i+t > pi/2 :
                # no obstruction if lune starts above the disk
                if qc[1] > qd or rdisk < 2**-9 : dl += arccos(qc[1]/rdprv)/pi
                # otherwise, disk determines the arc length
                elif qc[1] < qd and qd < rdprv : dl += arccos(qd/rdprv)/pi
    # if no disk or if disk radius < deprojected velocity radius, just consider the cones
    else:
        # if first intersection below x-axis,
        # include full deprojected velocity minus the angle subtended by (pc,qc)
        if qc[0] < 0 and qc[0] > -rdprv : dl += 1-arccos(qc[0]/rdprv)/pi
        # if first intersection above x-axis,
        # include angle subtended by (pc,qc)
        elif qc[0] < rdprv : dl += arccos(qc[0]/rdprv)/pi
        # if second intersection occurs,
        # this is the elliptical lune projected by the other cone
        # and is the arc length of the angle subtended by (pc,qc)
        if qc[1] <= rdprv and i+t > pi/2 : dl += arccos(qc[1]/rdprv)/pi
    return dl

# gauss-legendre quadrature nodes and weights for n=96
global lgndr_x,lgrndr_w,ell_zeros
lgndr_x,lgndr_w = array([
    [0.0162767448496030,3.2550614492363350e-02],\
    [0.0488129851360497,3.2516118713869058e-02],\
    [0.0812974954644255,3.2447163714064503e-02],\
    [0.1136958501106659,3.2343822568576118e-02],\
    [0.1459737146548969,3.2206204794030469e-02],\
    [0.1780968823676186,3.2034456231992893e-02],\
    [0.2100313104605672,3.1828758894411134e-02],\
    [0.2417431561638400,3.1589330770727397e-02],\
    [0.2731988125910491,3.1316425596861541e-02],\
    [0.3043649443544964,3.1010332586313926e-02],\
    [0.3352085228926254,3.0671376123669377e-02],\
    [0.3656968614723136,3.0299915420827796e-02],\
    [0.3957976498289086,2.9896344136328527e-02],\
    [0.4254789884073006,2.9461089958167899e-02],\
    [0.4547094221677431,2.8994614150555195e-02],\
    [0.4834579739205964,2.8497411065085482e-02],\
    [0.5116941771546677,2.7970007616848484e-02],\
    [0.5393881083243575,2.7412962726029329e-02],\
    [0.5665104185613972,2.6826866725591932e-02],\
    [0.5930323647775720,2.6212340735672534e-02],\
    [0.6189258401254685,2.5570036005349562e-02],\
    [0.6441634037849671,2.4900633222483839e-02],\
    [0.6687183100439161,2.4204841792364991e-02],\
    [0.6925645366421715,2.3483399085926417e-02],\
    [0.7156768123489676,2.2737069658329195e-02],\
    [0.7380306437444001,2.1966644438744458e-02],\
    [0.7596023411766475,2.1172939892191479e-02],\
    [0.7803690438674331,2.0356797154333719e-02],\
    [0.8003087441391408,1.9519081140145365e-02],\
    [0.8194003107379316,1.8660679627411546e-02],\
    [0.8376235112281871,1.7782502316045324e-02],\
    [0.8549590334346014,1.6885479864245576e-02],\
    [0.8713885059092965,1.5970562902561790e-02],\
    [0.8868945174024205,1.5038721026994580e-02],\
    [0.9014606353158523,1.4090941772314505e-02],\
    [0.9150714231208981,1.3128229566961599e-02],\
    [0.9277124567223087,1.2151604671088017e-02],\
    [0.9393703397527551,1.1162102099839300e-02],\
    [0.9500327177844377,1.0160770535007201e-02],\
    [0.9596882914487425,9.1486712307841406e-03],\
    [0.9683268284632642,8.1268769256997199e-03],\
    [0.9759391745851365,7.0964707911539303e-03],\
    [0.9825172635630147,6.0585455042361323e-03],\
    [0.9880541263296237,5.0142027429296021e-03],\
    [0.9925439003237626,3.9645543384451793e-03],\
    [0.9959818429872094,2.9107318179334802e-03],\
    [0.9983643758631817,1.8539607889386321e-03],\
    [0.9996895038832307,7.9679206555383390e-04],\
    ]).T
# add in negatives -- symmetric weights
lgndr_x = 0.5*(array([-lgndr_x[::-1],lgndr_x]).flatten()+1.)
lgndr_w = array([ lgndr_w[::-1],lgndr_w]).flatten()
# reference array of zeros
ell_zeros = zeros(len(lgndr_x))
