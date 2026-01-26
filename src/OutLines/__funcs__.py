from numpy import array,append,ones,zeros,absolute,where,max,min,sum,nansum,\
                    log,log10,exp,sqrt,cos,sin,arccos,arcsin,pi,inf,nan,isnan,\
                    isfinite,sign,arange,linspace,logspace,argmin,nanmin,\
                    argmax,nanmax,interp,diff,quantile,nanquantile,count_nonzero
from numpy.random import seed,default_rng
global uniform_sampler
uniform_sampler = default_rng().random
from scipy.optimize import brentq,newton,fminbound
from scipy.integrate import trapezoid,simpson,fixed_quad,cumulative_trapezoid
from scipy.special import wofz
from functools import partial
# possible density profiles
def dens_plaw(u,beta,vini,VF,alpha):
    return x[VF](u,beta,vini)**-alpha
def dens_lognorm(u,beta,vini,VF,x1,sigma):
    return exp( -0.5*( log10(x[VF](u,beta,vini)) - x1 )**2/sigma**2 )
def dens_norm(u,beta,vini,VF,xi,sigma):
    return exp( -0.5*( x[VF](u,beta,vini) - xi )**2/sigma**2 )
def dens_pulseslin(u,beta,vini,VF,sigma,dx,x0):
    return sum([dens_norm(u,beta,vini,VF,x,sigma) for x in arange(x0,x[VF](0.95,beta,vini)+dx,dx)],axis=0)
def dens_pulseslog(u,beta,vini,VF,sigma,dx,x0):
    return sum([dens_norm(u,beta,vini,VF,x,sigma) for x in 10**arange(log10(x0),log10(x[VF](0.95,beta,vini))+dx,dx)],axis=0)
def dens_pulsedamplin(u,beta,vini,VF,r,sigma,dx,x0):
    return dens_pulseslin(u,beta,vini,VF,sigma,dx,x0) * dens_exp(u,beta,vini,VF,r)
def dens_pulsedamplog(u,beta,vini,VF,r,sigma,dx,x0):
    return dens_pulseslog(u,beta,vini,VF,sigma,dx,x0) * dens_exp(u,beta,vini,VF,r)
def dens_dlogic(u,beta,vini,VF,k,x0):
    expon = exp(-k*(x[VF](u,beta,vini)-x0))
    return k*expon/(1+expon)**2
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
def dens_exp(u,beta,vini,VF,r):
    return exp(r-r*x[VF](u,beta,vini))
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
# possible velocity fields expressed as x = f(w)
# Pauldrach 1986 beta-law approximation to CAK 1975
def x_cak(u,beta,vini):
    return ( 1-((u-vini)/(1-vini))**(1/beta) )**-1
# Scarlata 2015 velocity power law
def x_vplaw(u,beta,vini,A=0.5):
    return (((u-vini)/(1-vini))/A)**(1/beta)+1
# Steidel 2010 acceleration power law
def x_aplaw(u,beta,vini):
    return ( 1  -  ((u-vini)/(1-vini))**2 )**(1/(1-beta))
# possible velocity fields expressed as w = f(x)
# Pauldrach 1986 beta-law approximation to CAK 1975
def w_cak(xv,beta,vini):
    return (1-vini)*(1-1/xv)**beta + vini
# Scarlata 2015 velocity power law
def w_vplaw(xv,beta,vini,A=0.5):
    return A*(1-vini)*(xv-1)**beta + vini
# Steidel 2010 acceleration power law
def w_aplaw(xv,beta,vini):
    return (1-vini)*(1-xv**(1-beta))**0.5 + vini
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
# normalized density and velocity profiles to scale emissivity
global n,x,v
n = {'PowerLaw':       dens_plaw,\
     'PowerLaw2':      dens_dplaw,\
     'Exponential':    dens_exp,\
     'DerivLogistic':  dens_dlogic,\
     'LogNormal':      dens_lognorm,\
     'Normal':         dens_norm,\
     'Shell':          dens_shell,\
     'FRED':           dens_fred,\
     'Pulses':         dens_pulseslin,\
     'DampedPulses':   dens_pulsedamplin,\
     'DampedPulsesLog':dens_pulsedamplog}
x = {'VelPlaw':   x_vplaw,\
     'AccPlaw':   x_aplaw,\
     'BetaCAK':   x_cak}
v = {'VelPlaw':   w_vplaw,\
     'AccPlaw':   w_aplaw,\
     'BetaCAK':   w_cak}
dxdv = {'VelPlaw':   dxdw_vplaw,\
        'AccPlaw':   dxdw_aplaw,\
        'BetaCAK':   dxdw_cak}
# Gaussian profile
def gauss(w,w0,sigv):
    return exp(-0.5*((w-w0)/(w0*sigv))**2)/sqrt(2*pi*(w0*sigv)**2)
# Voigt profile
def voigt(w,w0,sigv):
    return wofz((w-w0)/(w0*sigv)).real
# derivative of a Gaussian
def dog(w,w0,var=1e-7):
    return -w/var * gauss(w,w0,var)
# derivative of velocity fields
def dxdw(w,beta,vini,VF):
    return dxdv[VF](w,beta,vini)
# inverse of the radial velocity gradient via central finite difference method
#def dxdw(w,beta,VF,h=2**-10):
#    return (x[VF](w+h,beta)-x[VF](w,beta))/h
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
# monte carlo integrator, default n_samp=10^4 gives 1% error since E~sqrt(N)
def int_mc(func,w0,w1,n_samp=1e4):
    fun_mx = brentq(lambda w: -func(w),w0,w1,maxiter=100)
    w_mc   = uniform_sampler(n_samp)*(w1-w0)+w0
    fun_mc = uniform_sampler(n_samp)*(f_mx+2**-24)
    mc_rng = (w1-w0)*f_mx
    n_incl = count_nonzero(f_mc<=func(w_mc))
    return mc_rng*incl/n_samp
# pre-calculate geometry terms which do not depend on velocity
# to speed up the integral
def precalc_geometry(incl,tO,tC,vdisk):
    geo = []
    for theta in [tO,tC]:
        # spherical cap projected as ellipse located at f
        # with horizontal axis g and vertical axis h
        S = sin(incl)*cos(theta)
        g = sin(theta)
        h = g*cos(incl)
        # coefficients for deprojected velocity -- ellipse intersection
        if h == 0 :
            A = 1
            B = 0
            C = 0
        else:
            A = (g/h)**2-1
            B = 2*S*(g/h)**2
            C = ((S/h)**2-1)*g**2
        # store
        geo += [S,g,h,A,B,C,incl,theta]
        # if a disk term included, calculate radius and projection
        # in velocity space for comparison with deprojected velocity
        if vdisk > 2**-10 :
            geo += [vdisk,vdisk*sin(incl),sin(incl)**-2-1]
        else:
            geo += [0,0,1]
    return array(geo)
# compute angular length of deprojected velocity inscribed by cone
# S = displacement of ellipse
# g = axis of ellipse in x coordinate
# h = axis of ellipse in y coordinate
# A,B,C = coefficients for intersection of deprojected velocity with bicone ellipse
# rdisk, adisk = radius of disk
# adisk = sin(i) projection coefficient for intersection
def cone_inscr(u,S,g,h,A,B,C,i,t,rdisk,rdisk_sini,adisk,cone,w):
    # if w cannot project onto u
    # then return no circle
    if u > w :
        return 0.
    # radius of deprojected velocity
    rdprv = sqrt(1-(u/w)**2)
    # intersection points (p,q) of deprojected velocity with projected spherical cap
    C += rdprv**2
    qc = (B+array([-1,1])*sqrt(B**2-4*A*C))/(2*A)
    # intersection points (pd,qd) of deprojected velocity with projected disk
    qd = sqrt((rdprv**2-rdisk**2)/adisk)
    # pc = sqrt(rbnd**2-qc**2) # --> no need to calculate
    # if no intersection occurs, qc[0] (always smaller than qc[1]) > rdprv
    # or b^2 < 4ac, giving a NaN intersection
    # if u/w close to 1, then arc length -> 0 due to integral limit
    if qc[0] > rdprv or B**2 < 4*A*C or u/w > 1-2**-12 :
        return 0.
    # if disk present and obstructs all of the deprojected velocity
    elif rdisk_sini > rdprv and cone == 'post' :
        return 0.
    # if deprojected velocity is fully inscribed by the ellipse
    elif rdprv < h-S :
        # if dis is too small to obstruct deprojected velocity
        # or if cone is anterior to the disk, deprojected velocity is fully seen
        if rdisk < rdprv or cone == 'ante':
            return 1.
        # if disk present and large enough to obstruct some of the deprojected velocity,
        # and deprojected velocity is fully enclosed by the cone
        elif rdisk > rdprv and rdisk_sini < rdprv :
            return 2*arccos(qd/rdprv)/pi
    # if disk present and large enough to obstruct some of the deprojected velocity,
    # consider cones minus the arc length within cone
    # which is obstructed by disk using either
    # the first intersection (posterior -- cone or cone + lune obstructed)
    # or the second intersection (anterior -- lune is obstructed)
    elif rdisk > rdprv and rdisk > 2**-10 :
        # the cone is posterior to disk
        if cone == 'post' :
            dl = 0.
            # disk is below projected cone -- no obstruction
            if qc[0] > 0 and qc[0] > qd :
                dl += arccos(qc[0]/rdprv)/pi
            # disk obstructs lower part of projected cone
            elif (qc[0] > 0 and qc[0] < qd) or (qc[0] < 0 and qc[0] > -qd) :
                dl += arccos(qd/rdprv)/pi
            # disk obstructs middle part of projected cone
            elif qc[0] < 0 and qc[0] < -qd :
                dl += (2*arccos(qd/rdprv)-arccos(qc[0]/rdprv))/pi
            # projected lune from anterior cone -- unaffected by disk
            if qc[1] <= rdprv and i+t > pi/2 :
                dl += arccos(qc[1]/rdprv)/pi
        # if cone is anterior to disk
        elif cone == 'ante' :
            dl = 0.
            # if first intersection below x-axis,
            # include full deprojected velocity minus the angle subtended by (pc,qc)
            if qc[0] < 0 and qc[0] > -rdprv :
                dl += 1-arccos(qc[0]/rdprv)/pi
            # if first intersection above x-axis,
            # include angle subtended by (pc,qc)
            else :
                dl += arccos(qc[0]/rdprv)/pi
            # if second intersection occurs and inclination + theta > pi/2,
            # there is an elliptical lune projected by the posterior cone
            # and is the arc length of the angle subtended by (pc,qc)
            # but can be obstructed by disk
            if qc[1] < rdprv and i+t > pi/2 :
                # no obstruction if lune starts above the disk
                if qc[1] > qd or rdisk < 2**-9 :
                    dl += arccos(qc[1]/rdprv)/pi
                # otherwise, disk determines the arc length
                elif qc[1] < qd and qd < rdprv :
                    dl += arccos(qd/rdprv)/pi
    # if no disk or if disk radius < deprojected velocity radius, just consider the cones
    else:
        dl = 0.
        # if first intersection below x-axis,
        # include full deprojected velocity minus the angle subtended by (pc,qc)
        if qc[0] < 0 and qc[0] > -rdprv :
            dl += 1-arccos(qc[0]/rdprv)/pi
        # if first intersection above x-axis,
        # include angle subtended by (pc,qc)
        elif qc[0] < rdprv :
            dl += arccos(qc[0]/rdprv)/pi
        # if second intersection occurs,
        # this is the elliptical lune projected by the other cone
        # and is the arc length of the angle subtended by (pc,qc)
        if qc[1] <= rdprv and i+t > pi/2 :
            dl += arccos(qc[1]/rdprv)/pi
    return dl
