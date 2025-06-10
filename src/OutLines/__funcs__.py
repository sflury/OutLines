from numpy import trapz,array,append,ones,zeros,absolute,where,max,min,sum,\
                    log,log10,exp,sqrt,cos,sin,arccos,pi,inf,nan,isnan,arange
from scipy.optimize import brentq
from scipy.special import wofz
# possible density profiles
def dens_plaw(u,beta,VF,alpha):
    return x[VF](u,beta)**-alpha
def dens_lognorm(u,beta,VF,x1,sigma):
    return exp( -0.5*( log10(x[VF](u,beta)) - x1 )**2/sigma**2 )
def dens_norm(u,beta,VF,xi,sigma):
    return exp( -0.5*( x[VF](u,beta) - xi )**2/sigma**2 )
def dens_pulses(u,beta,VF,sigma,dx,x0):
    return sum([dens_norm(u,beta,VF,x+x0,sigma) for x in arange(0,11*dx+x0,dx)],axis=0)
def dens_pulsedamp(u,beta,VF,r,sigma,dx,x0):
    return dens_pulses(u,beta,VF,sigma,dx,x0) * dens_exp(u,beta,VF,r)
def dens_shell(u,beta,VF,xi,sigma):
    if hasattr(u,'__len__'):
        den = zeros(len(u))
        den[ absolute( x[VF](u,beta) - xi ) < sigma ] = 1
        return den
    else:
        if x[VF](u,beta) >  xi - sigma and x[VF](u,beta) < xi + sigma :
            return 1
        else:
            return 0
def dens_dplaw(u,beta,VF,a1,a2,x1):
    if hasattr(u,'__len__'):
        den = zeros(len(u))
        den[x[VF](u,beta) < x1]  = ( 1-u[x[VF](u,beta) < x1]**(1/beta) )**a1
        den[x[VF](u,beta) >= x1] = x1**(a2-a1) * ( 1-u[x[VF](u,beta) >= x1]**(1/beta) )**a2
        return den
    else:
        if x[VF](u,beta) < x1 :
            return ( 1-u**(1/beta) )**a1
        else:
            return x1**(a2-a1) * ( 1-u**(1/beta) )**a2
def dens_exp(u,beta,VF,r):
    return exp(r-r*x[VF](u,beta))
def dens_fred(u,beta,VF,r1,r2,x1):
    dx = x[VF](u,beta) - (x1-sqrt(r1*r2))
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
def x_cak(u,beta):
    return ( 1-u**(1/beta) )**-1
# Scarlata 2015 velocity power law
def x_vplaw(u,beta,A=0.03):
    return (u/A)**(1/beta)+1
# Steidel 2010 acceleration power law
def x_aplaw(u,beta):
    return ( 1  -  u**2 )**(1/(1-beta))
# possible velocity fields expressed as w = f(x)
# Pauldrach 1986 beta-law approximation to CAK 1975
def w_cak(xv,beta):
    return (1-1/xv)**beta
# Scarlata 2015 velocity power law
def w_vplaw(xv,beta,A=0.03):
    return A*(xv-1)**beta
# Steidel 2010 acceleration power law
def w_aplaw(xv,beta):
    return (1-xv**(1-beta))**0.5

# normalized density and velocity profiles to scale emissivity
global n,x,v
n = {'PowerLaw':       dens_plaw,\
     'LogNormal':      dens_lognorm,\
     'Normal':         dens_norm,\
     'Shell':          dens_shell,\
     'FRED':           dens_fred,\
     'PowerLaw2':      dens_dplaw,\
     'Exponential':    dens_exp,\
     'Pulses':         dens_pulses,\
     'DampedPulses':   dens_pulsedamp}
x = {'VelPlaw':   x_vplaw,\
     'AccPlaw':   x_aplaw,\
     'BetaCAK':   x_cak}
v = {'VelPlaw':   w_vplaw,\
     'AccPlaw':   w_aplaw,\
     'BetaCAK':   w_cak}
# Gaussian profile
def gauss(w,w0,sigv):
    return exp(-0.5*((w-w0)/(w0*sigv))**2)/sqrt(2*pi*(w0*sigv)**2)
# Voigt profile
def voigt(w,w0,sigv):
    return wofz((w-w0)/(w0*sigv)).real
# derivative of a Gaussian
def dog(w,w0,var=1e-7):
    return -w/var * gauss(w,w0,var)
# inverse of the radial velocity gradient via central finite difference method
def dxdw(w,beta,VF,h=2**-10):
    return (x[VF](w+h,beta)-x[VF](w,beta))/h
# difference between observed and predicted velocity
def diff_u(u,beta,yi,VF):
    return yi - u/x[VF](u,beta)*sqrt(x[VF](u,beta)**2-1)
# solve for minimum possible velocity
def solve_u(beta,VF,yi):
    # check if a sign change occurs
    if diff_u(yi,beta,yi,VF) *  diff_u(1-2**-24,beta,yi,VF) < 0 :
        # if so, root-find to determine the minimum velocity
        return brentq(diff_u,yi,1-2**-24,maxiter=100,args=(beta,yi,VF))
    # if not, field is close to terminal velocity and profile is ~ zero-valued
    else:
        return 1-2**-24
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
        # coefficients for band contour -- ellipse intersection
        A = (g/h)**2-1
        B = 2*S*(g/h)**2
        C = ((S/h)**2-1)*g**2
        # store
        geo += [S,g,h,A,B,C,incl+theta]
        # if a disk term included, calculate radius and projection
        # in velocity space for comparison with band contours
        if vdisk > 2**-10 :
            geo += [vdisk,vdisk*sin(incl),sin(incl)**-2-1]
        else:
            geo += [0,0,1]
    return array(geo)
# compute angular length of band contour enclosed by bicone
# S = displacement of ellipse
# g = axis of ellipse in x coordinate
# h = axis of ellipse in y coordinate
# A,B,C = coefficients for intersection of band contour with bicone ellipse
# rdisk, adisk = radius of disk
# adisk = sin(i) projection coefficient for intersection
def cone_bands(u,S,g,h,A,B,C,it,rdisk,rdisk_sini,adisk,cone,w):
    # if w cannot project onto u
    # then return no circle
    if u > w :
        return 0.
    # radius of band contour
    rband = sqrt(1-(u/w)**2)
    # intersection points (p,q) of band contour with projected spherical cap
    C += rband**2
    qc = (B+array([-1,1])*sqrt(B**2-4*A*C))/(2*A)
    # intersection points (pd,qd) of band contour with projected disk
    qd = sqrt((rband**2-rdisk**2)/adisk)
    # pc = sqrt(rbnd**2-qc**2) # --> no need to calculate
    # if no intersection occurs, qc[0] (always smaller than qc[1]) > rband
    # or b^2 < 4ac, giving a NaN intersection
    # if u/w close to 1, then band length -> 0 due to integral limit
    if qc[0] > rband or B**2 < 4*A*C or u/w > 1-2**-12 :
        return 0.
    # if disk present and obstructs all of the band
    elif rdisk_sini > rband and cone == 'post' :
        return 0.
    # if band contour is fully inscribed by the ellipse
    elif rband < h-S :
        # if dis is too small to obstruct band
        # or if cone is anterior to the disk, band is fully seen
        if rdisk < rband or cone == 'ante':
            return 1.
        # if disk present and large enough to obstruct some of the band,
        # and band is fully enclosed by the cone
        elif rdisk > rband and rdisk_sini < rband :
            return 2*arccos(qd/rband)/pi
    # if disk present and large enough to obstruct some of the band,
    # consider cones minus the band length within cone
    # which is obstructed by disk using either
    # the first intersection (posterior -- cone or cone + lune obstructed)
    # or the second intersection (anterior -- lune is obstructed)
    elif rdisk > rband and rdisk > 2**-10 :
        # the cone is posterior to disk
        if cone == 'post' :
            dl = 0.
            # disk is below projected cone -- no obstruction
            if qc[0] > 0 and qc[0] > qd :
                dl += arccos(qc[0]/rband)/pi
            # disk obstructs lower part of projected cone
            elif (qc[0] > 0 and qc[0] < qd) or (qc[0] < 0 and qc[0] > -qd) :
                dl += arccos(qd/rband)/pi
            # disk obstructs middle part of projected cone
            elif qc[0] < 0 and qc[0] < -qd :
                dl += (2*arccos(qd/rband)-arccos(qc[0]/rband))/pi
            # projected lune from anterior cone -- unaffected by disk
            if qc[1] <= rband and it > pi/2 :
                dl += arccos(qc[1]/rband)/pi
        # if cone is anterior to disk
        elif cone == 'ante' :
            dl = 0.
            # if first intersection below x-axis,
            # include full band minus the angle subtended by (pc,qc)
            if qc[0] < 0 and qc[0] > -rband :
                dl += 1-arccos(qc[0]/rband)/pi
            # if first intersection above x-axis,
            # include angle subtended by (pc,qc)
            else :
                dl += arccos(qc[0]/rband)/pi
            # if second intersection occurs and inclination + theta > pi/2,
            # there is an elliptical lune projected by the posterior cone
            # and is the arc length of the angle subtended by (pc,qc)
            # but can be obstructed by disk
            if qc[1] < rband and it > pi/2 :
                # no obstruction if lune starts above the disk
                if qc[1] > qd or rdisk < 2**-9 :
                    dl += arccos(qc[1]/rband)/pi
                # otherwise, disk determines the arc length
                elif qc[1] < qd and qd < rband :
                    dl += arccos(qd/rband)/pi
    # if no disk or if disk radius < band radius, just consider the cones
    else:
        dl = 0.
        # if first intersection below x-axis,
        # include full band minus the angle subtended by (pc,qc)
        if qc[0] < 0 and qc[0] > -rband :
            dl += 1-arccos(qc[0]/rband)/pi
        # if first intersection above x-axis,
        # include angle subtended by (pc,qc)
        elif qc[0] < rband :
            dl += arccos(qc[0]/rband)/pi
        # if second intersection occurs,
        # this is the elliptical lune projected by the other cone
        # and is the arc length of the angle subtended by (pc,qc)
        if qc[1] <= rband and it > pi/2 :
            dl += arccos(qc[1]/rband)/pi
    return dl
