from OutLines.__funcs__ import *
# velocity integral limits
def calc_limits(w,w0,vinf,vini,vapr,beta,VF,Inflow=False):
    # error check for length of w
    if not hasattr(w,'__len__'):
        w = array([w])
    # observed velocity from wavelengths, in vinf/c units
    y = absolute(w-w0)/w0 / vinf
    ind = where(y<1)[0]
    # u = v/vinf (not y since y here is used for vobs/v)
    # with no occultation, the minimum velocity is just the observed velocity
    umin = y[ind]
    # primary cone type -- start with all as anterior
    cone = array(len(umin)*['ante']).astype(str)
    # when occultation by source occurs, need to solve
    # for minimum contributing velocity to observed y
    if Inflow :
        umin[w[ind]<w0] = array(list(map(partial(solve_u0,beta,vini,VF),y[ind][w[ind]<w0])))
        cone[w[ind]<w0] = 'post'
    else:
        umin[w[ind]>w0] = array(list(map(partial(solve_u0,beta,vini,VF),y[ind][w[ind]>w0])))
        cone[w[ind]>w0] = 'post'
    # by definition, velocity cannot exceed terminal velocity
    # since u = v/vinf, max value is simply 1 unless aperture effects occur
    umax = ones(len(ind))
    if vapr < vinf :
        wapr = array(list(map(partial(solve_u1,vapr/vinf,beta,vini,VF),y[ind])))
        umax = min([umax,wapr],axis=0)
    if vini > 0 :
        umin[umin<vini/vinf] = vini/vinf
    return umin,umax,cone
# resonant or fluorescent emission
def resfluor(w,u,beta,vini,geo,cone,par,VF,DP):
    if geo[1] == 1. :
        return n[DP](w,beta,vini,VF,*par)*dxdw(w,beta,vini,VF)
    else:
        ell  = array(list(map(partial(cone_bands,u,*geo[:11],cone),w)))
    if geo[12] > 0. :
        ell -= array(list(map(partial(cone_bands,u,*geo[11:],cone),w)))
    ell = nanmax([ell,zeros(256)],axis=0)
    return ell*n[DP](w,beta,vini,VF,*par)*dxdw(w,beta,vini,VF)
# integral over column densities for range of allowed velocities
def phi_int(beta,incl,tO,tC,vdisk,vini,par,VF,DP,u,umin,umax,cone):
    # if cone projection excludes some velocities, then limit the integral
    if incl+tO < pi/2 :
        umax = min([umax,u/cos(incl+tO)])
    if incl > tO :
        umin   = max([umin,u/cos(incl-tO)])
    # if minimum exceeds maximum due to geometry limits,
    # no velocity bands are contained
    if umin >= umax:
        return 0
    else:
        # set up geometry terms that do not depend on velocity
        geo = precalc_geometry(incl,tO,tC,vdisk)
        # return the integral
        return fixed_quad(resfluor,umin,umax,args=(u,beta,vini,geo,cone,par,VF,DP),n=256)[0]
# calculate unnormalized profile for a sphere or bicone
def calc_phi(w,w0,vinf,beta,incl,tO,tC,xdisk,vini,vapr,*par,VF='BetaCAK',DP='PowerLaw',Geometry='Sphere'):
    # obtain velocity limits for integral
    umin,umax,cone = calc_limits(w,w0,vinf,vini,vapr,beta,VF)
    u = absolute(w-w0)/w0 / vinf
    # check disk radius and convert to velocity
    if xdisk <= 1 :
        vdisk = 2**-10
    else:
        vdisk = max([v[VF](xdisk,beta,vini),2**-10])
    # profile integrated at each velocity
    phi = zeros(len(w))
    phi[ u < 1 ] = array( list( map( \
            partial(phi_int,beta,incl,tO,tC,vdisk,vini,par,VF,DP),\
            u[ u < 1 ], umin, umax, cone ) ) )
    return phi

'''
Name:
    build_profile_model

Purpose:
    Construct a unique function based on user-specified model options. This
    function will then compute and return the normalized line profiles.

Keyword Arguments:
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
    :profile (*function*): a function which calculates the line profile for a
                        galactic outflow for the user-provided model settings

'''
def build_profile_model(VelocityField='BetaCAK',DensityProfile='PowerLaw',Geometry='Sphere',Aperture=False,Disk=False,FromRest=True):
    kwargs = dict(VF=VelocityField,DP=DensityProfile)
    if FromRest:
        if 'spher' in Geometry.lower():
            if Aperture :
                def profile(w,w0,vinf,beta,vapr,*par):
                    phi = calc_phi(w,w0,vinf,beta,0,pi/2,0,0,0,vapr,*par,**kwargs)
                    return phi/trapezoid(phi,x=w)
            else:
                def profile(w,w0,vinf,beta,*par):
                    phi = calc_phi(w,w0,vinf,beta,0,pi/2,0,0,0,1,*par,**kwargs)
                    return phi/trapezoid(phi,x=w)
        elif 'filled' in Geometry.lower():
            if Aperture and Disk :
                def profile(w,w0,vinf,beta,incl,tO,xdisk,vapr,*par):
                    phi = calc_phi(w,w0,vinf,beta,incl,tO,0,xdisk,0,vapr,*par,**kwargs)
                    return phi/trapezoid(phi,x=w)
            elif Disk and not Aperture :
                def profile(w,w0,vinf,beta,incl,tO,xdisk,*par):
                    phi = calc_phi(w,w0,vinf,beta,incl,tO,0,xdisk,0,1,*par,**kwargs)
                    return phi/trapezoid(phi,x=w)
            elif Aperture and not Disk :
                def profile(w,w0,vinf,beta,incl,tO,vapr,*par):
                    phi = calc_phi(w,w0,vinf,incl,tO,0,0,0,vapr,*par,**kwargs)
                    return phi/trapezoid(phi,x=w)
            else:
                def profile(w,w0,vinf,beta,incl,tO,*par):
                    phi = calc_phi(w,w0,vinf,beta,incl,tO,0,0,0,1,*par,**kwargs)
                    return phi/trapezoid(phi,x=w)
        elif 'hollow' in Geometry.lower() or 'cavity' in Geometry.lower() or 'open' in Geometry.lower():
            if 'fix' in Geometry.lower():
                if Aperture and Disk :
                    def profile(w,w0,vinf,beta,incl,tO,xdisk,vapr,*par):
                        phi = calc_phi(w,w0,vinf,beta,incl,tO,tO-0.174533,xdisk,0,vapr,*par,**kwargs)
                        return phi/trapezoid(phi,x=w)
                elif Disk and not Aperture :
                    def profile(w,w0,vinf,beta,incl,tO,xdisk,*par):
                        phi = calc_phi(w,w0,vinf,beta,incl,tO,tO-0.174533,xdisk,0,1,*par,**kwargs)
                        return phi/trapezoid(phi,x=w)
                elif Aperture and not Disk :
                    def profile(w,w0,vinf,beta,incl,tO,vapr,*par):
                        phi = calc_phi(w,w0,vinf,beta,incl,tO,tO-0.174533,0,0,vapr,*par,**kwargs)
                        return phi/trapezoid(phi,x=w)
                else:
                    def profile(w,w0,vinf,beta,incl,tO,*par):
                        phi = calc_phi(w,w0,vinf,beta,incl,tO,tO-0.174533,0,0,1,*par,**kwargs)
                        return phi/trapezoid(phi,x=w)
            else:
                if Aperture and Disk :
                    def profile(w,w0,vinf,beta,incl,tO,tC,xdisk,vapr,*par):
                        phi = calc_phi(w,w0,vinf,beta,incl,tO,tC,xdisk,0,vapr,*par,**kwargs)
                        return phi/trapezoid(phi,x=w)
                elif Disk and not Aperture :
                    def profile(w,w0,vinf,beta,incl,tO,tC,xdisk,*par):
                        phi = calc_phi(w,w0,vinf,beta,incl,tO,tC,xdisk,0,1,*par,**kwargs)
                        return phi/trapezoid(phi,x=w)
                elif Aperture and not Disk :
                    def profile(w,w0,vinf,beta,incl,tO,tC,vapr,*par):
                        phi = calc_phi(w,w0,vinf,beta,incl,tO,tC,0,0,vapr,*par,**kwargs)
                        return phi/trapezoid(phi,x=w)
                else:
                    def profile(w,w0,vinf,beta,incl,tO,tC,*par):
                        phi = calc_phi(w,w0,vinf,beta,incl,tO,tC,0,0,1,*par,**kwargs)
                        return phi/trapezoid(phi,x=w)
        else:
            print('Geometry not recognized. Options are \'Sphere\',\n'+
                '\'FilledCones\', \'HollowCones\', or \'HollowConesFixedCavity\'.')
    else:
        if 'spher' in Geometry.lower():
            if Aperture :
                def profile(w,w0,vinf,beta,vini,vapr,*par):
                    phi = calc_phi(w,w0,vinf,beta,0,pi/2,0,0,vini,vapr,*par,**kwargs)
                    return phi/trapezoid(phi,x=w)
            else:
                def profile(w,w0,vinf,beta,vini,*par):
                    phi = calc_phi(w,w0,vinf,beta,0,pi/2,0,0,vini,1,*par,**kwargs)
                    return phi/trapezoid(phi,x=w)
        elif 'filled' in Geometry.lower():
            if Aperture and Disk :
                def profile(w,w0,vinf,beta,incl,tO,xdisk,vini,vapr,*par):
                    phi = calc_phi(w,w0,vinf,beta,incl,tO,0,xdisk,vini,vapr,*par,**kwargs)
                    return phi/trapezoid(phi,x=w)
            elif Disk and not Aperture :
                def profile(w,w0,vinf,beta,incl,tO,xdisk,vini,*par):
                    phi = calc_phi(w,w0,vinf,beta,incl,tO,0,xdisk,vini,1,*par,**kwargs)
                    return phi/trapezoid(phi,x=w)
            elif Aperture and not Disk :
                def profile(w,w0,vinf,beta,incl,tO,vini,vapr,*par):
                    phi = calc_phi(w,w0,vinf,incl,tO,0,0,vini,vapr,*par,**kwargs)
                    return phi/trapezoid(phi,x=w)
            else:
                def profile(w,w0,vinf,beta,incl,tO,vini,*par):
                    phi = calc_phi(w,w0,vinf,beta,incl,tO,0,0,vini,1,*par,**kwargs)
                    return phi/trapezoid(phi,x=w)
        elif 'hollow' in Geometry.lower() or 'cavity' in Geometry.lower() or 'open' in Geometry.lower():
            if 'fix' in Geometry.lower():
                if Aperture and Disk :
                    def profile(w,w0,vinf,beta,incl,tO,xdisk,vini,vapr,*par):
                        phi = calc_phi(w,w0,vinf,beta,incl,tO,tO-0.174533,xdisk,vini,vapr,*par,**kwargs)
                        return phi/trapezoid(phi,x=w)
                elif Disk and not Aperture :
                    def profile(w,w0,vinf,beta,incl,tO,xdisk,vini,*par):
                        phi = calc_phi(w,w0,vinf,beta,incl,tO,tO-0.174533,xdisk,vini,1,*par,**kwargs)
                        return phi/trapezoid(phi,x=w)
                elif Aperture and not Disk :
                    def profile(w,w0,vinf,beta,incl,tO,vini,vapr,*par):
                        phi = calc_phi(w,w0,vinf,beta,incl,tO,tO-0.174533,0,vini,vapr,*par,**kwargs)
                        return phi/trapezoid(phi,x=w)
                else:
                    def profile(w,w0,vinf,beta,incl,tO,vini,*par):
                        phi = calc_phi(w,w0,vinf,beta,incl,tO,tO-0.174533,0,vini,1,*par,**kwargs)
                        return phi/trapezoid(phi,x=w)
            else:
                if Aperture and Disk :
                    def profile(w,w0,vinf,beta,incl,tO,tC,xdisk,vini,vapr,*par):
                        phi = calc_phi(w,w0,vinf,beta,incl,tO,tC,xdisk,vini,vapr,*par,**kwargs)
                        return phi/trapezoid(phi,x=w)
                elif Disk and not Aperture :
                    def profile(w,w0,vinf,beta,incl,tO,tC,xdisk,vini,*par):
                        phi = calc_phi(w,w0,vinf,beta,incl,tO,tC,xdisk,vini,1,*par,**kwargs)
                        return phi/trapezoid(phi,x=w)
                elif Aperture and not Disk :
                    def profile(w,w0,vinf,beta,incl,tO,tC,vini,vapr,*par):
                        phi = calc_phi(w,w0,vinf,beta,incl,tO,tC,0,vini,vapr,*par,**kwargs)
                        return phi/trapezoid(phi,x=w)
                else:
                    def profile(w,w0,vinf,beta,incl,tO,tC,vini,*par):
                        phi = calc_phi(w,w0,vinf,beta,incl,tO,tC,0,vini,1,*par,**kwargs)
                        return phi/trapezoid(phi,x=w)
        else:
            print('Geometry not recognized. Options are \'Sphere\',\n'+
                '\'FilledCones\', \'HollowCones\', or \'HollowConesFixedCavity\'.')
    return profile
'''
Name:
    profile

Purpose:
    Calculate the spectral line profile for a galactic outflow
    under the Sobolev approximation for the user-specified geometry,
    density profile, and velocity field.

Arguments:
    :wave   (*np.ndarray*) : array of observed wavelengths
    :wave0  (*float*) : rest-frame central wavelength
    :vinf   (*float*) : terminal velocity in c units
    :beta   (*float*) : power-law index of velocity field
    :incl   (*float*) : inclination of the cone
                        (only if filled or hollow conical Geometry)
    :thetaO (*float*) : opening angle of the filled cone
                        (only if filled or hollow conical Geometry)
    :thetaC (*float*) : opening angle of cavity in cone
                        (only if hollow cone Geometry)
    :xdisk  (*float*) : scaled disk radius (only provide if Disk)
    :vapr   (*float*) : scaled aperture limit on velocity
                        (only if correcting for a circular aperture)
    :p1,p2,...: (*floats*) : density field parameters, must match selected
                        'DP' -- see documentation

Returns:
    :profile (*function*) : function which computes the normalized line profile
                        in units of wavelength^(-1) for the given options
'''
