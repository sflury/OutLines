from scipy.integrate import fixed_quad
from functools import partial
from OutLines.__funcs__ import *
#from funcs import *
# velocity integral limits
def calc_limits(w,w0,vinf,vapr,beta,VelocityField,Inflow=False):
    # error check for length of w
    if not hasattr(w,'__len__'):
        w = array([w])
    # observed velocity from wavelengths
    y = absolute(w-w0)/w0 /vinf
    # primary cone type
    if Inflow :
        ind = where((y<1)&(w>w0))[0]
        cone = array(len(ind)*['ante'])
    else:
        ind = where((y<1)&(w<w0))[0]
        cone = array(len(ind)*['ante'])
    # u = v/vinf (not y since y here is used for vobs/v)
    # with no occultation, the minimum velocity is just the observed velocity
    umin = y[ind]
    # by definition, velocity cannot exceed terminal velocity
    # since u = v/vinf, max value is simply 1
    # umax = ones(len(ind))
    # but must account for backlighting by source as maximum possible velocity
    umax = array(list(map(partial(solve_u,beta,VelocityField),y[ind])))
    wapr = sqrt( y[ind]**2 + (vapr/vinf)**2 )
    umax = min([umax,wapr],axis=0)
    return umin,umax,cone
# absorption per unit velocity
# cos t = v_obs / v_inf = u / v_inf = u if u in v_inf units
# u = cos t
# sin t = sqrt(1-cos^2 t) = sqrt(1-u^2)
def absorp(w,u,beta,geo,cone,par,VelocityField,DensityProfile):
    if geo[1] == 1.:
        bl = 1.
    else:
        bl  = array(list(map(partial(cone_bands,u,*geo[:10],cone),w)))
    if geo[9] > 0 :
        bl = bl - array(list(map(partial(cone_bands,u,*geo[10:],cone),w)))
        bl = max([bl,zeros(256)],axis=0)
    return bl*n[DensityProfile](w,beta,VelocityField,*par)*dxdw(w,beta,VelocityField)
# integral over emissivities for range of allowed velocities
def phi_int(beta,incl,tO,tC,vdisk,par,VelocityField,DensityProfile,u,umin,umax,cone):
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
        return fixed_quad(absorp,umin,umax,args=(u,beta,geo,cone,par,VelocityField,DensityProfile),n=256)[0]
# calculate unnormalized profile for a sphere or bicone
def calc_phi(w,w0,vinf,beta,incl,tO,tC,xdisk,vapr,*par,VelocityField='BetaCAK',DensityProfile='PowerLaw'):
    # obtain velocity limits for integral
    umin,umax,cone = calc_limits(w,w0,vinf,vapr,beta,VelocityField)
    u = absolute(w-w0)/w0 / vinf
    # check disk radius and convert to velocity
    if xdisk <= 1 :
        vdisk = 2**-10
    else:
        vdisk = max([v[VelocityField](xdisk,beta),2**-10])
    # profile integrated at each velocity
    phi = zeros(len(w))
    phi[ (u < 1) & (w < w0) ] = array(list(map( \
        partial(phi_int,beta,incl,tO,tC,vdisk,par,VelocityField,DensityProfile), \
        u[ (u < 1) & (w < w0) ], umin, umax, cone )))
    return max([phi,zeros(len(w))],axis=0)

'''
Name:
    build_profile_model

Purpose:
    Construct a unique function based on user-specified model options. This
    function will then compute and return the normalized line profiles.

Keyword Arguments:
    :VelocityField  (*str*): string indicating the desired velocity field. Options are
    :DensityProfile (*str*): string indicating the desired density profile. Options are
                        \'PowerLaw\',
                        \'PowerLaw2\',
                        \'Exponential\',
                        \'Normal\',
                        \'LogNormal\',
                        \'FRED\',
                        \'Shell\'.
                        Default is \'PowerLaw\'.
    :Geometry (*str*): string indicating the desired outflow Geometry. Options are
                        \'Sphere\'
                        \'FilledCones\'
                        \'HollowCones\' or \'CavityCones\'
    :Aperture (*bool*): boolean indicating whether an aperture correction
                        is to be applied to the velocity limits of the integrals
    :Disk     (*bool*): boolean indicating whether an obstructing disk is to be
                        included in calculation of the cone projections

Returns:
    :profile (*function*): a function which calculates the line profile for a
                        galactic outflow for the user-provided model settings

'''
def build_profile_model(VelocityField='BetaCAK',DensityProfile='PowerLaw',Geometry='Sphere',Aperture=False,Disk=False):
    kwargs = dict(VelocityField=VelocityField,DensityProfile=DensityProfile)
    if 'spher' in Geometry.lower():
        if Aperture :
            def profile(w,w0,vinf,beta,vapr,*par):
                phi = calc_phi(w,w0,vinf,beta,0,pi/2,0,0,vapr,*par,**kwargs)
                return phi/trapz(phi,x=w)
        else:
            def profile(w,w0,vinf,beta,*par):
                phi = calc_phi(w,w0,vinf,beta,0,pi/2,0,0,1,*par,**kwargs)
                return phi/trapz(phi,x=w)
    elif 'filled' in Geometry.lower():
        if Aperture and Disk :
            def profile(w,w0,vinf,beta,incl,tO,xdisk,vapr,*par):
                phi = calc_phi(w,w0,vinf,beta,incl,tO,0,xdisk,vapr,*par,**kwargs)
                return phi/trapz(phi,x=w)
        elif Disk and not Aperture :
            def profile(w,w0,vinf,beta,incl,tO,xdisk,*par):
                phi = calc_phi(w,w0,vinf,beta,incl,tO,0,xdisk,1,*par,**kwargs)
                return phi/trapz(phi,x=w)
        elif Aperture and not Disk :
            def profile(w,w0,vinf,beta,incl,tO,vapr,*par):
                phi = calc_phi(w,w0,vinf,incl,tO,0,0,vapr,*par,**kwargs)
                return phi/trapz(phi,x=w)
        else:
            def profile(w,w0,vinf,beta,incl,tO,*par):
                phi = calc_phi(w,w0,vinf,beta,incl,tO,0,0,1,*par,**kwargs)
                return phi/trapz(phi,x=w)
    elif 'hollow' in Geometry.lower() or 'cavity' in Geometry.lower():
        if 'fix' in Geometry.lower():
            if Aperture and Disk :
                def profile(w,w0,vinf,beta,incl,tO,xdisk,vapr,*par):
                    phi = calc_phi(w,w0,vinf,beta,incl,tO,tO-0.174,xdisk,vapr,*par,**kwargs)
                    return phi/trapz(phi,x=w)
            elif Disk and not Aperture :
                def profile(w,w0,vinf,beta,incl,tO,tC,xdisk,*par):
                    phi = calc_phi(w,w0,vinf,beta,incl,tO,tO-0.174,xdisk,1,*par,**kwargs)
                    return phi/trapz(phi,x=w)
            elif Aperture and not Disk :
                def profile(w,w0,vinf,beta,incl,tO,tC,vapr,*par):
                    phi = calc_phi(w,w0,vinf,beta,incl,tO,tO-0.174,0,vapr,*par,**kwargs)
                    return phi/trapz(phi,x=w)
            else:
                def profile(w,w0,vinf,beta,incl,tO,tC,*par):
                    phi = calc_phi(w,w0,vinf,beta,incl,tO,tO-0.174,0,1,*par,**kwargs)
                    return phi/trapz(phi,x=w)
        else:
            if Aperture and Disk :
                def profile(w,w0,vinf,beta,incl,tO,tC,xdisk,vapr,*par):
                    phi = calc_phi(w,w0,vinf,beta,incl,tO,tC,xdisk,vapr,*par,**kwargs)
                    return phi/trapz(phi,x=w)
            elif Disk and not Aperture :
                def profile(w,w0,vinf,beta,incl,tO,tC,xdisk,*par):
                    phi = calc_phi(w,w0,vinf,beta,incl,tO,tC,xdisk,1,*par,**kwargs)
                    return phi/trapz(phi,x=w)
            elif Aperture and not Disk :
                def profile(w,w0,vinf,beta,incl,tO,tC,vapr,*par):
                    phi = calc_phi(w,w0,vinf,beta,incl,tO,tC,0,vapr,*par,**kwargs)
                    return phi/trapz(phi,x=w)
            else:
                def profile(w,w0,vinf,beta,incl,tO,tC,*par):
                    phi = calc_phi(w,w0,vinf,beta,incl,tO,tC,0,1,*par,**kwargs)
                    return phi/trapz(phi,x=w)
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
    :thetaI (*float*) : opening angle of cavity in cone
                        (only if hollow cone Geometry)
    :xdisk  (*float*) : scaled disk radius (only provide if Disk)
    :vapr   (*float*) : scaled aperture limit on velocity
                        (only if correcting for a circular aperture)
    :p1,p2,...: (*floats*) : density field parameters, must match selected
                        'DensityProfile' -- see documentation

Returns:
    :profile (*function*) : function which computes the normalized line profile
                        in units of wavelength^(-1) for the given options
'''
