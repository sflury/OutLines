from OutLines.__funcs__ import *
# velocity integral limits
def calc_limits(ww0,vinf,vini,vapr,beta,VF,Inflow=False):
    # error check for length of w
    if not hasattr(ww0,'__len__'):
        w = array([ww0])
    # observed velocity from wavelengths, in vinf/c units
    y = absolute(((ww0)**2-1)/((ww0)**2+1)) / vinf
    ind = where(y<1)[0]
    # u = v/vinf (not y since y here is used for vobs/v)
    # with no occultation, the minimum velocity is just the observed velocity
    umin = y[ind]
    # primary cone type -- start with all as anterior
    cone = array(len(umin)*['ante']).astype(str)
    # when occultation by source occurs, need to solve
    # for minimum contributing velocity to observed y
    if Inflow :
        umin[ww0[ind]<1] = array(list(map(partial(solve_u0,beta,vini,VF),y[ind][ww0[ind]<1])))
        cone[ww0[ind]<1] = 'post'
    else:
        umin[ww0[ind]>1] = array(list(map(partial(solve_u0,beta,vini,VF),y[ind][ww0[ind]>1])))
        cone[ww0[ind]>1] = 'post'
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
def resfluor(w,u,vinf,beta,vini,geo,cone,par,VF,DP):
    # polar geometry
    omega  = vinf*w                                    # relative velocity
    Lrnzt1 = sqrt(1-omega**2)                          # Lorenzt factor
    DopRel = square(1-omega)*sqrt((1+omega)/(1-omega)) # relative Doppler shift
    dTheta = Lrnzt1/(omega*DopRel)                     # total polar term
    # azimuthal geometry
    # sphere
    if geo[1] >= 1-2**-20 and geo[8] <= 2**-20 :
        ell = 1
    # hemisphere
    elif geo[1] >= 1-2**-20 and geo[8] >= 1-2**-20 :
        ell = array(list(map(partial(hemi_inscr,u,vinf,*geo[6:11],cone),w)))
    # cones
    elif geo[1] < 1-2**-20 and geo[8] < 2**-20 :
        ell  = array(list(map(partial(cone_inscr,u,vinf,*geo[:11],cone),w)))
        # cavity in cone
        if geo[12] > 0. :
            ell -= array(list(map(partial(cone_inscr,u,vinf,*geo[11:],cone),w)))
    # cones with disk
    else:
        ell  = array(list(map(partial(cndk_inscr,u,vinf,*geo[:11],cone),w)))
        # cavity in cone
        if geo[12] > 0. :
            ell -= array(list(map(partial(cndk_inscr,u,vinf,*geo[11:],cone),w)))
        ell = nanmax([ell,ell_zeros],axis=0)
    return ell*n[DP](w,beta,vini,VF,*par)*dxdw(w,beta,vini,VF)* dTheta
# integral over column densities for range of allowed velocities
def phi_int(vinf,beta,incl,tO,tC,vdisk,vini,par,VF,DP,u,umin,umax,cone):
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
        return fixed_quad(resfluor,umin,umax,args=(u,vinf,beta,vini,geo,cone,par,VF,DP))
# calculate unnormalized profile for a sphere or bicone
def calc_phi(ww0,vinf,beta,incl,tO,tC,xdisk,vini,vapr,*par,VF='BetaCAK',DP='PowerLaw',Pulse='Normal'):
    # for ensembles, call recursively for each pulse
    if 'Pulse' in DP :
        phi_sum = zeros(len(ww0))
        check   = 1
        scale   = 1
        if 'Damp' in DP : xi = par[3]
        else:             xi = par[2]
        # while loop for emission since largest shells have smallest flux
        # only keep computing if the flux is > 0.001, v < 99.9% terminal
        # --> or, for damped pulses, if > 7 e-foldings, >99.9% of total reached
        while ( check > 2**-10 and xi < x[VF](0.999,beta,vini) ) or xi < 5 :
            if 'Damp' in DP :
                par1 = [xi,par[1]]
                scale = exp(-par[0]*xi)
                xi += par[2]
            else:
                par1 = [xi,par[0]]
                xi += par[1]
            phi_pls = scale*calc_phi(ww0,vinf,beta,incl,tO,tC,xdisk,vini,vapr,\
                                                    *par1,VF=VF,DP=Pulse)
            phi_sum += phi_pls
            check = max(phi_pls)
        return phi_sum
    # obtain velocity limits for integral
    umin,umax,cone = calc_limits(ww0,vinf,vini,vapr,beta,VF)
    # improve precision for bubble/shell integrals by limiting the radial range
    if 'LogNormal' in DP :
        umin = max([umin,len(umin)*[v[VF](10**(par[0]-3*par[1]),beta,vini)]],axis=0)
        umax = min([umax,len(umax)*[v[VF](10**(par[0]+3*par[1]),beta,vini)]],axis=0)
    # 99.9 percentile for normal distribution
    elif 'Normal' in DP or 'Shell' in DP :
        umin = nanmax([umin,len(umin)*[v[VF](par[0]-3*par[1],beta,vini)]],axis=0)
        umax = nanmin([umax,len(umax)*[v[VF](par[0]+3*par[1],beta,vini)]],axis=0)
    # line of sight velocity with relativistic corrections
    u = absolute(((ww0)**2-1)/((ww0)**2+1)) / vinf
    # check disk radius and convert to velocity
    if xdisk <= 1 :
        vdisk = 2**-20
    elif xdisk == inf :
        vdisk = 1
    else:
        vdisk = max([v[VF](xdisk,beta,vini),2**-20])
    # profile integrated at each velocity
    phi = zeros(len(ww0))
    phi[ u < 1 ] = array( list( map( \
            partial(phi_int,vinf,beta,incl,tO,tC,vdisk,vini,par,VF,DP),\
            u[ u < 1 ], umin, umax, cone ) ) )
    return phi
