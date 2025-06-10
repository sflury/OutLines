from numpy import pi, inf
# dictionary of all parameter names
ProfName = {'Static':
            {'Emission':['FluxStatic'],\
            'Absorption':['LogColumnStatic'],\
            'PCygni':['LogColumnStatic','FluxStatic']},\
            'Outflow':
            {'Emission':['FluxOutflow'],\
            'Absorption':['LogColumnOutflow'],\
            'PCygni':['LogColumnOutflow','FluxOutflow']}}
GeomName = {'Spherical':[],\
            'FilledCones':['Inclination','OpeningAngle'],\
            'HollowConesFixedCavity':['Inclination','OpeningAngle'],\
            'HollowCones':['Inclination','OpeningAngle','CavityAngle']}
DensName = {'PowerLaw':['PowerLawIndex'],\
            'Exponential':['DecayRate'],\
            'PowerLaw2':['PowerLawIndex1','OuterPowerLawIndex2',\
                              'Inflection'],\
            'Normal':['Peak','Dispersion'],\
            'LogNormal':['LogPeak','LogDispersion'],\
            'Shell':['Location','Width'],\
            'FRED':['RiseRate','DecayRate','Peak'],\
            'Pulses':['PulseWidth','PulseInterval','PulsePhase'],\
            'DampedPulses':['DecayRate','PulseWidth','PulseInterval','PulsePhase']}
# dictionary of all parameter labels
ProfLabs = {'Static':
            {'Emission':[r'$F_s$'],\
            'Absorption':[r'$\log N_s$'],\
            'PCygni':[r'$\log N_s$',r'$F_s$']},\
            'Outflow':
            {'Emission':[r'$F_o$'],\
            'Absorption':[r'$\log N_o$'],\
            'PCygni':[r'$\log N_o$',r'$F_o$']}}
GeomLabs = {'Spherical':[],\
            'FilledCones':[r'$i$',r'$\theta_o$'],\
            'HollowConesFixedCavity':[r'$i$',r'$\theta_o$'],\
            'HollowCones':[r'$i$',r'$\theta_o$',r'$\theta_c$']}
DensLabs = {'PowerLaw':[r'$\alpha$'],\
            'Exponential':[r'$\gamma$'],\
            'PowerLaw2':[r'$\alpha_1$',r'$\alpha_2$',r'$x_1$'],\
            'Normal':[r'$x_1$',r'$\sigma_x$'],\
            'LogNormal':[r'$\log x_1$',r'$\log\sigma_x$'],\
            'Shell':[r'$x_1$',r'$\sigma_x$'],\
            'FRED':[r'$r_1$',r'$r_2$',r'$x_1$'],\
            'Pulses':['$\sigma_x$','$\Delta x$','$x_0$'],\
            'DampedPulses':[r'$\gamma$','$\sigma_x$','$\Delta x$','$x_0$']}
# dictionary of all parameter initial guesses
BetaPars = {'BetaCAK':[1],\
            'AccPlaw':[2],\
            'VelPlaw':[2]}
ProfPars = {'Static':
            {'Emission':[1],\
            'Absorption':[14],\
            'PCygni':[14,1]},\
            'Outflow':
            {'Emission':[1],\
            'Absorption':[14],\
            'PCygni':[14,1]}}
GeomPars = {'Spherical':[],\
            'FilledCones':[pi/6,pi/6],\
            'HollowConesFixedCavity':[pi/6,pi/6],\
            'HollowCones':[pi/3,pi/4,pi/6]}
DensPars = {'PowerLaw':[2],\
            'Exponential':[0.1],\
            'PowerLaw2':[1,2,3],\
            'Normal':[3,1],\
            'LogNormal':[0.5,0.1],\
            'Shell':[3,1],\
            'FRED':[1,1,3],\
            'Pulses':[0.2,1,2],\
            'DampedPulses':[0.1,0.2,1,2]}
# dictionary of all recommended parameter bounds
BetaBounds = {'BetaCAK':[[0.1],[5]],\
            'AccPlaw':[[1.001],[5]],\
            'VelPlaw':[[0.1],[5]]}
ProfBounds = {'Static':
            {'Emission':[[0],[inf]],\
            'Absorption':[[10],[25]],\
            'PCygni':[[10,0],[25,inf]]},\
            'Outflow':
            {'Emission':[[0],[inf]],\
            'Absorption':[[10],[25]],\
            'PCygni':[[10,0],[25,inf]]}}
GeomBounds = {'Spherical':[[],[]],\
            'FilledCones':[[0,0],[pi/2,pi/2]],\
            'HollowConesFixedCavity':[[0,0.174],[pi/2,pi/2]],\
            'HollowCones':[[0,0,0],[pi/2,pi/2,pi/2]]}
DensBounds = {'PowerLaw':[[0],[10]],\
            'Exponential':[[0],[10]],\
            'PowerLaw2':[[0,0,1],[5,10,10]],\
            'Normal':[[1,0.01],[100,10]],\
            'LogNormal':[[0,0.01],[2,1]],\
            'Shell':[[1,0.01],[100,10]],\
            'FRED':[[0.01,0.01,1],[100,100,100]],\
            'Pulses':[[0,0,0],[1,3,10]],\
            'DampedPulses':[[0,0,0,0],[10,1,3,10]]}
AperBounds = [[0],[1]]
DiskBounds = [[1],[inf]]
