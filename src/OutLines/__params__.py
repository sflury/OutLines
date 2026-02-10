from numpy import pi, inf
# dictionary of all parameter names
ProfName = {'Static':
            {'Nebular':['FluxStatic'],\
            'Absorption':['LogColumnStatic'],\
            'Resonant':['LogColumnStatic'],\
            'Fluorescent':['LogColumnStatic'],\
            'PCygni':['LogColumnStatic']},\
            'Outflow':
            {'Nebular':['FluxOutflow'],\
            'Absorption':['LogColumnOutflow'],\
            'Resonant':['LogColumnOutflow'],\
            'Fluorescent':['LogColumnOutflow'],\
            'PCygni':['LogColumnOutflow']}}
GeomName = {'Spherical':[],\
            'Hemisphere':['Inclination'],\
            'FilledCones':['Inclination','OpeningAngle'],\
            'HollowConesFixedCavity':['Inclination','OpeningAngle'],\
            'HollowCones':['Inclination','OpeningAngle','CavityAngle']}
DensName = {'PowerLaw':['PowerLawIndex'],\
            'Exponential':['DecayRate'],\
            'PowerLaw2':['PowerLawIndex1','PowerLawIndex2','Inflection'],\
            'Normal':['Radius','Width'],\
            'LogNormal':['LogRadius','LogWidth'],\
            'DLogistic':['Radius','GrowthRate'],\
            'Shell':['Radius','Width'],\
            'FRED':['RiseRate','DecayRate','Radius'],\
            'Pulses':['PulseWidth','PulseInterval','InitRadius'],\
            'DampedPulses':['DecayRate',\
                      'PulseWidth','PulseInterval','InitRadius'],\
            'PacketPulses':['PacketWidth','PacketRadius',\
                      'PulseWidth','PulseInterval','InitRadius'],\
            'PulsesLog':['PulseWidth','PulseInterval','InitRadius'],\
            'DampedPulsesLog':['DecayRate',\
                      'PulseWidth','PulseInterval','InitRadius'],\
            'PacketPulsesLog':['PacketWidth','PacketRadius',\
                      'PulseWidth','PulseInterval','InitRadius'],\
            }
# dictionary of all parameter labels
ProfLabs = {'Static':
            {'Nebular':[r'$F_s$'],\
            'Absorption':[r'$\log N_s$'],\
            'Resonant':[r'$\log N_s$'],\
            'Fluorescent':[r'$\log N_s$'],\
            'PCygni':[r'$\log N_s$']},\
            'Outflow':
            {'Nebular':[r'$F_o$'],\
            'Absorption':[r'$\log N_o$'],\
            'Resonant':[r'$\log N_o$'],\
            'Fluorescent':[r'$\log N_o$'],\
            'PCygni':[r'$\log N_o$']}}
GeomLabs = {'Spherical':[],\
            'Hemisphere':[r'$i$'],\
            'FilledCones':[r'$i$',r'$\theta_o$'],\
            'HollowConesFixedCavity':[r'$i$',r'$\theta_o$'],\
            'HollowCones':[r'$i$',r'$\theta_o$',r'$\theta_c$']}
DensLabs = {'PowerLaw':[r'$\alpha$'],\
            'Exponential':[r'$\gamma$'],\
            'PowerLaw2':[r'$\alpha_1$',r'$\alpha_2$',r'$x_1$'],\
            'Normal':[r'$x_1$',r'$\sigma_x$'],\
            'LogNormal':[r'$\log x_1$',r'$\log\sigma_x$'],\
            'DLogistic':['$x_1$','$\alpha$'],\
            'Shell':[r'$x_1$',r'$\sigma_x$'],\
            'FRED':[r'$r_1$',r'$r_2$',r'$x_1$'],\
            'Pulses':[r'$\sigma_x$',r'$\Delta x$',r'$x_0$'],\
            'DampedPulses':[r'$\gamma$',\
                      r'$\sigma_x$',r'$\Delta x$',r'$x_0$'],\
            'PacketPulses':[r'$\sigma_p$','$x_p$',\
                      r'$\sigma_x$',r'$\Delta x$',r'$x_0$'],\
            'PulsesLog':[r'$\sigma_x$',r'$\Delta \log x$',r'$\log x_0$'],\
            'DampedPulsesLog':[r'$\gamma$',\
                      r'$\sigma_x$',r'$\Delta \log x$',r'$\log x_0$'],\
            'PacketPulsesLog':[r'$\sigma_p$','$x_p$',\
                      r'$\sigma_x$',r'$\Delta \log x$',r'$\log x_0$'],\
            }
# dictionary of all parameter initial guesses
BetaPars = {'BetaCAK':[1],\
            'AccPlaw':[2],\
            'VelPlaw':[2]}
ProfPars = {'Static':
            {'Nebular':[1],\
            'Absorption':[14],\
            'Resonant':[14],\
            'Fluorescent':[14],\
            'PCygni':[14]},\
            'Outflow':
            {'Nebular':[1],\
            'Absorption':[14],\
            'Resonant':[14],\
            'Fluorescent':[14],\
            'PCygni':[14]}}
GeomPars = {'Spherical':[],\
            'Hemisphere':[pi/4],\
            'FilledCones':[pi/6,pi/6],\
            'HollowConesFixedCavity':[pi/6,pi/6],\
            'HollowCones':[pi/3,pi/4,pi/6]}
DensPars = {'PowerLaw':[2],\
            'Exponential':[0.1],\
            'PowerLaw2':[1,2,3],\
            'Normal':[3,1],\
            'LogNormal':[0.5,0.1],\
            'DLogistic':[3,1],\
            'Shell':[3,1],\
            'FRED':[1,1,3],\
            'Pulses':[0.2,1,2],\
            'DampedPulses':[0.1,0.2,1,2],\
            'PacketPulses':[1,4,0.2,1,2],\
            'PulsesLog':[0.2,0.1,0.3],\
            'DampedPulsesLog':[0.1,0.2,0.1,0.3],\
            'PacketPulsesLog':[0.1,0.5,0.2,0.1,0.3],\
            }
# dictionary of all recommended parameter bounds
BetaBounds = {'BetaCAK':[[0.1],[5]],\
            'AccPlaw':[[1.001],[5]],\
            'VelPlaw':[[0.1],[5]]}
ProfBounds = {'Static':
            {'Nebular':[[0],[inf]],\
            'Absorption':[[10],[25]],\
            'Resonant':[[10],[25]],\
            'Fluorescent':[[10],[25]],\
            'PCygni':[[10],[25]]},\
            'Outflow':
            {'Nebular':[[0],[inf]],\
            'Absorption':[[10],[25]],\
            'Resonant':[[10],[25]],\
            'Fluorescent':[[10],[25]],\
            'PCygni':[[10],[25]]}}
GeomBounds = {'Spherical':[[],[]],\
            'Hemisphere':[[0],[pi/2]],\
            'FilledCones':[[0,0],[pi/2,pi/2]],\
            'HollowConesFixedCavity':[[0,0.174],[pi/2,pi/2]],\
            'HollowCones':[[0,0,0],[pi/2,pi/2,pi/2]]}
DensBounds = {'PowerLaw':[[0],[10]],\
            'Exponential':[[0],[10]],\
            'PowerLaw2':[[0,0,1],[5,10,10]],\
            'Normal':[[1,0.01],[100,10]],\
            'LogNormal':[[0,0.01],[2,1]],\
            'DLogistic':[[1,0],[100,10]],\
            'Shell':[[1,0.01],[100,10]],\
            'FRED':[[0.01,0.01,1],[100,100,100]],\
            'Pulses':[[0,0,0],[1,3,10]],\
            'DampedPulses':[[0,0,0,0],[10,1,3,10]],\
            'PacketPulses':[[0,0,0,0,0],[3,10,1,3,10]],\
            'PulsesLog':[[0,0,0],[1,1,1]],\
            'DampedPulsesLog':[[0,0,0,0],[10,1,0.5,1]],\
            'PacketPulsesLog':[[0,0,0,0],[0.7,1,1,0.5,1]]
            }
AperBounds = [[0],[1]]
DiskBounds = [[1],[inf]]
