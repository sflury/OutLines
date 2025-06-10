'''
name:
    OutLines
dependencies:
    python     3.12.9
    numpy      1.26.3
    scipy      1.14.1
    matplotlib 3.10.0
author:
    Sophia Flury 2025.05.07
'''
# OutLines functions
import OutLines.__absorption__
import OutLines.__emission__
import OutLines.__funcs__
import OutLines.__params__
import OutLines.properties
import OutLines.profiles
import OutLines.visualize
# add class calls directly to OutLines for convenience
OutLines.Emission            = OutLines.profiles.Emission
OutLines.Absorption          = OutLines.profiles.Absorption
OutLines.Properties          = OutLines.properties.Properties
OutLines.PlotGeometry        = OutLines.visualize.PlotGeometry
OutLines.PlotConeProjection  = OutLines.visualize.PlotConeProjection
# quote from Shelly's "Ode to the West Wind"
def Quote():
    print('\nWild Spirit, which are moving everywhere;\n'+\
          'Destroyer and preserver\n'+
          'on whose stream, mid the steep sky\'s commotion,\n'+
          'Loose clouds Shook from the tangled boughs of Heaven.\n'+\
          12*' '+' ~ Percy Bysshe Shelley, \"Ode to the West Wind\"')
