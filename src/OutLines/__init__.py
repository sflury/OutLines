'''
name:
    OutLines
dependencies:
    python     3.12.9
    numpy      2.3.5
    scipy      1.14.1
    matplotlib 3.10.0
author:
    Sophia Flury 2026.01.20
'''
# OutLines functions
import OutLines.__absorption__
import OutLines.__nebular__
import OutLines.__resonfluor__
import OutLines.__funcs__
import OutLines.__params__
import OutLines.__properties__
import OutLines.__profiles__
import OutLines.__visualize__
# add class calls directly to OutLines for convenience
OutLines.Nebular             = OutLines.__profiles__.Nebular
OutLines.Absorption          = OutLines.__profiles__.Absorption
OutLines.Resonant            = OutLines.__profiles__.Resonant
OutLines.Fluorescent         = OutLines.__profiles__.Fluorescent
OutLines.Properties          = OutLines.__properties__.Properties
OutLines.PlotGeometry        = OutLines.__visualize__.PlotGeometry
OutLines.PlotConeProjection  = OutLines.__visualize__.PlotConeProjection
OutLines.PlotIsoContours     = OutLines.__visualize__.PlotIsoContours
# quote from Shelly's "Ode to the West Wind"
def Quote():
    print('\nWild Spirit, which art moving everywhere;\n'+\
          'Destroyer and preserver\n'+
          'on whose stream, mid the steep sky\'s commotion,\n'+
          'Loose clouds Shook from the tangled boughs of Heaven.\n'+\
          12*' '+' ~ Percy Bysshe Shelley, \"Ode to the West Wind\"')
