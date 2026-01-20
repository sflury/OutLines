<img width="512" alt="OutLines Logo" src="https://github.com/sflury/OutLines/blob/main/docs/logo.png">

OutLines models spectral line profiles from winds, bubbles, and
outflows following the formalism in [Flury 2025](https://ui.adsabs.harvard.edu/abs/2025arXiv251210650F/abstract) 
(see also [Flury, Moran, & Eleazer 2023](https://ui.adsabs.harvard.edu/abs/2023MNRAS.525.4231F)
for an earlier version). A primary goal of OutLines is to remain agnostic to the
underlying physics while also drawing on physical motivations for the geometry,
velocity, and gas density distributions. A cartoon of the model
illustrating the observation of spherical outflows is shown below,
depecting the Doppler shift (colored arrows) of light emitted (yellow and
orange) and absorbed (orange) by gas in the outflow.

<img width="512" alt="image of a model outflow" src="https://github.com/sflury/OutLines/assets/42982705/9af5bf13-d2ce-441b-b429-294833ae5edc">

Physically justifiable assumptions include the density profile, including a
variety of continuous or shell-like gas density distributions,
and the velocity field, including the so-called beta law from approximations
to CAK theory (e.g., Castor et al. 1975, Barlow et al. 1977) and other power law
solutions, under the Sobolev approximation that small-scale ("local")
gas velcities contribute negligibly to the net velocity field.

Emission and absorption profiles are computed in velocity space for the
specified wavelength(s). Following the equation of radiative transfer, emission
line profiles should be added in flux space while absorption line profiles
should be added in optical depth space. One or multiple lines can be computed
for each set of outflow properties, which is highly recommended in the case of
multiple features from the same species such as the \[O III\] doublet or Si II
UV absorption features. However, simultaneous fitting of lines from different
phases of the interstellar medium or even different ionization zones is
cautioned as different phases may not share the same wind, bubble, or
outflow properties.

OutLines currently supports line profile models for absorption lines and
nebular, resonant, and fluorescent emission lines.

## Installation

While it is possible to download this repository, OutLines is readily accessible
via `pip install`, which will automatically ensure the appropriate dependencies
are also installed. An example call in the terminal command line is shown below.

``` bash
$ pip3 install SpecOutLines
```

## Line Profile Classes

| CLASS       | LINE FEATURE         |
|-------------|----------------------|
| Absorption  | resonant absorption  |
| Nebular     | nebular emission     |
| Resonant    | resonant emission    |
| Fluorescent | fluorescent emission |


## Example Usage -- \[O III\] 4959,5007 Profiles
``` python
import OutLines as OL
from numpy import linspace
import matplotlib.pyplot as plt
model = OL.Nebular([4958.911,5006.843],Geometry='HollowCones',AddStatic=True,Disk=True)
model.update_params(['TerminalVelocity'],[500,30,15,45])
model.update_params(['OpeningAngle','CavityAngle','Inclination'],[500,30,15,45])
model.update_params(['FluxOutflow1','FluxStatic1'],[1/2.98,1/2.98])
wave = linspace(4950,5015,651)wave = linspace(4950,5015,651)
plt.plot(wave,model.get_profile(wave),lw=2,color='C3')
plt.plot(wave,model.get_outflow(wave),dashes=[3,3],lw=2,color='C3')
plt.plot(wave,model.get_static(wave),':',lw=2,color='C3')
plt.show()
model.print_settings()
model.print_params()
props = OL.Properties(model)
props.print_props()
```
<img width="480" alt="image of predicted \[O III\] doublet profile" src="https://github.com/sflury/OutLines/blob/main/examps/o_iii.png">

``` text
----------------------------------
|          MODEL SETTINGS          |
----------------------------------
|     Line   1     : 4958.911      |
|     Line   2     : 5006.843      |
|          Profile : Nebular       |
|    VelocityField : BetaCAK       |
|   DensityProfile : PowerLaw      |
|         Geometry : HollowCones   |
|  StaticComponent : Yes           |
|         Aperture : No            |
|             Disk : Yes           |
----------------------------------
 ----------------------------------
|         MODEL PARAMETERS         |
 ----------------------------------
|     DopplerWidth :    8.994 km/s |
| TerminalVelocity :  500.000 km/s |
|    VelocityIndex :    1.000      |
|      FluxStatic1 :    0.336      |
|      FluxStatic2 :    1.000      |
|     FluxOutflow1 :    0.336      |
|     FluxOutflow2 :    1.000      |
|      Inclination :   45.000°     |
|     OpeningAngle :   30.000°     |
|      CavityAngle :   15.000°     |
|       DiskRadius :    2.000      |
|    PowerLawIndex :    2.000      |
 ----------------------------------
 -------------------------------------------
|             MODEL PROPERTIES              |
 -------------------------------------------
|        x.out :    1.500                   |
|        v.out :  166.667  km s^-1          |
|         Mdot :    1.057  Msun yr^-1       |
|         pdot :    0.111  10^34 dyne       |
|         Edot :    0.009  10^42 erg s^-1   |
|        v.esc :    1.667                   |
|     pdot.esc :    0.093                   |
|     Edot.esc :    0.154                   |
 -------------------------------------------
| Mdot, pdot, Edot / R0^2 n0 [kpc^2 cm^-3]  |
|  pdot.est, Edot.esc / v0 [100 km s^-1]    |
 -------------------------------------------
```


## Example Usage -- Si II 1260 Profile
``` python
import OutLines as OL
from numpy import linspace
import matplotlib.pyplot as plt
kwargs = dict(Geometry='FilledCones',DensityProfile='LogNormal',AddStatic=True)
model = OL.Absorption(1260.4221,1.18,**kwargs)
wave = linspace(1259.25,1260.75,1001)
plt.plot(wave,model.get_profile(wave),lw=2,color='C3')
plt.plot(wave,model.get_outflow(wave),dashes=[3,3],lw=2,color='C3')
plt.plot(wave,model.get_static(wave),':',lw=2,color='C3')
plt.show()
```
<img width="480" alt="image of predicted Si II 1260 absorption profile" src="https://github.com/sflury/OutLines/blob/main/examps/si_ii.png">

## Referencing `OutLines`

While this code is provided publicly, it did require substantial effort to
develop and document. Any use thereof must be cited in any publications in which
this code is used. The BibTeX reference for the [Flury (2025)](https://arxiv.org/abs/2512.10650) paper which 
presents the models and code is below; however, a GitHub CCF is also provided 
for convenience.

``` bibtex
@ARTICLE{Flury2025,
       author = {{Flury}, Sophia R.},
        title = "{OutLines: Modeling Astrophysical Winds, Bubbles, and Outflows}",
      eprint={2512.10650},
      archivePrefix={arXiv},
      primaryClass={astro-ph.GA},
      url={https://arxiv.org/abs/2512.10650}, 
         year = {2025},
        month = {dec} }
```

I developed and implemented the spherical geometry, power law
density, CAK approximation model for analysis of broad [O III] lines
observed in Mrk 462. That model was presented in
[Flury, Moran, & Eleazer (2023) MNRAS 525, 4231](https://ui.adsabs.harvard.edu/abs/2023MNRAS.525.4231F)
The BibTeX reference is below.

``` bibtex
@ARTICLE{Flury2023,
       author = {{Flury}, Sophia R. and {Moran}, Edward C. and {Eleazer}, Miriam},
        title = "{Galactic outflow emission line profiles: evidence for dusty, radiatively driven ionized winds in Mrk 462}",
      journal = {\mnras},
         year = 2023,
        month = nov,
       volume = {525},
       number = {3},
        pages = {4231-4242},
          doi = {10.1093/mnras/stad2421} }
```

## Licensing
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.
