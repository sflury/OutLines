# OutLines
Computes spectral line profiles arising from galactic outflows following 
formalism and examples in Flury, Moran, & Eleazer (2023) MNRAS while 
remaining agnostic to the underlying physics.
Script currently supports nebular emission lines which do not undergo 
self-absorption (and are therefore always optically thin, i.e. "pure" 
emission) and absorption lines without infilling effects (although line 
infilling can be done manually to produce
P-Cygni line profiles).

Physically justifiable assumptions include the density profile
$$n \propto r^{-\alpha}$$
and velocity profile (from approximations to CAK theory)
$$v \propto (1-r^{-1})^{\beta}$$
under the Sobolev approximation that small-scale ("local")
gas velcities contribute negligibly to the net velocity field.

While this code is provided publicly, I request that any use thereof be 
cited in any publications in which this code is used. I developed and 
implemented this script for Flury, Moran, & Eleazer (2023) MNRAS with
example applications to the \[O III\] line in Mrk 462.

## Example Usage -- \[O III\] Profile for Mrk 462
```
from numpy import arange
from OutLines import *
# speed of light in km/s
c    = 2.99792458e5
# rest-frame wavelengths
wr = arange(4900,5050,0.25)
# predict line profiles for both [O III] doublet transitions using Mrk 462 results
oiii_outflow = 0.332*f5007 * phi_out(wr,5007-47.932,745/c,1,1.3) + f5007 * phi_out(wr,5007,745/c,1.122,1.369)
```
![image of predicted \[O III\] doublet profile](oiii_examp.png "[OIII]4959,5007 profile")


## BibTex
*pending*

## Licensing
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.
