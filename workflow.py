
from matplotlib import pyplot as plt
import numpy as np

from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table

import lines
import ism
import worst_case
import ffd

from starprop import StarProp
from probbop import FlProb


## Need alternate example for something manual and/or where simbad has missing info (fix missing info situation within starprop.py)

"""AD Leo used as example."""

leo = StarProp("AD Leo") # Based on simbad query


"""UV Line Fluxes are Known"""

flux_unit = u.Unit('erg s-1 cm-2')
line_fluxes = dict(
    CII = 220e-15 * flux_unit,
    SiIII = 140e-15 * flux_unit,
    SiIV = 160e-15 * flux_unit,
)

leo.setLines(line_fluxes)


"""UV Line Fluxes are Not Known""" ## Method requires a good radius

leo.setPeriod(2.2270*u.d)
leo.populateLines()

#leo.setAge(250*u.Myr)
#leo.populateLines()


print(leo.all_lines) # see assumed quiescent line fluxes




"""Overlight Probability is Based on FFDs"""

uvFFD = ffd.dfuv130
uvFFD._set_norm(ffd.totflarerate) ## this stuff can be moved under the hood, it's here for now

leo.setFFDuv(uvFFD)


## will want to move more of this under the hood (e.g., choosing contrast of indivudal lines matching chosen observation mode), here for now until observatory class is more mature
fl_correlations = Table.read('tables/flare_line_fits.ecsv')
_select = ( (fl_correlations['xparam'] == "equiv dur") &
                (fl_correlations['yparam'] == 'contrast') &
                (fl_correlations['xband'] == 'fuv130') &
                (fl_correlations['percentile'] ==  "50"))

leo.setContrast(fl_correlations[_select])


## with proper obsclass, determination of amplitude H will be computed from obs mode in relation to quiescent flux
obsclass = None

flareBOP = FlProb(leo,obsclass)

flareBOP.runDiag()

## Example Spectra ---

# in order to resolve the lines, do not use a grid step below 0.05 AA
wavegrid = np.arange(1100, 1800, 0.01) * u.AA

