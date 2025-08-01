
from matplotlib import pyplot as plt
import numpy as np
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u


#%% settings

verbose = True


#%% input stellar parameters

"""AD Leo used as example."""

coords = SkyCoord((154.90117008*u.deg, 19.8700029*u.deg))
# example alternate option copying in simbad coordinates
# coords = SkyCoord("10 19 36.2808181226 +19 52 12.010446571",
#                   unit=(u.hourangle, u.deg))
# coordinates are used to estimate ISM velocity and do not need to be exact

distance = 4.97*u.pc

"""
If unknown, keep these as the defaults shown. if known, update.
For the ISM values, if multiple components have been identified, use the dominant component.
"""
stellar_rv = None
ism_rv = 'local ism'
log_Nh = 17 * u.Unit('dex(cm-2)')


#%% Stellar activity parameters

"""
There are a range of options for specifiying the stellar activity. 
Each will ultimately be used to infer other quiescent line fluxes of the star, which is the basis of the
relationships from Loyd+ 2018. 
They are listed in order of most to least accurate (modulo some nuances).
Fill out only the topmost line for which you can find information on the star, *no others*.
"""

# direct measurement
# F_si4 = None
F_si4 = 160e-15 * u.Unit('erg s-1 cm-2') # example for AD Leo

# line-line scaling, Pineda & Loyd 2026
F_si3 = None
# F_si3 = 140e-15 * u.Unit('erg s-1 cm-2') # example for AD Leo

# Pineda+ 2021 if Teff <= 3500, Loyd+ 2021 if Teff > 3500
rotation_period, mass, Teff = None, None, None

# less accurate FUV line-line scalings, Pineda & Loyd 2026
Fc2 = None
# F_c2 = 220e-15 * u.Unit('erg s-1 cm-2') # example for AD Leo
Fn5 = None
# F_n5 = 140e-15 * u.Unit('erg s-1 cm-2') # example for AD Leo
Fn5 = None

# Melbourne+ 2020
CaII_S = None
CaII_Rhk = None
LHa_Lbol = None
Ha_EW = None # will need a relationship to convert this to LHa_Lbol

# even less accurate line-line scaling, Pineda & Loyd 2026
F_lya_intrinsic = None

_value_list = [F_si4, F_si3, rotation_period, Fc2, Fn5, CaII_S, CaII_Rhk, Ha_EW, F_lya_intrinsic]
_is_provided = [x is not None for x in _value_list]
if sum(_is_provided) > 1:
    raise ValueError('Values provided for multiple lines in this cell. You must provide only one.')


#%% compute general parameters for worst-case flare


