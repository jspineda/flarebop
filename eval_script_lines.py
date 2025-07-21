
from matplotlib import pyplot as plt
import numpy as np
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u


#%%: input stellar parameters

"""AD Leo used as example."""

coords = SkyCoord((154.90117008*u.deg, 19.8700029*u.deg))
# example alternate option copying in simbad coordinates
# coords = SkyCoord("10 19 36.2808181226 +19 52 12.010446571",
#                   unit=(u.hourangle, u.deg))
# coordinates are used to estimate ISM velocity and do not need to be exact

distance = 4.97*u.pc

F = dict(
    o1 = None,
    c2 = None,
    si3 = None,
    si4 = None,
    n5 = None
)


#%% --OR-- input the rotation period
#todo
#%% --OR-- input Ca II and/or Ha measurements
#todo
"""Use melbourne. It would be helpful to have a tool to estimate LHa/Lbol based on EW here."""

#%%