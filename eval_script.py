from pathlib import Path

from matplotlib import pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u

import lines
import ism
import worst_case

#%% settings

# in order to resolve the lines, do not use a grid step below 0.05 AA
wavegrid = np.arange(1100, 1800, 0.01) * u.AA


#%% input stellar parameters

"""AD Leo used as example."""

coords = SkyCoord(154.90117008*u.deg, 19.8700029*u.deg)
# example alternate option copying in simbad coordinates
# coords = SkyCoord("10 19 36.2808181226 +19 52 12.010446571",
#                   unit=(u.hourangle, u.deg))
# coordinates are used to estimate ISM velocity and do not need to be exact

distance = 4.97*u.pc

"""
If unknown, keep these as the defaults shown. if known, update.
For the ISM values, if multiple components have been identified, use the dominant component.
"""
rv_star = None
rv_ism = 'local ism'
Nh = 17 * u.Unit('dex(cm-2)') # note dex units mean this is the log10 of the value in cm-2


#%% Stellar activity parameters

"""
For now the only option is to specify some number of known UV line fluxes,
but in the future we should add the ability to specify a rotation period,
R'hk, S, or Ha activity metric if the star has no UV measurements.
I'm tempted to make a correlation with GALEX mags too, but that
could be a rabit hole.
"""

flux_unit = u.Unit('erg s-1 cm-2')
line_fluxes = dict(
    CII = 220e-15 * flux_unit,
    SiIII = 140e-15 * flux_unit,
    SiIV = 160e-15 * flux_unit,
)


#%% fill any reamining line fluxes based on line-line correlations

dist_to_1AU = (distance / u.au) ** 2
dist_to_1AU = dist_to_1AU.to_value('')
line_fluxes_1AU = {line : flux * dist_to_1AU for line, flux in line_fluxes.items()}
all_line_fluxes_1AU = lines.fill_quiescent_lines(**line_fluxes_1AU)
all_line_fluxes = {line : flux / dist_to_1AU for line, flux in all_line_fluxes_1AU.items()}

#%% optional: make a quiescent spectrum
"""This could be useful for SNR estimates."""

qspec = lines.spectrum(wavegrid, **all_line_fluxes)

#%% make a spectrum for a worst case flare

fluxes_worst_case = {}
for line, quiescent_flux in all_line_fluxes.items():
    flare_contrast = worst_case.line_contrasts[line]
    flare_flux = quiescent_flux * flare_contrast
    fluxes_worst_case[line] = flare_flux
fspec = lines.spectrum(wavegrid, **fluxes_worst_case)


#%% add ISM absorption

if rv_star is None:
    print('No ISM absorption added because stellar RV not provided.')
else:
    if rv_ism == 'local ism':
        rv_ism = ism.local_ism_rv(coords)
    rv_diff = rv_ism - rv_star
    transmission = ism.transmission(wavegrid, rv_diff, Nh.to('cm-2'))
    fspec *= transmission
    try:
        qspec *= transmission
    except NameError:
        pass


#%% plot spectrum, if desired

def setup_plot():
    plt.figure()
    plt.xlabel('Wavelength (Å)')
    plt.ylabel('Flux Density (erg s-1 cm-2 Å-1)')

try:
    _ = qspec # just to check tat qspec is defined
    setup_plot()
    plt.plot(wavegrid, qspec, color='k')
    plt.title('quiescent')
except NameError:
    pass

setup_plot()
plt.plot(wavegrid, fspec, color='C3', label='flare')
plt.title('flare')


#%% optional save quiescent spectrum in an ETC uploadable format

qdata = np.vstack((wavegrid.to_value('AA'),
                   qspec.to_value('erg s-1 cm-2')))
filepath = 'ad_leo_quiescent_spectrum_for_etc.dat'
np.savetxt(filepath, qdata.T)

#%% save flare spectrum in an ETC uploadable format

fdata = np.vstack((wavegrid.to_value('AA'),
                   fspec.to_value('erg s-1 cm-2')))
filepath = 'ad_leo_flare_spectrum_for_etc.dat'
np.savetxt(filepath, fdata.T)