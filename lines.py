from math import pi

from astropy import table
Table = table.Table
import numpy as np
from astropy.modeling.models import Voigt1D
from astropy import units as u
from astropy import constants as const
from astropy.coordinates import SkyCoord

import worst_case
import ism

correlations = Table.read('tables/quiescent_line_fits.ecsv')
correlations.add_index('xband')
correlations.add_index('yband')

line_properties = Table.read('tables/fuv_line_list.ecsv')
linekeys = np.unique(line_properties['qkey']).tolist()

def fill_quiescent_lines(**input_line_fluxes):
    fluxdict = input_line_fluxes.copy()
    input_lines = list(fluxdict.keys())
    xlines = input_lines[:]
    if 'OI' in xlines:
        xlines.remove('OI')
    xcorrelations = correlations.loc['xband', xlines]

    to_fill = set(linekeys) - set(input_lines)
    for line in to_fill:
        # FIXME once we have O I from Fernando
        if line == 'OI':
            not_lya = set(input_lines) - {'Lya'}
            not_lya = list(not_lya)
            fluxdict[line] = fluxdict[not_lya[0]]
            continue

        ycorrelations = xcorrelations.loc['yband', line]

        # pick the correlation with the smallest scatter
        ipick = np.argmin(ycorrelations['scatter'])
        correlation = ycorrelations[ipick]

        # calculate the flux of the line to fill
        x = fluxdict[correlation['xband']]
        x = x.to_value('erg s-1 cm-2')
        lx = np.log10(x)
        alpha = correlation['alpha']
        beta = correlation['beta']
        pivot = correlation['pivot']
        ly = alpha * (lx - pivot) + beta
        y = 10**ly * u.Unit('erg s-1 cm-2')

        fluxdict[line] = y

    return fluxdict


def fwhm_wave_from_velocity(w0, fwhm_velocity):
    fwhm_wave =  fwhm_velocity / const.c * w0
    return fwhm_wave.to('AA')


def voigt_profile(wavegrid, w0, A, fwhm_velocity):
    fwhm_l = w0**2/(2*pi)/const.c * A
    fwhm_g = fwhm_wave_from_velocity(w0, fwhm_velocity)
    voigt = Voigt1D(
        x_0=w0,
        amplitude_L=1,
        fwhm_L=fwhm_l,
        fwhm_G=fwhm_g
    )
    y = voigt(wavegrid)
    y /= np.max(y)
    y = y.to_value('') / wavegrid.unit
    return y


def standard_line_profiles(wavegrid, line_table):
    profiles = []
    for line in line_table:
        w0 = line['wave'] * u.AA
        fwhm_n_kms = line['FWHM narrow'] * u.km / u.s
        A = line['A'] / u.s
        f = voigt_profile(wavegrid, w0, A, fwhm_n_kms)
        if not np.ma.is_masked(line['FWHM broad']):
            fwhm_b_kms = line['FWHM broad'] * u.km / u.s
            ratio = line['broad narrow ratio']
            f += ratio * voigt_profile(wavegrid, w0, A, fwhm_b_kms)
        f *= line['component intensity']
        profiles.append(f)
    return profiles


def spectrum(wavegrid, **line_fluxes):
    """

    Parameters
    ----------
    wavegrid : array
        Wavelength grid with units.
    line_fluxes :
        Line keys and fluxes given as keywords, such as SiIV=10, CII=3, etc. Units optional.

    Returns
    -------

    """
    inrange = (line_properties['wave'] > wavegrid[0]) & (line_properties['wave'] < wavegrid[-1])
    linetable = line_properties.copy()[inrange]

    # check that fluxes are given for all lines in the wavegid range
    input_lines = set(line_fluxes.keys())
    inrange_lines = set(linetable['qkey'])
    missing_lines = inrange_lines - input_lines
    if missing_lines:
        raise ValueError(f'No fluxes provided for these lines that are in the '
                         f'{wavegrid[0]}–{wavegrid[-1]} range of the wave grid:'
                         f'\n\t{missing_lines} ')

    # get profiles for all line components
    profiles = standard_line_profiles(wavegrid, linetable)
    linetable['profiles'] = profiles # units get lost here
    profile_unit = profiles[0].unit
    linetable.add_index('qkey')

    # group by line, normalize to the input flux, and sum
    spec = np.zeros_like(wavegrid.value)
    for line in inrange_lines:
        components = linetable.loc[line]
        fs = np.vstack(components['profiles'])
        axis = 1 if fs.shape[0] == len(wavegrid) else 0
        f = np.sum(fs, axis=axis) * profile_unit
        F = np.trapz(f, wavegrid)
        input_flux = line_fluxes[line]
        normfac = input_flux / F
        f = f * normfac
        spec = spec + f

    return spec


def automated_line_spectrum(
        wavegrid,
        known_fluxes,
        ra,
        dec,
        distance,
        rv_star=None,
        rv_ism='local ism',
        Nh=17,
):

    coords = SkyCoord(ra, dec)

    dist_to_1AU = (distance / u.au) ** 2
    dist_to_1AU = dist_to_1AU.to_value('')
    line_fluxes_1AU = {line: flux * dist_to_1AU for line, flux in known_fluxes.items()}
    all_line_fluxes_1AU = fill_quiescent_lines(**line_fluxes_1AU)
    all_line_fluxes = {line: flux / dist_to_1AU for line, flux in all_line_fluxes_1AU.items()}

    fluxes_worst_case = {}
    for line, quiescent_flux in all_line_fluxes.items():
        flare_contrast = worst_case.line_contrasts[line]
        flare_flux = quiescent_flux * flare_contrast
        fluxes_worst_case[line] = flare_flux
    fspec = spectrum(wavegrid, **fluxes_worst_case)

    if rv_star is None:
        print('No ISM absorption added because stellar RV not provided.')
    else:
        if rv_ism == 'local ism':
            rv_ism = ism.local_ism_rv(coords)
        rv_diff = rv_ism - rv_star
        transmission = ism.transmission(wavegrid, rv_diff, Nh.to('cm-2'))
        fspec *= transmission

    return fspec