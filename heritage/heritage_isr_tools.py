import numpy as np
from astropy import units as u
from astropy.modeling.models import Gaussian1D

import synphot

#%% tools to generate spectra for input to etc

def gaussian_amp(flux, sigma):
    return flux / np.sqrt(2 * np.pi) / sigma

def sig_from_fwhm_wave(fwhm):
    sig = fwhm/2/np.sqrt(2*np.log(2))
    return sig.to(fwhm.unit)

def flare_line(w, w0, flux, fwhm):
    sigma = sig_from_fwhm_wave(fwhm)
    amp = gaussian_amp(flux, sigma)
    gauss = Gaussian1D(amp, w0, sigma)
    return gauss(w)

vegaspec = synphot.SourceSpectrum.from_vega()
u_band = synphot.SpectralElement.from_filter('johnson_u')
def normed_bb(w, Umag):
    Tbb = 9000*u.K
    bb = synphot.models.BlackBody1D(temperature=Tbb)
    bb = synphot.SourceSpectrum(bb)
    bb_normed = bb.normalize(Umag*synphot.units.VEGAMAG, band=u_band, vegaspec=vegaspec)
    y = bb_normed(w, flux_unit='FLAM')
    return y.to('erg s-1 cm-2 AA-1')

def heritage_isr_spec(wavegrid, Fc4, Fsi4, Flya, Uflare):
    w = wavegrid
    yc4 = flare_line(w, 1548.2*u.AA, Fc4, 0.2*u.AA)
    ysi4 = flare_line(w, 1393.8*u.AA, Fsi4, 0.2*u.AA)
    ylya = flare_line(w, 1215.7*u.AA, Flya, 0.5*u.AA)
    ybb = normed_bb(w, Uflare)
    y = ybb + yc4 + ysi4 + ylya
    return y