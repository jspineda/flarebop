import numpy as np
from astropy import units as u
from astropy import constants as const
from astropy import coordinates as coord
from scipy.special import wofz


Tism = 10000*u.K
wlab_H = 1215.67*u.AA
wlab_D = 1215.34*u.AA
D_H = 1.5e-5  # from bourrier+ 2017 (kepler 444)
f = 4.1641e-01
A = 6.2649e8 / u.s
mH, mD = 1 * u.u, 2 * u.u


def local_ism_rv(ra, dec):
    # see http://lism.wesleyan.edu/ColoradoLIC.html and
    # lallement & bertin 1992
    v_sun = 26*u.km/u.s
    v_sun_direction = coord.SkyCoord(74.5*u.deg, -6*u.deg)
    v_sun_vector = v_sun * v_sun_direction.represent_as('cartesian').xyz

    target_directions = coord.SkyCoord(ra, dec)
    target_vectors = target_directions.represent_as('cartesian').xyz

    v_sun_component = np.dot(v_sun_vector.to_value(), target_vectors) * v_sun_vector.unit
    return v_sun_component


def doppler_shift(w, velocity):
    return (1 + velocity/const.c)*w


def voigt_xsection(w, w0, f, gamma, T, mass, b=None):
    """
    Compute the absorption cross section using hte voigt profile for a line.

    Parameters
    ----------
    w : astropy quantity array or scalar
        Scalar or vector of wavelengths at which to compute cross section.
    w0: quanitity
        Line center wavelength.
    f: scalar
        Oscillator strength.
    gamma: quantity
        Sum of transition rates (A values) out of upper and lower states. Just Aul for a resonance line where only
        spontaneous decay is an issue.
    T: quantity
        Temperature of the gas. Can be None if you provide a b value instead.
    mass: quantity
        molecular mass of the gas
    b : quantity
        Doppler b value (in velocity units) of the line
    Returns
    -------
    x : quantity
        Cross section of the line at each w.
    """

    nu = const.c / w
    nu0 = const.c / w0
    if T is None:
        sigma_dopp = b/const.c*nu0/np.sqrt(2)
    else:
        sigma_dopp = np.sqrt(const.k_B*T/mass/const.c**2) * nu0
    dnu = nu - nu0
    gauss_sigma = sigma_dopp.to(u.Hz).value
    lorentz_FWHM = (gamma/2/np.pi).to(u.Hz).value
    phi = voigt(dnu.to(u.Hz).value, gauss_sigma, lorentz_FWHM) * u.s
    x = np.pi*const.e.esu**2/const.m_e/const.c * f * phi
    return x.to('cm2')


def voigt(x, gauss_sigma, lorentz_FWHM):
    """
    Return the Voigt line shape at x with Lorentzian component HWHM gamma
    and Gaussian component HWHM alpha.

    """
    sigma = gauss_sigma
    gamma = lorentz_FWHM/2.0
    return np.real(wofz((x + 1j*gamma)/sigma/np.sqrt(2))) / sigma /np.sqrt(2*np.pi)


def transmission(w, rv, Nh, T):
    w0s = doppler_shift((wlab_H, wlab_D)*u.AA, rv)
    xsections = [voigt_xsection(w, w0, f, A, T, m) for
                 w0, m in zip(w0s, (mH, mD))]
    tau = xsections[0]*Nh + xsections[1]*Nh*D_H
    return np.exp(-tau)