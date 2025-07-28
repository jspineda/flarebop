from matplotlib import pyplot as plt
import numpy as np
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u

from settings import risk_tolerance, characteristic_duration


def eqd(rate):
    eqd = (rate/beta)**(-1/alpha) * eqd_pivot
    return eqd.to('s')


def rate(eqd):
    rate = beta * (eqd/eqd_pivot)**-alpha
    return rate


def flare_constrast(eqd):
    Fp_Fq = 0.34 * eqd.to_value('s')**0.59
    return Fp_Fq


