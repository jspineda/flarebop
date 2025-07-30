from matplotlib import pyplot as plt
import numpy as np
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u

from settings import risk_tolerance, characteristic_duration


def amplitude(rate, alpha, beta, pivot):
    eqd = (rate/beta)**(-1/alpha) * pivot
    return eqd.to('s')


def rate(amplitude, alpha, beta, pivot):
    rate = beta * (amplitude/pivot)**-alpha
    return rate


