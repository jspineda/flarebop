
from scipy import special
from scipy.stats import binned_statistic
from scipy.interpolate import make_interp_spline
from scipy.special import erf
import numpy as np
import matplotlib.pyplot as plt

from astropy import units as u
import astropy.constants as const

import ffd
import pdb

### from Mendoza et al. 2022 (based on broadband optical flares):

def flare_eqn(t,ampl):
    '''
    The equation that defines the shape for the Continuous Flare Model
    '''
    #Values were fit & calculated using MCMC 256 walkers and 30000 steps

    A,B,C,D1,D2,f1 = [0.9687734504375167,-0.251299705922117,0.22675974948468916,
                      0.15551880775110513,1.2150539528490194,0.12695865022878844]

    # We include the corresponding errors for each parameter from the MCMC analysis

    A_err,B_err,C_err,D1_err,D2_err,f1_err = [0.007941622683556804,0.0004073709715788909,0.0006863488251125649,
                                              0.0013498012884345656,0.00453458098656645,0.001053149344530907 ]

    f2 = 1-f1

    eqn = ((1 / 2) * np.sqrt(np.pi) * A * C * f1 * np.exp(-D1 * t + ((B / C) + (D1 * C / 2)) ** 2)
                        * special.erfc(((B - t) / C) + (C * D1 / 2))) + ((1 / 2) * np.sqrt(np.pi) * A * C * f2
                        * np.exp(-D2 * t+ ((B / C) + (D2 * C / 2)) ** 2) * special.erfc(((B - t) / C) + (C * D2 / 2)))

    return eqn * ampl

### from Mendoza et al. 2022 (based on broadband optical flares):
def flare_model(t,tpeak, fwhm, ampl, upsample=False, uptime=10):
    '''
    The Continuous Flare Model evaluated for single-peak (classical) flare events.
    Use this function for fitting classical flares with most curve_fit
    tools.

    References
    --------------
    Davenport et al. (2014) http://arxiv.org/abs/1411.3723
    Jackman et al. (2018) https://arxiv.org/abs/1804.03377

    Parameters
    ----------
    t : 1-d array
        The time array to evaluate the flare over

    tpeak : float
        The center time of the flare peak

    fwhm : float
        The Full Width at Half Maximum, timescale of the flare

    ampl : float
        The amplitude of the flare


    Returns
    -------
    flare : 1-d array
        The flux of the flare model evaluated at each time

        A continuous flare template whose shape is defined by the convolution of a Gaussian and double exponential
        and can be parameterized by three parameters: center time (tpeak), FWHM, and ampitude
    '''

    t_new = (t-tpeak)/fwhm

    if upsample:
        dt = np.nanmedian(np.diff(np.abs(t_new)))
        timeup = np.linspace(min(t_new) - dt, max(t_new) + dt, t_new.size * uptime)

        flareup = flare_eqn(timeup,ampl)

        # and now downsample back to the original time...

        downbins = np.concatenate((t_new - dt / 2.,[max(t_new) + dt / 2.]))
        flare,_,_ = binned_statistic(timeup, flareup, statistic='mean',bins=downbins)
    else:

        flare = flare_eqn(t_new,ampl)

    return flare

flare_peak = flare_model(0,0,1,1) # value for profile comparison


## below is for working with profile and computing flare durations at given threshold -- may not be necessary?
class FlareThresh():
    ""
    
    def __init__(self,starclass):
        ""

        ## pulls thresholds from starclass, then sets up for calculations
        self.__setup__()
        
        ##test:

    def __setup__(self):
        ""
        ## Loading data for initial det up
        ## flare durations
        fldat = np.genfromtxt("./fits/output tables/flare_prof_dur.dat",skip_header=2,delimiter=',')
        cent = np.where(fldat[:,0] >= 0.5)[0][0] # central pivot value
        
        lowp = make_interp_spline( np.log(fldat[:-1,0]), np.log(fldat[:-1,1]) , k=1)
        highp = make_interp_spline( fldat[:,0] , fldat[:,1] , k=1)
        
        self._fldur_ext = {"Data":fldat, "Pivot":fldat[cent,0], "ExtLow":lowp, "ExtHigh":highp}
        # want to take care of extrapolations beyond initial data, uses relative time units, and relative amplitude
    
    
    def fldur(self,tr_in,fwhm):
        "For computing in relative time the duration that the flare exceeds the input thresholds in relative amplitude"
        
        if np.any(tr_in <= 0):
            raise ValueError(f'Input thresholds are negative, '
                             f'quiescence exceeds rate limits, check inputs..'  )
        
        reldur = np.zeros_like(tr_in)
    
        # pivot change for interpolation regimes
        indhigh = np.where(tr_in >= self._fldur_ext["Pivot"] )
        indlow = np.where(tr_in < self._fldur_ext["Pivot"] )
        
        reldur[indlow] = np.exp( self._fldur_ext["ExtLow"]( np.log(tr_in[indlow]) ) )
        reldur[indhigh] = self._fldur_ext["ExtHigh"](tr_in[indhigh])

        iover = np.where(tr_in >= flare_peak)[0]

        if len(iover)>0:
            reldur[iover] = 0
    
        return reldur * fwhm

