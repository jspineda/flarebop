from astropy import units as u

# cumulative FFD for rate of flares at least as big as amplitude
class FFD(object):
    def __init__(self, alpha, rate_constant, pivot_amplitude):
        self.alpha = alpha
        self.rate_constant = rate_constant
        self.pivot_amplitude = pivot_amplitude

    def amplitude(self, rate):
        ratio = rate/self.rate_constant
        eqd = (ratio.to_value(''))**(-1/self.alpha) * self.pivot_amplitude
        return eqd.to('s')

    def rate(self, amplitude):
        ratio = amplitude/self.pivot_amplitude
        rate = self.rate_constant * (ratio.to_value(''))**-self.alpha
        return rate


## add class for differential FFD object (alpha above is not same as alpha below...)
class dFFD(object):
    def __init__(self,alpha,delt_min,Norm=1):
        self.alpha = alpha
        self.d_m = delt_min
        self.norm = Norm    # units of number of flares, or flares per time = number of observed flares/ total time

    def _set_norm(self,norm):
        "To change normalization from pdf to differential FFD, set to unity for proper pdf, amplitude above delt_min"
        self.norm = norm

    def dist(self,amplitude):
        ratio = amplitude / self.d_m
        dist= (self.alpha - 1 )*( ratio.to_value('') )**-self.alpha  / self.d_m
        return dist * self.norm
    
  
  ## add in here Jim's flare evolution model?
  
  
    
# result from Loyd+ 2018b analysis with more sig figs
# to mitigate numerical errors for large extrapolations
fuv130 = FFD(
    alpha = 0.7607041138106101,
    rate_constant = 10**(0.5653602074246712)/u.d, ## i think this should be 10**0.56 ....
    pivot_amplitude = 1000*u.s
)

totflarerate = fuv130.rate(47*u.s) # rhat

# differential FFD based on Loyd+ 2018; set
# Norm = 37.62710710560971 / u.d  # based on delt_min = 47 * u.s ,
# in order to make dFFD consistent with cumulative FFD fuv130
#

dfuv130 = dFFD(
    alpha = 1.7607041138106101, delt_min = 47*u.s
)
