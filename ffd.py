from astropy import units as u


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

# result from Loyd+ 2018b analysis with more sig figs
# to mitigate numerical errors for large extrapolations
fuv130 = FFD(
    alpha = 0.7607041138106101,
    rate_constant = 0.5653602074246712/u.d,
    pivot_amplitude = 1000*u.s
)