from astropy import units as u

characteristic_duration = 1*u.h
risk_tolerance = 1e-4

# result from Loyd+ 2018b analysis with more sig figs
# to mitigate numerical errors for large extrapolations
fuv130_ffd = dict(
    alpha = 0.7607041138106101,
    beta = 0.5653602074246712,
    pivot = 1000
)