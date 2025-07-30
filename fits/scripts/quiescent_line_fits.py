from math import nan
from pathlib import Path

from matplotlib import pyplot as plt
plt.ion()
import numpy as np
from astropy import table
Table = table.Table
from sklearn.linear_model import QuantileRegressor

np.set_printoptions(precision=3)


#%% load in flare catalogs

cat = Table.read('fits/data/all_lines_compiled.ecsv')

for colname in cat.colnames:
    if not hasattr(cat[colname], 'mask'):
        cat[colname] = table.MaskedColumn(cat[colname])


#%% data retrieval functions

def get_data(fitdef):
    xband, yband = fitdef['xband'], fitdef['yband']
    x = cat['F-' + xband]
    y = cat['F-' + yband]
    return x, y


#%% fitting function

def perform_fit(fitdef):
    x, y = get_data(fitdef)
    log = np.log10

    # mask missing data
    use = ~np.isnan(x.filled(nan)) & ~np.isnan(y.filled(nan))
    with np.errstate(invalid='ignore', divide='ignore'):
        lx, ly = list(map(log, (x, y)))

    masks = [np.ones(len(ly), bool)]
    while True:
        mask = masks[-1] & use
        lxpivot = np.mean(lx[mask])
        coeffs, cov = np.polyfit(lx[mask] - lxpivot, ly[mask], 1, cov=True)
        oc = ly - np.polyval(coeffs, lx - lxpivot)
        std = np.std(oc[use])
        masks.append(oc.filled(np.inf)/std < 3.)
        if np.all(masks[-1] == masks[-2]):
            break
    rho = cov[1,0]/np.sqrt(cov[0,0]*cov[1,1])

    alpha, beta = coeffs
    xpivot = 10**lxpivot

    result = dict(alpha=alpha,
                  beta=beta,
                  pivot=xpivot,
                  scatter=std,
                  correlation=rho,
                  covariance=cov,
                  mask=mask,
                  npts=sum(mask))
    return result


#%% define fits we want to produce

lines = 'CIII SiIII Lya NV CII SiIV CIV HeII'.split(' ')
fits = []
for i, xline in enumerate(lines):
    for yline in lines[i+1:]:
        fit = dict(xband=xline, yband=yline)
        fits.append(fit)


#%% carry out fits and fill in the table

fitparams = 'alpha beta pivot npts scatter correlation covariance'.split()
fitrows = []
maskrows = []
for fit in fits:
    result = perform_fit(fit)
    fitrow = fit.copy()
    for param in fitparams:
        if result[param] is not None:
            fitrow[param] = result[param]
    fitrows.append(fitrow)

    maskrow = fit.copy()
    maskrow['mask'] = ~result['mask'] # inverting to define the more conventional way (masked values are *not* used)
    maskrows.append(maskrow)


#%% construct tables

fits = Table(rows=fitrows, masked=True)
fits.meta['notes'] = ('All fits are of the form y = beta * (x / pivot) ** alpha. The covariance column gives the '
                      'covariance matrix for alpha and beta from the fit and correlation gives the correlation coefficient'
                      'between alpha and beta. Correlation values near zero mean the best fit can be extrapolated '
                      'further with less fear of numerical error (but extrapolation is still dangerous from a '
                      'scientific persective). Scatter gives the standard deviation of residuals.')
fits['alpha'].format = '.3g'
fits['beta'].format = '.3g'
fits['pivot'].format = '.3g'
fits['scatter'].format = '.2g'
fits['correlation'].format = '.2g'

masks = Table(rows=maskrows)


#%% add fit for O I from cruz-aguirre

"""I couldn't find a table of line fluxes to recreat this fit. Note that it is for K and M stars.

Some math to convert their Lbol ratios to 1 AU:
C = 1/(4 pi (1 AU)**2)
Fline(1AU) = Lline/Lbol * Lbol / C
log(FOI * C / Lbol) = a * log(Fuv * C / Lbol) + b
log(FOI) + log(C/Lbol) = a * (log(Fuv) + log(C/Lbol)) + b
log(FOI) = a (log(Fuv) + (a - 1) * log(C/Lbol) + b

Ah so it doesn't work. I emailed Fernando to ask if he can send me a table of O I fluxes. 
"""


#%% save fit and mask tables

for name in fits.colnames:
    mask = fits[name] == None
    try:
        mask |= np.isnan(fits[name].filled(nan))
    except:
        pass
    fits[name].mask = mask

fits.write('fits/output tables/quiescent_line_fits.ecsv', overwrite=True)
masks.write('fits/output tables/quiescent_line_fits_point_masks.ecsv', overwrite=True)


#%% prep to make plots

fits = Table.read('fits/output tables/quiescent_line_fits.ecsv')
masks = Table.read('fits/output tables/quiescent_line_fits_point_masks.ecsv')

# setup to be able to find the right mask
loc_keys = ('xband', 'yband')
for key in loc_keys:
    masks.add_index(key)
def get_mask(fitdef):
    slim = masks[:]
    for key in loc_keys:
        if isinstance(slim, table.row.Row):
            break
        slim = slim.loc[key, fitdef[key]]
    return ~slim['mask']


#%% make plots

folder = Path('fits/plots/quiescent line fits/')

for fitdef in fits:
    x, y = get_data(fitdef)
    fitmask = get_mask(fitdef)

    fig, ax = plt.subplots(1, 1, figsize=[4, 3])
    plt.xlabel(fitdef['xband'])
    plt.ylabel(fitdef['yband'])

    plt.plot(x, y, 'o', color='C0', fillstyle='none')
    plt.plot(x[fitmask], y[fitmask], 'o', color='C0')

    plt.xscale('log')
    plt.yscale('log')

    xln = np.min(x[fitmask]), np.max(x[fitmask])
    xln = np.array(xln)
    a, b, xpivot = fitdef['alpha'], fitdef['beta'], fitdef['pivot']
    lyln = np.polyval([a, b], np.log10(xln/xpivot))
    yln = 10**lyln
    plt.plot(xln, yln, 'k-')

    eqn = f'$y = {b:.2g}\\ (x/{xpivot:.2g})^{{{a:.2f}}}$'
    plt.annotate(eqn, xy=(0.02,0.98), xycoords='axes fraction', ha='left', va='top', fontsize='small')
    plt.tight_layout()

    filename = (f'{fitdef['xband']}.vs.{fitdef['yband']}.pdf')
    plt.savefig(folder / filename)

