from math import nan
from pathlib import Path

from matplotlib import pyplot as plt
plt.ion()
import numpy as np
from astropy import table
Table = table.Table
from sklearn.linear_model import QuantileRegressor

np.set_printoptions(precision=3)

#%% units

unit_map = {
    'equiv dur' : 's',
    'contrast' : '',
    'cmltv fwhm' : 's'
}

#%% load in flare catalogs

cat = Table.read('fits/data/all_flares_compiled.ecsv')

for colname in cat.colnames:
    if not hasattr(cat[colname], 'mask'):
        cat[colname] = table.MaskedColumn(cat[colname])


#%% data retrieval functions

def get_colname(parameter, band, err=False):
    if err:
        parameter += '_err'
    return '-'.join((parameter, band))


def get_data(fitdef):
    sets = ((fitdef['xparam'], fitdef['xband']),
            (fitdef['yparam'], fitdef['yband']))
    vals, errs = [], []
    for param, band in sets:
        val = cat[get_colname(param, band, False)]
        errname = get_colname(param, band, True)
        err = np.abs(cat[errname]) if errname in cat.colnames else None
        # the np.abs() fixes a mistake from flare papers where I used an error estimate for a ratio of
        # error = a/b * sqrt((err_a/a)**2 + (err_b/b)**2) where a/b should have been abs(a/b)
        vals.append(val)
        errs.append(err)

    return *vals, *errs


#%% fitting function

def perform_fit(fitdef):
    x, y, xerr, yerr = get_data(fitdef)
    log = np.log10

    # mask low snr data
    use = ~np.isnan(x.filled(nan)) & ~np.isnan(y.filled(nan))
    if xerr is not None:
        use &= x.filled(0) / xerr.filled(1) > 3
    if yerr is not None:
        use &= y.filled(0) / yerr.filled(1) > 3

    with np.errstate(invalid='ignore', divide='ignore'):
        lx, ly = list(map(log, (x, y)))

    if fitdef['percentile'] != 'n/a':
        quantile = float(fitdef['percentile'])/100
        qr = QuantileRegressor(quantile=quantile, alpha=0)
        mask = use
        lxpivot = np.mean(lx[mask])
        xfit = lx[mask] - lxpivot
        yfit = ly[mask]
        qr.fit(xfit[:,None], yfit)
        coeffs = qr.coef_[0], qr.intercept_
        oc = ly - np.polyval(coeffs, lx)
        std = np.std(oc)
        cov = None
        rho = None
    else:
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

fit_defintion_order = 'xparam, xband, yparam, yband, percentile'.split(', ')
fits = [
    ['equiv dur', 'fuv130', 'contrast', 'fuv130', 'n/a'],
    ['equiv dur', 'fuv130', 'contrast', 'fuv130', '50'],
    ['equiv dur', 'fuv130', 'contrast', 'fuv130', '95'],
    ['equiv dur', 'fuv130', 'cmltv fwhm', 'fuv130', 'n/a'],
    ['equiv dur', 'fuv130', 'cmltv fwhm', 'fuv130', '50'],
    ['equiv dur', 'fuv130', 'cmltv fwhm', 'fuv130', '95'],
    ['contrast', 'fuv130', 'cmltv fwhm', 'fuv130', 'n/a'],
    ['contrast', 'fuv130', 'cmltv fwhm', 'fuv130', '50'],
    ['contrast', 'fuv130', 'cmltv fwhm', 'fuv130', '95'],
]
lines = 'c3, si3, lya wings, n5, o1, c2, si4, c4, he2'.split(', ')
for line in lines:
    line_fits = [
        ['equiv dur', 'fuv130', 'contrast', line, 'n/a'],
        ['equiv dur', 'fuv130', 'contrast', line, '95'],
        ['contrast', 'fuv130', 'contrast', line, 'n/a'],
        ['contrast', 'fuv130', 'contrast', line, '95'],
    ]
    fits.extend(line_fits)


#%% carry out fits and fill in the table

fitparams = 'alpha beta pivot npts scatter correlation covariance'.split()
fitrows = []
maskrows = []
for fit in fits:
    fitrow = dict(zip(fit_defintion_order, fit))
    result = perform_fit(fitrow)
    for param in fitparams:
        if result[param] is not None:
            fitrow[param] = result[param]
    fitrow['xunit'] = unit_map[fitrow['xparam']]
    fitrow['yunit'] = unit_map[fitrow['yparam']]
    fitrows.append(fitrow)

    maskrow = dict(zip(fit_defintion_order, fit))
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

#%% save fit and mask tables

for name in fits.colnames:
    mask = fits[name] == None
    try:
        mask |= np.isnan(fits[name].filled(nan))
    except:
        pass
    fits[name].mask = mask

fits.write('fits/output tables/flare_line_fits.ecsv', overwrite=True)
masks.write('fits/output tables/flare_line_fits_point_masks.ecsv', overwrite=True)


#%% prep to make plots

fits = Table.read('fits/output tables/flare_line_fits.ecsv')
masks = Table.read('fits/output tables/flare_line_fits_point_masks.ecsv')

# setup to be able to find the right mask
loc_keys = ('xparam', 'xband', 'yparam', 'yband', 'percentile')
for key in loc_keys:
    masks.add_index(key)
def get_mask(fitdef):
    slim = masks[:]
    for key in loc_keys:
        slim = slim.loc[key, fitdef[key]]
    return ~slim['mask']

def makelabel(fitdef, axis='x'):
    param = fitdef[f'{axis}param']
    band = fitdef[f'{axis}band']
    unit = fitdef[f'{axis}unit']
    label = f'{param} - {band}'
    if unit:
        label += f' ({unit})'
    return label


def safe_errplot(x, xerr, y, yerr, mask, **pltkws):
    x, y = x[mask], y[mask]
    if xerr is not None:
        xerr = xerr[mask]
    if yerr is not None:
        yerr= yerr[mask]
    return plt.errorbar(x, y, yerr, xerr, **pltkws)


#%% make plots

folder = Path('fits/plots/flare line fits/')

for fitdef in fits:
    x, y, xerr, yerr = get_data(fitdef)
    fitmask = get_mask(fitdef)

    fig, ax = plt.subplots(1, 1, figsize=[4, 3])
    plt.xlabel(makelabel(fitdef, 'x'))
    plt.ylabel(makelabel(fitdef, 'y'))

    sets = (
        ('inactive', dict(marker='.', color='C0')),
        ('active', dict(marker='d', color='C1')),
        ('40Myr', dict(marker='s', color='C2'))
    )
    for sample, kws in sets:
        splmask = cat['sample'].filled('') == sample
        safe_errplot(x, xerr, y, yerr, splmask, alpha=0.5, ls='none', **kws)
        mask = splmask & fitmask
        safe_errplot(x, xerr, y, yerr, mask, ls='none', elinewidth=1, **kws)

    plt.xscale('log')
    plt.yscale('log')

    xln = np.min(x[fitmask]), np.max(x[fitmask])
    xln = np.array(xln)
    a, b, xpivot = fitdef['alpha'], fitdef['beta'], fitdef['pivot']
    lyln = np.polyval([a, b], np.log10(xln/xpivot))
    yln = 10**lyln
    plt.plot(xln, yln, 'k-')

    eqn = f'$y = {b:.2g}\\ (x/{xpivot:.0f})^{{{a:.2f}}}$'
    pct = f'percentile = {fitdef['percentile']}'
    plt.annotate('\n'.join((eqn,pct)), xy=(0.02,0.98), xycoords='axes fraction', ha='left', va='top', fontsize='small')
    plt.tight_layout()

    filename = (f'{fitdef['xparam'].replace(' ', '')}-{fitdef['xband']}.vs.'
                f'{fitdef['yparam'].replace(' ', '')}-{fitdef['yband']}')
    if fitdef['percentile'] != 'n/a':
        filename += f'.{fitdef['percentile']}th-percentile'
    filename += '.pdf'
    plt.savefig(folder / filename)

