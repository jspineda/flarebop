from math import nan

from matplotlib import pyplot as plt
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


def get_data(xparam, xband, yparam, yband):
    sets = ((xparam, xband),
            (yparam, yband))
    vals, errs = [], []
    for param, band in sets:
        val = cat[get_colname(param, band, False)]
        errname = get_colname(param, band, True)
        err = cat[errname] if errname in cat.colnames else None
        vals.append(val)
        errs.append(err)

    return *vals, *errs


#%% fitting function

def perform_fit(xparam, xband, yparam, yband, quantile=None):
    x, y, xerr, yerr = get_data(xparam, xband, yparam, yband)
    log = np.log10

    # mask low snr data
    use = ~np.isnan(x.filled(nan)) & ~np.isnan(y.filled(nan))
    if xerr is not None:
        use &= x.filled(0) / xerr.filled(1) > 3
    if yerr is not None:
        use &= y.filled(0) / yerr.filled(1) > 3

    with np.errstate(invalid='ignore', divide='ignore'):
        lx, ly = list(map(log, (x, y)))

    if quantile:
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
            oc = ly - np.polyval(coeffs, lx)
            std = np.std(oc[use])
            masks.append(oc.filled(0)/std < 3.)
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

fit_defintion_order = 'xparam, xband, yparam, yband, quantile'.split(', ')
fits = [
    ['equiv dur', 'fuv130', 'contrast', 'fuv130', None],
    ['equiv dur', 'fuv130', 'contrast', 'fuv130', 0.5],
    ['equiv dur', 'fuv130', 'contrast', 'fuv130', 0.95],
    ['equiv dur', 'fuv130', 'cmltv fwhm', 'fuv130', None],
    ['equiv dur', 'fuv130', 'cmltv fwhm', 'fuv130', 0.5],
    ['equiv dur', 'fuv130', 'cmltv fwhm', 'fuv130', 0.95],
    ['contrast', 'fuv130', 'cmltv fwhm', 'fuv130', None],
    ['contrast', 'fuv130', 'cmltv fwhm', 'fuv130', 0.5],
    ['contrast', 'fuv130', 'cmltv fwhm', 'fuv130', 0.95],
]
lines = 'si4, c3, si3, lya wings, n5, o1, c2, si4, c4, he2'.split(', ')
for line in lines:
    line_fits = [
        ['equiv dur', 'fuv130', 'contrast', line, None],
        ['equiv dur', 'fuv130', 'contrast', line, 0.95],
        ['contrast', 'fuv130', 'contrast', line, None],
        ['contrast', 'fuv130', 'contrast', line, 0.95],
    ]
    fits.extend(line_fits)


#%% carry out fits and fill in the table

fitparams = 'alpha beta pivot npts scatter correlation covariance'.split()
fitrows = []
maskrows = []
for fit in fits:
    fitrow = dict(zip(fit_defintion_order, fit))
    result = perform_fit(**fitrow)
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

#%% fit to the FUV 130 FFD
"""
This fitting process was complex, so I am not reproducing it here. I'm simply entering a higher precision result
from the values I stored from the Loyd+ 2018b analysis.
"""
row = dict(
    xparam = 'equiv dur',
    xband = 'fuv130',
    xunit = 's',
    yparam = 'rate',
    yband = 'fuv130',
    yunit = 'd-1',
    quantile = np.ma.masked,
    alpha = 0.7607041138106101,
    beta = 0.5653602074246712,
    pivot = 1000,
    npts = np.ma.masked
)
fits.add_row(row)

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


#%% plotting function

def plot_fit(xparam, xeqn, yparam, yeqn, row, fitmask):
    x, xerr, y, yerr = get_data(xparam, row['xband'], yparam, row['yband'])

    fig, ax = plt.subplots(1, 1, figsize=[3, 2.5])
    plt.xscale('log')
    plt.yscale('log')

    sets = (
        ('inactive', dict(marker='.', color='C0')),
        ('active', dict(marker='d', color='C1')),
        ('40Myr', dict(marker='s', color='C2'))
    )
    for sample, kws in sets:
        splmask = cat['sample'] == sample
        plt.errorbar(x[splmask], y[splmask], yerr[splmask], xerr[splmask], alpha=0.5, ls='none', **kws)
        mask = splmask & fitmask
        plt.errorbar(x[mask], y[mask], yerr[mask], xerr[mask], ls='none', elinewidth=1, **kws)

    xln = np.min(x[fitmask]), np.max(x[fitmask])
    xln = np.array(xln)
    a, b, xpivot = row['alpha'], row['beta'], row['pivot']
    lyln = np.polyval([a, b], np.log10(xln/xpivot))
    yln = 10**lyln
    plt.plot(xln, yln, 'k-')

    xunit, yunit = row['xunit'], row['yunit']
    xunit_txt = f'\\ [\\mathrm{{{xunit}}}]' if xunit else ''
    yunit_txt = f'\\ [\\mathrm{{{yunit}}}]' if xunit else ''
    eqn = f'${yeqn} = {b:.2g}{yunit_txt}\\ ({xeqn}/{xpivot:.0f}{xunit_txt})^{a:.2f}$'
    plt.annotate(eqn, xy=(0.02,0.98), xycoords='axes fraction', ha='left', va='top', fontsize='small')

    return fig


#%% save plots





#%% peak flux, fwm, eqd for fuv130

fitkws = dict(xparam='PEW', xband='fuv130',
              yparam='peak_ratio', yband='fuv130', quantile=0.95)
result = fit(**fitkws)
add_row(**fitkws,
        xpband='fuv130', xplbl='equiv dur', xunit='s',
        ypband='fuv130', yplbl='contrast', yunit='')
fig = plot_fit(fitkws['xparam'], '\\delta', fitkws['yparam'], 'F_{pk}/F_q', row=fits[-1], fitmask=result[-1])




#%% line peak flux vs fuv130 eqd


#%% line peak flux vs fuv130 peak flux


#%% fits between peak flux, fwm, and eqd

fig, axs = plt.subplots(1,2, figsize=[8,4])
for ax in axs:
    ax.set_xscale('log')
    ax.set_yscale('log')
fig.supxlabel('Equivalent Durations (s)')

fitcat = cat.copy()
keep = fitcat['PEW'] / fitcat['PEW_err'] > 3.
fitcat = fitcat[keep]

pew = fitcat['PEW']
xln = pew[[0, -1]]

def shapefit(ax, ykey, ylbl, xpivot):
    y = np.log10(fitcat[ykey])
    x = np.log10(pew)
    good = ~np.isnan(x) & ~np.isnan(y)
    ax.plot(pew, fitcat[ykey], '.')
    ax.set_ylabel(ylbl)
    (m, b), cov = np.polyfit(x[good] - np.log10(xpivot), y[good], 1, cov=True)
    eqn = f'{ylbl} = {10**b:.5g}(eqd/{xpivot} s)^{m:.5f}'
    yln = 10 ** b * (xln/xpivot) ** m
    ln, = ax.plot(xln, yln, color='0.6')
    print(eqn)
    print(cov)
    print(cov[0,1]/np.sqrt(cov[0,0]*cov[1,1]))
    return m, b, cov

pivot1 = 251
m1, b1, cov1 = shapefit(axs[0], 'peak_ratio', 'Fp/Fq', pivot1)
row = 'equiv dur, fuv130, s, contrast, fuv130, '.split(', ') + [m1, b1, pivot1, cov1]
fits.add_row(row)

pivot2 = 255
m2, b2, cov2 = shapefit(axs[1], 'FWHM', 'FWHM', pivot2)
row = 'equiv dur, fuv130, s, fhwm, fuv130, s'.split(', ') + [m2, b2, pivot2, cov2]
fits.add_row(row)


#%% line peak ratio fits


eqd130 = cat['']

for key in linekeys:
