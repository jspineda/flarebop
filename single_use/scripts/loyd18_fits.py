from matplotlib import pyplot as plt
import numpy as np
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u

#%%
cat = Table.read('single_use/data/loyd18_flares.ecsv')
cat.sort('PEW') # equivalent duration (photometric equivalent width)


#%%

fits = Table(names='x xband xunit y yband yunitt alpha beta pivot covariance'.split(),
             dtype='str str str str str str float float float object'.split(),
             masked=True)
fits.meta['notes'] = ('All fits are of the form y = beta * (x / pivot) ** alpha. The covariance column gives the '
                      'covariance matrix for alpha and beta from the fit.')
fits['alpha'].format = '.3g'
fits['beta'].format = '.3g'
fits['pivot'].format = '.3g'


#%% fit to the FUV 130 FFD
"""
This fitting process was complex, so I am not reproducing it here, simply entering a higher precision result
from the values I stored for the original paper.
"""
pivotffd = 1000
alphaffd = 0.7607041138106101
betaffd = 0.5653602074246712
row = 'equiv dur, fuv130, s, nu, fuv130, d-1'.split(', ')
row += [alphaffd, betaffd, pivotffd, None]
mask = np.zeros(len(row), bool)
mask[-1] = True
fits.add_row(row, mask)


#%% fitting function
def fit(line):
    log = np.log10
    y130, e130 = cat['PEW'], cat['PEW_err']
    y, e = cat[f'peak_ratio{line}'], cat[f'peak_ratio_err-line']
    use = (y4 / e4 > 3) & (y / e > 3)
    if np.sum(use) < 5:
        raise ValueError
    y4, e4, y, e = [x[use] for x in (y4, e4, y, e)]
    ly4, ly = list(map(log, (y4, y)))
    l4e, le = [ee/yy/np.log(10) for ee,yy in ((e4,y4), (e,y))]
    assert np.all(le > 0)

    def line_func(p, x):
        m, b = p
        return m*x + b
    # model = odr.Model(line_func)
    # data = odr.RealData(ly4, ly, sx=l4e, sy=le)
    # coeffs = np.polyfit(ly4, ly, 1)
    # fitter = odr.ODR(data, model, beta0=coeffs)
    # out = fitter.run()
    # cov = out.cov_beta
    # if return_cov:
    #     return out.beta, out.cov_beta
    # rho = cov[1,0]/np.sqrt(cov[0,0]*cov[1,1])
    # std = np.std(ly - line_func(out.beta, ly4))
    # return out.beta, out.sd_beta, rho, std

    use = [np.ones(len(ly), bool)]
    while True:
        coeffs, cov = np.polyfit(ly4[use[-1]], ly[use[-1]], 1, cov=True)
        oc = ly - line_func(coeffs, ly4)
        std = np.std(oc)
        use.append(oc/std < 3.)
        if np.all(use[-1] == use[-2]):
            break
    if return_cov:
        return coeffs, cov
    else:
        rho = cov[1,0]/np.sqrt(cov[0,0]*cov[1,1])
        return coeffs, np.sqrt(np.diag(cov,0)), rho, std


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
    eqn = f'{ylbl} = {10**b:.5g}(eqd/{xpivot} s)^{m:.5f}'.format(ylbl, 10 ** b, m)
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

keys = 'si4 fuv130 c3 si3 lya_broad n5 o1 c2 si4 c4 he2'.split()
labels = dict(zip(keys, keys))
labels['lya_broad'] = 'lya wings'
eqd130 = cat['']

for key in keys:
