import numpy as np
from astropy import table
Table = table.Table

np.set_printoptions(precision=3)


#%% load in flare catalogs

Mcat = Table.read('fits/data/loyd18b_flares_inactive.ecsv')
Fcat = Table.read('fits/data/loyd18b_flares_active.ecsv')
Hcat = Table.read('fits/data/loyd18a_flares_40Myr_fuv130.ecsv')
Mcat['sample'] = 'inactive'
Fcat['sample'] = 'active'
Hcat['sample'] = '40Myr'
cat = table.vstack((Mcat, Fcat, Hcat))

#%% standardize column names so there aren't special cases to handle later
replacement_map = {
    'lya_broad' : 'lya wings',
    'PEW' : 'equiv dur',
    'peak_ratio' : 'contrast',
    'FWHM' : 'cmltv fwhm',
}
for colname in cat.colnames:
    newname = colname
    if '-' not in colname and colname != 'sample':
        newname = f'{newname}-fuv130'
    for key, val in replacement_map.items():
        newname = newname.replace(key, val)
    cat.rename_column(colname, newname)


#%% correct contrast values
# original constrast values were Fpk/Fq + 1. The +1 was a mistake. This work will define contrast as
# Fpk/Fq - 1, so that the minimum flare contrast is 0.
#
# this also messed up the uncertainties. Some math to correct those too:
# r = p/q + 1 %original contrast definition, r is contrast, p peak flux, q quiescent flux
#
# s_r = r * √((s_p/p)2 + (s_q/q)2)  %what Loyd 2018a,b used to compute uncty, wrong in several ways (sigh)
#
# r' = p/q - 1 % how we're defining contrast here
# s_r' = √((s_p/q)2 + (s_q*p/q2)2) = p/q * √((s_p/p)2 + (s_q/q)2) = (r-1) * s_r/r = (r-1)/r * s_r

contrast_names = [name for name in cat.colnames if 'contrast' in name and 'err' not in name]
for name in contrast_names:
    r = cat[name]
    ename = name.replace('contrast', 'contrast_err')
    e = cat[ename]
    rr = r - 2
    ee = (r - 1) / r * e
    cat[name] = rr
    cat[ename] = ee


#%% sort and save
cat.sort('equiv dur-fuv130') # equivalent duration (photometric equivalent width)

cat.write('fits/data/all_flares_compiled.ecsv', overwrite=True)