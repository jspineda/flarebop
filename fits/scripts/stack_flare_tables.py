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

# standardize column names so there aren't special cases to handle later
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


cat.sort('equiv dur-fuv130') # equivalent duration (photometric equivalent width)

cat.write('fits/data/all_flares_compiled.ecsv', overwrite=True)