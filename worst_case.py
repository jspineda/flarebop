from math import log10

from astropy import units as u
from astropy.table import Table
import ffd

characteristic_duration = 1*u.h
risk_tolerance = 1e-4
rate = risk_tolerance / characteristic_duration
rate = rate.to('d-1')
eqd_fuv130 = ffd.fuv130.amplitude(rate)

line_list = Table.read('tables/fuv_line_list.ecsv')

flare_correlations = Table.read('tables/flare_line_fits.ecsv')
_select = (
    (flare_correlations['xparam'] == 'equiv dur') &
    (flare_correlations['xband'] == 'fuv130') &
    (flare_correlations['yparam'] == 'contrast') &
    (flare_correlations['yband'] != 'fuv130') &
    (flare_correlations['percentile'] == '95')
)
flare_correlations = flare_correlations[_select]
flare_correlations.add_index('yband')

# treat all of lya like the wings
_lya_row = flare_correlations.loc['lya wings']
_lya_row['yband'] = 'lya'
flare_correlations.add_row(_lya_row)

key_map = dict(zip(line_list['fkey'], line_list['qkey']))
line_contrasts = {}
for row in flare_correlations:
    fkey = row['yband']
    qkey = key_map.get(fkey, fkey)
    alpha = row['alpha']
    beta = row['beta']
    pivot = row['pivot']
    lx = log10(eqd_fuv130.to_value('s'))
    ly = alpha * (lx - pivot) + beta
    y = 10**ly
    line_contrasts[qkey] = y