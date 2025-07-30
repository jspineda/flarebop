import warnings
from pathlib import Path
import re

from astropy import table
from astropy import units as u
import numpy as np
from matplotlib import pyplot as plt
from astroquery.simbad import Simbad

folder = Path('fits/data')


#%% simbad tic id query function

def get_simbad_info(names):
    simquery = Simbad()
    simquery.add_votable_fields('ids', 'plx_value')
    simbad = simquery.query_objects(list(names))
    dist = 1000/simbad['plx_value']
    simbad['dist'] = dist
    return simbad[['main_id', 'dist']]

#%% function for nice table formatting

def format_numeric_cols(tbl, fmt='.2g'):
    for name in tbl.colnames:
        if tbl[name].dtype == float:
            tbl[name].format = fmt

#%% function to join tables and keep values from priority table

def merge_tables(tbl1, tbl2, priority):
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', message='Cannot merge meta')
        merged = table.join(tbl1, tbl2, 'main_id', 'outer')

    for name in merged.colnames:
        if name.endswith('_1'):
            name2 = name[:-2] + '_2'

            if priority == 2:
                replace = (~merged[name].mask) & (~merged[name2].mask)
                merged[name][replace] = merged[name2][replace]

            move = (merged[name].mask) & (~merged[name2].mask)
            merged[name][move] = merged[name2][move]

            merged.rename_column(name, name[:-2])

    remove_names = [name for name in merged.colnames if name.endswith('_2')]
    merged.remove_columns(remove_names)

    return merged


#%% load france tables

france1 = table.Table.read(folder / 'france2018_table1_SpTs.csv')
france2 = table.Table.read(folder / 'france2018_table2_SpTs.csv')
france3 = table.Table.read(folder / 'france2018_table3.csv')
france4 = table.Table.read(folder / 'france2018_table4.csv')
france_hosts = table.join(france1, france3)
france_nonhosts = table.join(france2, france4)
france = table.vstack((france_hosts, france_nonhosts))
format_numeric_cols(france)

#%% filter to just M stars

france_Ms = ((np.char.count(france['SpT'], 'M') > 0)
             & (np.char.count(france['SpT'], '+') == 0))
france = france[france_Ms]

#%% get simbad id and distance

france_names = france['name'].tolist()
france_info = get_simbad_info(france_names)
france = table.hstack((france, france_info))

#%% convert fluxes to 1 AU and standardize line flux names

line_names = 'C II,Si III,Si IV,N V'.split(',')
line_names = np.char.add(line_names, ' (1e-15)')
col_names = line_names.tolist() + 'Fbol (1e-7),Feuv (1e-14)'.split(',')
dfac = (france['dist']*u.pc/u.au)**2
dfac = dfac.to_value('')
for name in col_names:
    fac, = re.findall(r'\de-\d+', name)
    fac = float(fac)
    france[name] *= fac * dfac
    newname = re.sub(r' \(\de-\d+\)', '', name)
    if 'F' in name:
        newname = newname.replace('F', 'F-')
    else:
        newname = 'F-' + newname.replace(' ', '')
    france.rename_column(name, newname)

#%% filter france columns

keep_france_cols = ['main_id', 'SpT'] + [name for name in france.colnames if 'F-' in name]
france = france[keep_france_cols]


#%% load linsky table

linsky = table.Table.read(folder / 'linsky2020_table1.csv')
format_numeric_cols(linsky)

#%% filter for M stars

linsky_Ms = np.char.count(linsky['SpT'], 'M') > 0
linsky = linsky[linsky_Ms]

#%% get simbad id and distance

linsky_names = linsky['alternate name'].tolist()
for i, name in enumerate(linsky_names):
    if not name or len(re.findall('GJ|HD|HIP', name)) == 0:
        linsky_names[i] = linsky['name'][i]
linsky_simbad_map = {'HD 79210J':'GJ 338A'}
linsky_names = [linsky_simbad_map.get(name, name) for name in linsky_names]
linsky_info = get_simbad_info(linsky_names)
linsky = table.hstack((linsky, linsky_info))

#%% correct for distance errors

correction_factor = (linsky['dist']/linsky['distance'])**2
linsky['Flya'] *= correction_factor


#%% filter linsky columns

linsky.rename_column('Flya', 'F-Lya')
linsky = linsky[['main_id', 'F-Lya', 'SpT']]


#%% join france and linsky

cat = merge_tables(france, linsky, priority=1)


#%% load melbourne table

melbourne = table.Table.read(folder / 'melbourne2020_data.mrt', format='ascii.mrt')
format_numeric_cols(melbourne)

#%% get simbad info

melbourne_names = melbourne['Name'].tolist()
mel_info = get_simbad_info(melbourne_names)
melbourne = table.hstack((melbourne, mel_info))


#%% convert all melbourne luminosities to 1 au

dfac = 1/(4*np.pi*u.au**2)
dfac = dfac.to_value('cm-2')

line_names = 'CII,MgII,HeII,SiII,SiIII,CIV,NV,Lya'.split(',')
col_names = np.char.add('L-', line_names)

for name in col_names:
    fluxes = melbourne[name] * dfac
    newname = name.replace('L-','F-')
    melbourne[newname] = fluxes


#%% filter melbourne columns

melbourne.rename_column('SpType', 'SpT')
keep_mel_cols = ['main_id', 'SpT'] + [name for name in melbourne.colnames if 'F-' in name]
melbourne = melbourne[keep_mel_cols]


#%% join melbourne

cat = merge_tables(cat, melbourne, priority=2)


#%% load pineda table

pineda = table.Table.read(folder / 'pineda2021_table5.csv')
pineda_lines = 'Lyalpha,He II,C II,C III,C IV,Si III,Si IV,N V'.split(',')
format_numeric_cols(pineda)

#%% get simbad info

pineda_names = pineda['Name'].tolist()
pineda_info = get_simbad_info(pineda_names)
pineda = table.hstack((pineda, pineda_info))


#%% convert fluxes to 1 AU, standardize column names, mask limits

pineda_fac = 1e-14
dfac = (pineda['dist']*u.pc/u.au)**2
dfac = dfac.to_value('')
for name in pineda.colnames[1:-2]:
    if 'limit' in name:
        continue
    pineda[name] *= dfac * pineda_fac

#%% mask limits
for name in pineda.colnames[1:-2]:
    if 'limit' in name:
        mask = pineda[name].filled('') != ''
        fluxname = name.replace(' limit', '')
        if not hasattr(pineda[fluxname], 'mask'):
            pineda[fluxname] = table.MaskedColumn(pineda[fluxname], mask=False)
        pineda[fluxname].mask |= mask
        continue

#%% standardize column names
for name in pineda_lines:
    if 'Lya' in name:
        newname = name.replace('Lyalpha', 'F-Lya')
    else:
        newname = 'F-' + name.replace(' ', '')
    pineda.rename_column(name, newname)

#%% filter pineda columns

keep_pineda_cols = ['main_id'] + [name for name in pineda.colnames if 'F-' in name]
pineda = pineda[keep_pineda_cols]

#%% join pineda table

cat = merge_tables(cat, pineda, priority=2)

#%% load loyd table

loyd = table.Table.read(folder / 'loyd2021_line_fluxes.ecsv')
format_numeric_cols(loyd)

#%% get simbad info

loyd_names = loyd['name'].tolist()
loyd_info = get_simbad_info(loyd_names)
loyd = table.hstack((loyd, loyd_info))


#%% convert fluxes to 1 AU, standardize column names, mask limits

dfac = (loyd['dist']*u.pc/u.au)**2
dfac = dfac.to_value('')
for name in loyd.colnames:
    if loyd[name].dtype == float:
        loyd[name] *= dfac

#%% standardize column names
roman_map = {1:'I', 2:'II', 3:'III', 4:'IV', 5:'V'}
for name in loyd.colnames:
    if 'Ftot' in name:
        (element, number, err), = re.findall(r'(\w+)(\d)_Ftot(_?\w*)', name)
        if err == '':
            roman = roman_map[int(number)]
            newname = f'F-{element.title()}{roman}'
            loyd.rename_column(name, newname)

#%% filter loyd columns

keep_loyd_cols = ['main_id'] + [name for name in loyd.colnames if 'F-' in name]
loyd = loyd[keep_loyd_cols]


#%% join loyd table

cat = merge_tables(cat, loyd, priority=2)


#%% record units
for name in cat.colnames:
    if cat[name].dtype == float:
        cat[name].unit = u.Unit('erg s-1 cm-2')


#%% save merged tables

cat.write(folder / 'all_lines_compiled.ecsv', overwrite=True)