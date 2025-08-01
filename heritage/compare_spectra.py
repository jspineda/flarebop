from astropy import units as u
import numpy as np
from matplotlib import pyplot as plt

import lines
from heritage import heritage_isr_tools as isr2017

wavegrid = np.arange(1100, 1800, 0.05) * u.AA

#%%
flxunt = u.Unit('erg s-1 cm-2')
examples = {
    'AU Mic active early-M': dict(
        known_fluxes=dict(CIII=(2.34+1.57+0.96+3.98+2.05+1.80)*1e-14 * flxunt, # pagano+ 2000
                          SiIII=7.90e-14 * flxunt,
                          NV=(4.52+2.11)*1e-14*flxunt,
                          OI=(1.92+4.89+4.99)*1e-14*flxunt,
                          CII=(6.89+14.39)*1e-14*flxunt,
                          SiIV=(4.93+2.92)*1e-14 * flxunt,
                          CIV=(21.36+11.34)*1e-14*flxunt,
                          HeII=26.35*1e-14*flxunt),
        ra=311.2897187392163 * u.deg,
        dec=-31.3408994138447 * u.deg,
        distance=9.722100 * u.pc,
        rv_star=-6.31 * u.km / u.s,
        rv_ism='local ism',
        Nh=17 * u.Unit('dex(cm-2)'),
        category='early M',
        activity_class='active',
        isr2017_fluxes=dict(Uflare=8.44,
                            Fc4=3.22E-11 * flxunt,
                            Fsi4=2.97E-11 * flxunt,
                            Flya=1.19E-09 * flxunt)
    ),

    'AD Leo active mid-M' : dict(
        known_fluxes=dict(CII=220e-15*flxunt, # france+ 2018
                          SiIII=140e-15*flxunt,
                          SiIV=160e-15*flxunt),
        ra=154.90117008*u.deg,
        dec=19.8700029*u.deg,
        distance=4.97*u.pc,
        rv_star=12.43*u.km/u.s,
        rv_ism='local ism',
        Nh=17*u.Unit('dex(cm-2)'),
        category='mid M',
        activity_class='active',
        isr2017_fluxes=dict(Uflare=4.042,
                            Fc4=1.34E-09*flxunt,
                            Fsi4=2.56E-09*flxunt,
                            Flya=4.95E-08*flxunt)
    ),

    'TRAPPIST-1 active late-M': dict(
        known_fluxes=dict(CIII=101e-18 * flxunt, # wilson+ 2021
                          SiIII=25e-18 * flxunt,
                          Lya=1.4e-14 * flxunt,
                          NV=(36+16)*1e-18*flxunt,
                          CII=(30+55)*1e-18*flxunt,
                          SiIV=(25+14)*1e-18 * flxunt,
                          CIV=(100+53)*1e-18*flxunt),
        ra=346.6223687285788 * u.deg,
        dec=-05.0413992505183 * u.deg,
        distance=12.43 * u.pc,
        rv_star=-52 * u.km / u.s,
        rv_ism='local ism',
        Nh=17 * u.Unit('dex(cm-2)'),
        category='late M',
        activity_class='active',
        isr2017_fluxes=dict(Uflare=19.63,
                            Fc4=2.46E-15 * flxunt,
                            Fsi4=3.54E-16 * flxunt,
                            Flya=9.09E-14 * flxunt)
    ),

    'GJ 667 C inactive early-M': dict(
        known_fluxes=dict(SiIII=0.51e-15 * flxunt,
                          NV=0.72e-15 * flxunt,
                          CII=0.65e-15 * flxunt,
                          SiIV=0.83e-15 * flxunt),
        ra=259.7451137490841 * u.deg,
        dec=-34.9968368537172 * u.deg,
        distance=6.97107 * u.pc,
        rv_star=6.368 * u.km / u.s,
        rv_ism='local ism',
        Nh=17 * u.Unit('dex(cm-2)'),
        category='early M',
        activity_class='inactive',
        isr2017_fluxes=dict(Uflare=10.66,
                            Fc4=4.91E-12 * flxunt,
                            Fsi4=3.13E-12 * flxunt,
                            Flya=1.82E-10 * flxunt)
    ),

    'GJ 436 inactive mid-M': dict(
        known_fluxes=dict(SiIII=0.52e-15 * flxunt,
                          NV=0.96e-15 * flxunt,
                          CII=1.09e-15 * flxunt,
                          SiIV=0.68e-15 * flxunt),
        ra=175.5462222295741 * u.deg,
        dec=26.7065696618828 * u.deg,
        distance=9.75321 * u.pc,
        rv_star=9.59 * u.km / u.s,
        rv_ism='local ism',
        Nh=17 * u.Unit('dex(cm-2)'),
        category='mid M',
        activity_class='inactive',
        isr2017_fluxes=dict(Uflare=10.941,
                            Fc4=3.87E-12 * flxunt,
                            Fsi4=2.36E-12 * flxunt,
                            Flya=1.43E-10 * flxunt)
    )
    # inactive late-Ms might not exist. Even Prox Cen with its 80+ day rotation period has
    # an Ha equivalent width > 1 Å in emission
}

#%%
for key, example in examples.items():
    oldfluxes = example['isr2017_fluxes']
    oldspec = isr2017.heritage_isr_spec(
        wavegrid,
        oldfluxes['Fc4'],
        oldfluxes['Fsi4'],
        oldfluxes['Flya'],
        oldfluxes['Uflare']
    )
    newspec = lines.automated_line_spectrum(
        wavegrid,
        example['known_fluxes'],
        example['ra'],
        example['dec'],
        example['distance'],
        example['rv_star'],
        example['rv_ism'],
        example['Nh']
    )
    newspec = lines.automated_line_spectrum(wavegrid, example['known_fluxes'], example['ra'], example['dec'], example['distance'], example['rv_star'], example['rv_ism'], example['Nh'])
    
    oldspec = oldspec.to_value('erg s-1 cm-2 AA-1')
    newspec = newspec.to_value('erg s-1 cm-2 AA-1')

    plt.figure()
    plt.plot(wavegrid, oldspec, lw=3, label='2017 Procedure')
    plt.plot(wavegrid, newspec, label='New Procedure')
    plt.xlabel('Wavelength (Å)')
    plt.ylabel('Flux Density (erg s-1 cm-2 Å-1)')
    plt.ylim(np.min(oldspec)/2, 2*np.max([oldspec, newspec]))
    plt.yscale('log')
    plt.title(key)
    plt.tight_layout()

    path = f'heritage/comparison_plots/{key}.pdf'
    plt.savefig(path)