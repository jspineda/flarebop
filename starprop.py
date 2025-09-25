


## file for dealing with stellar properties...

## need different file on the profile shape and probability of overlight


# flow -> input star things, then compute stuff

# minimum should be distance, K band, rotation period


import numpy as np
from scipy.interopolate import make_interp_spline

import astropy.units as u
import astropy.constants as const
from astropy.coordinates import SkyCoord  # for use with ISM stuff

from astroquery.simbad import Simbad
from astroquery.gaia import Gaia

import pdb

Simbad.reset_votable_fields()
Simbad.add_votable_fields('parallax')
Simbad.add_votable_fields('typed_id')
Simbad.add_votable_fields('id')
Simbad.add_votable_fields('rv_value')
Simbad.add_votable_fields('sptype')
Simbad.add_votable_fields('flux(K)')
#Simbad.add_votable_fields('flux(i)')
#Simbad.add_votable_fields('flux(g)')
#Simbad.add_votable_fields('flux(G)')


### utility functions first

def absMag(apparent,distance):
    "Compute absolute magnitude, assumes distance in PC"
    try:
        distance.cgs
    except:
        distance = distance * u.pc
        
    return apparent - 5*np.log10(  (distance/ (10*u.pc) ).cgs.value)


## From Mann et al. 2019 -- use their code for full consideration of posterior
def MassMK_Mann19(MK):
    "Converts Absolute Ks-band magnitude to stellar mass, at average metallicity"
    coefs = [-2.13e-4,1.42e-4,7.87e-3,-8.43e-4,-0.208,-0.642]
    mass = np.power(10,np.poly1d(coefs)(MK-7.5))
    err = mass * 0.02

    return mass,err

## as a class prevent frequent recalculation of splines
class GPSeq():
    ""
    def __init__(self):
        ## applicable for general M-dwarfs
        ## For returning bolometric luminosity from input mass, based on GP relation of Pineda~2021b
        datLM = Table.read("./tables/LMGPSequence.csv")
        datTM = Table.read("./tables/TMGPSequence.csv")

        self.lsun = const.L_sun.cgs
        self.funcLM = make_interp_spline(datLM['Mass'],datLM['LogL'])
        self.funcTM = make_interp_spline(datTM['Mass'],datTM['Tefg'])

    def getLum(self,mass):
        ""
        ## input mass in solar masses
        return np.power(10, self.funcLM(mass) )*self.lsun

    def getTeff(self,teff):
        ""
        return self.funcTM(teff)*u.K

GPs = GpSeq()

def photRadius(Gin,BpRp):
    ""
    # fill in with appropriate function...
    return 0.4*const.R_sun


def convTime(massin,method="Wright2018"):
    "Compute convective turn over time from input mass, for Rossby calculation; defualt relation from Wright et al. 2018"
    #input mass in solar masses
    # maybe add percriptions from Gossage?
    try:
        mass = (massin.to(u.g) / const.M_sun).to('')
    except:
        raise ValueError('Input Mass needs to have compatible mass units')
    
    if method == "Wright2011":
        coef = [-0.54,-1.49,1.16]
        return np.power(10,np.poly1d(coef)(np.log10(mass))) * u.d # output in days
    elif method == "Wright2018":
        coef = [0.31,-1.50,2.33]
        return np.power(10,np.poly1d(coef)(mass)) * u.d # output in days
    elif method == "Nunez2015":
        coef = [0.56,-1.63,1.24]
        return np.power(10,np.poly1d(coef)(np.log10(mass))) * u.d # output in days
    else:
        print("Check method keyword...")
        return None


def Rossby(mass,period,method="Wright2018"):
    "Compute Rossby Number, period input in days, see convTime for different relations" # Rossby = period / tau
    tau = convTime(mass,method=method)
    return (period / tau).to('')


## stellar properties class for user usage
class StarProp (object):
    
    def __init__(self, infoin , Distance=None, Kband=None, Mass=None, Period=None, RV=None, Gmag=None, Age=None, Nh=17, Verbose=True):
        ""
        #Infoin is either a string for simbad query name, or a tuple of coordinates
        #input keywords override query results
        
        self.input = infoin
        self.verb = Verbose
        # these cannot come from sql queries (at the moment)
        if Age is None:
            self.age = None
        else:
            self.age = Age.to(u.yr) # units?

        if Period is None:
            self.period = None
        else:
            self.period = Period.to(u.d) # units?
        
        self.Nh = Nh * u.Unit('dex(cm-2)') # default value of log10 Nh /cm2 of 17


        if isinstance(infoin,str):
            
            simtab = Simbad.query_objects([infoin])

            if simtab['MAIN_ID'] == '':
                print("Error: SIMBAD query failed")
                print("Use manual inputs..??")
                return None

            self.simbad = simtab
            self.spt = simtab['SP_TYPE'][0]
            ra = simtab['RA'][0].split(' ')
            dec = simtab['DEC'][0].split(' ')
            self.coords = SkyCoord("{0}h{1}m{2}s {3:+}d{4}m{5}s".format(ra[0],ra[1],ra[2],int(dec[0]),dec[1],dec[2]) )

            # what if values are missing?
            rv_sim = simtab['RV_VALUE'][0] * simtab['RV_VALUE'].unit # store a radial velocity
            self.dist = 1000*u.pc/simtab['PLX_VALUE'][0]           # store distance in pc
            self.MK = absMag(simtab['FLUX_K'][0],self.dist)
            self.mass, _ = MassMK_Mann19(self.MK)
            self._mass = self.mass * const.M_sun
            #self.gi = simtab['flux(g)'] - simtab['flux(i)']
            
            self.GaiaID = [i for  i in  simtab['ID'][0].split(',') if 'Gaia DR3' in i][0].strip(' ')
            
            if len(self.GaiaID) != 0:
                adql_query = """
                    SELECT source_id, ra, dec, designation,phot_g_mean_mag,phot_bp_mean_mag,phot_rp_mean_mag,bp_rp
                    FROM gaiadr3.gaia_source
                    WHERE designation = '{}' """.format(self.GaiaID)
                
                job = Gaia.launch_job(adql_query)
                gtab = job.get_results()
                
                if len(gtab) != 1:
                    print("Gaia Query failed..")
                self.G = gtab['phot_g_mean_mag']
                self.Bp_Rp = gtab['bp_rp']
            else:
                print("Gaia ID not found in SIMBAD..")
            
            self.photrad = photRadius(self.G,self.self.Bp_Rp) # need to define phot radius relation for use with loyd rot-act stuff
            self._photrad = self.photrad * const.R_sun
            self.LumBol = GPs(self.mass)
            
        elif isinstance(infoin,tuple):
            
            self.coords = SkyCoord(infoin[0], infoin[1])

        else:
            print("Check Input; needs to be star name (SIMBAD) or coordinates")
            return None
        
        ## over ride query results if keywords are given
        if Distance is not None:
            self.dist = Distance.to(u.pc)
        
        if (Mass is None) & (Kband is not None):
            self.MK = absMag(Kband,self.dist) # k_s from 2MASS
            self.mass = MassMK_Mann19(self.MK)
        
        if Gmag is not None:
            self.G = Gmag[0]
            self.Bp_Rp = Gmag[1]

        if Mass is not None:
            self.mass = Mass # units?

        if RV is not None:
            self.rv = RV.to('km s-1')

        self.state()
            

        # if age is known can use loyd relationships, requires radius to convert from surface flux (photometric radii based on Gaia color?)
        # if known early M + period known can also use loyd relationships, requires radius (photometric radii based on Gaia color?)
        # if mid to late dwarf + period known -> uses Pineda, requires 2MASS Kband to get mass (uses mas-lum relationship; typical metallicity)
        # if early M period unknown -- assume active from loyd, (photometric radii based on Gaia color?)
        # if mid to late period unknown -- assume active from Pineda (still need Kband to get mass)
        # if no period use optical relationships

        # coords is a tuple of RA and DEC
        # use K-band to get mass
        # rest of properties defined by sequences from Pineda et al. 2021

        # first thing could be a simbad query of the coordinates or name to cross check a spectral type (M3 is early otherwise late)
        

    def populateLines(self):
        ""
        #use to autotmatically update emission line quiescent values based on available inputs
        #have separate routine to take advantage of use inputs?

    def setLines(self,input):
        ""
        self.inlines = input
        
        # use rest of method to convert inputs to standard format, and populate missing lines
        

    def setFFDuv(self,dFFDin):
        ""
        ## check if dFFDin is normalized, should be not unity

        self.alpha = dFFDin.alpha
        self.delta_min = dFFDin.d_m # has units
        self.flarerate = dFFDin.norm # has units

        ## separate routine to take stellar properties and turn into optical FFD?

    def setContrast(self,contrast_fit):
        ""
        # maybe need to have this for a list of lines or pull from catalog of lines based on observatory input
        self.ampCurve = contrast_fit

    
    
    def setPeriod(self, perin):
        ""
        self.period = perin.to(u.d)

    def setAge(self, agein):
        ""
        self.age = agein.to(u.yr)

    def state(self):
        ""

        print("Input set as `{}' ".format(self.input))
        print("Target Coordinates are {}".format(self.coords) )
        
        print("Object Mass is set to {} solar masses".format(self.mass))
        
        if self.age is not None:
            print("Object Age is set to {}".format(self.age))

        if self.period is not None:
            print("Object Period is set to {}".format(self.period))

        if (self.period is None) & (self.age is None):
            print("Will need to set Age or Period if quiescent fluxes unknown, use ____ to set quiescent lines")

## if age is unknown what to do about optical FFDs

## Davenport+ 19 has model of optical FFDs as function of age and mass....
## look at his conversion from g-i to mass, and from period to age (mamajek and hillenbrand.. )

class RotAct(object):
    
    def __init__(self,starclass):
        ""
        
        if starclass.mass is not None and starclass.period is not None:
            self.rossby = Rossby(starclass.mass,starclass.period)
        else:
            self.rossby = None

        # takes as input starclass
        # examines stellar properties in class, to choose which of the rotation activity relations to utilize?
        # returns method defined by those functions and annotations of that choice
        
        acttab = Table.read('tables/rot_act.dat')
        
        # place holder
        _select = (
    (acttab['xparam'] == 'Rossby') &
    (acttab['xband'] == 'Pineda') &
    (acttab['yband'] == 'HeII')  )

        prop = acttab[_select]


#        _select = (
#    (acttab['xparam'] == 'Rossby') &
#    (acttab['xband'] == 'Pineda') &
#    (flare_correlations['yparam'] == 'contrast') &
#    (flare_correlations['yband'] != 'fuv130') &
#    (flare_correlations['percentile'] == '95')
#)

    def func(self,varin):
    
        return 0
    
    
## can then use relationships from pineda et al 2021b? R, Lum --> Teff ?? don't really need these may need luminosity actually for normalized

