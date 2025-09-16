


## file for dealing with stellar properties...

## need different file on the profile shape and probability of overlight


# flow -> input star things, then compute stuff

# minimum should be distance, K band, rotation period


import numpy as np

import astropy.units as u
import astropy.constants as const
from astropy.coordinates import SkyCoord  # for use with ISM stuff



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
    
    def __init__(self,coords, distance,Kband=None,Mass=None,Period=None,RV=None,Gaia=None,Nh=17):
        ""
        # use K-band to get mass
        # rest of properties defined by sequences from Pineda et al. 2021
        
        self.dist = distance.to(u.pc)
        self.coords = SkyCoord(coords[0], coords[1])
        
        if (Mass is None) & (Kband is None)
            print('Must set 2MASS Ks band with distance or else specify mass in solar units')
            return None
        
        if (Mass is None) & (Kband is not None):
            self.MK = absMag(Kband,distance) # k_s from 2MASS
            self.mass = MassMK_Mann19(self.MK)

        if RV is not None:
            self.rv = RV.to('km s-1')

        self.Nh = Nh * u.Unit('dex(cm-2)') # default value of log10 Nh /cm2 of 17

        # if age is known can use loyd relationships, requires radius to convert from surface flux (photometric radii based on Gaia color?)
        # if known early M + period known can also use loyd relationships, requires radius (photometric radii based on Gaia color?)
        # if mid to late dwarf + period known -> uses Pineda, requires 2MASS Kband to get mass (uses mas-lum relationship; typical metallicity)
        # if early M period unknown -- assume active from loyd, (photometric radii based on Gaia color?)
        # if mid to late period unknown -- assume active from Pineda (still need Kband to get mass)
        # if no period use optical relationships

    def populateLines(self):
        ""
        #use to autotmatically update emission line quiescent values based on available inputs
        #have separate routine to take advantage of use inputs?

    def setLines(self,input):
        ""
        self.inlines = input
        
        # use rest of method to convert inputs to standard format, and populate missing lines
        

    def setFFD(self,dFFDin):
        ""
        ## check if dFFDin is normalized, should be not unity

        self.alpha = dFFDin.alpha
        self.delta_min = dFFDin.d_m # has units
        self.flarerate = dFFDin.norm # has units


    def setContrast(self,contrast_fit):
        ""
        # maybe need to have this for a list of lines or pull from catalog of lines based on observatory input
        self.ampCurve = contrast_fit




## can then use relationships from pineda et al 2021b? R, Lum --> Teff ?? don't really need these may need luminosity actually for normalized

