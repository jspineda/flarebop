#

import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt


from scipy.special import erf
from scipy.optimize import root




class FlProb(object):
    ""
    def __init__(self,starclass,obsclass,Dt=1*u.hr,refprob=1e-4):
        ""
        
        # starclass is meant to be 'StarProp 'class in starprop.py; and obsclass is whatever we do for the observatory stuff 
        
        self.Dt = Dt
        self.ampH = 13906.603563778966  # 3581.92876456 is the standard reference for fuv130 FFD, and flare contrast 50%
        self.refprob = refprob # store reference probability for comparison

        ## replace with actual compute H; also check consistency of contrast curve choice with observatory selections
    
        #self.star = starclass
        
        self.alpha = starclass.alpha # slope in differential ffd
        self.d_m = starclass.delta_min.to(u.s) # minimun in units of seconds
        self.rate = starclass.flarerate # total flare rate with units
        
        self.sig = starclass.ampCurve['scatter'][0] # this will need to match however flare contrast curves as stored in starprop
        self.aslpe = starclass.ampCurve['alpha'][0]
        self.boffset = starclass.ampCurve['beta'][0] - self.aslpe * starclass.ampCurve['pivot'][0] # data tables use pivot as part of offset

        self.dcrit = self.matchED(self.ampH)
        self.dpeak = self.peakED(self.ampH)
        self.probability = self.totprob(self.ampH,self.Dt)
        self.Hlimit = self.maxThresh(self.refprob,self.Dt)
        self.percentilDP = self.percAmp(self.ampH,self.dpeak)

    def logAmplitude(self, delta):
        "Returns Log10 Flare Constrast form input Equivalent Duration (s)"
        return self.aslpe*np.log10(delta.to('s').value) + self.boffset

    def _logAmplitude(self, delta):
        ""
        # dimensionless utility
        return self.aslpe*np.log10(delta)  + self.boffset

    def amplitude(self, delta):
        "Returns Flare Constrast from input Equivalent Duration (s)"
        return np.power(10,self.logAmplitude(delta.to(u.s).value))


    def eta(self, delta, H, time):
        ""
        prefac = 0.5*(self.alpha - 1)*(1 - np.exp( - (self.rate*time).to('') )) * np.power( (delta / self.d_m).to('') , -self.alpha ) / self.d_m
        f1 = 1 - erf( (np.log10(H) - self.logAmplitude(delta) ) / np.sqrt(2) / self.sig  )
        return f1*prefac
        

    def peakED(self,H_in):
        "Compute Equivalent Duration at Peak Probability Density"
        
        maxdiv = lambda x: self.aslpe*np.exp( -( (np.log10(H_in) -  self._logAmplitude( x ) )/(self.sig*np.sqrt(2))  )**2
                                            )*(self.sig*np.log(10)*np.sqrt(2*np.pi) )**-1 + self.alpha*0.5*erf( (np.log10(H_in) -  self._logAmplitude(x) )/(self.sig*np.sqrt(2) ) ) - 0.5*self.alpha

        outm = root(maxdiv, self.dcrit.to(u.s).value )

        return outm.x * u.s
    

    def matchED(self,H_in):
        ""
        ## return equiv duration matching input amplitude
        return np.power(10, (np.log10(H_in) - self.boffset) / (self.aslpe) )*u.s


    def totprob(self,H_in,Dt):
        ""
        ## do calculation for total probability
    
        zm = ( np.log10(H_in) - self.aslpe*np.log10(self.d_m.to(u.s).value) - self.boffset) / self.sig / np.sqrt(2)
        xi = (self.alpha-1) * self.sig * np.sqrt(2) * np.log(10) / self.aslpe

        rscale = 0.5*(1 - np.exp( -(self.rate*Dt).to('') ) )
        
        f_H =  np.power(self.d_m.to(u.s).value,self.alpha-1 ) * np.exp(
                                (self.alpha-1)*np.log(10)*self.boffset/self.aslpe) * np.power(H_in, (1- self.alpha)/ self.aslpe)
    
        integ = erf(zm)*np.exp(xi*zm) - np.exp(xi**2 / 4)*(1 + erf(zm - xi*0.5) )
    
        return rscale*( 1 - f_H*integ)


    def maxThresh(self,prob, Dt ):
        ""
        ## inverts probability equation using input time to threshold associated with probability, i.e. H needs to be big enough so all events are less probable than input reference.
    
        newf = lambda x,dt: self.totprob(x,dt) - prob

        outH = root(newf,1000,args=(Dt)) # shouldn't be too much issue with input guess, function is pretty monotonic

        return outH.x[0]
        
        
    def percAmp(self,H_in,delta):
        ""
        ## use for evaluating the percentiles at H_in given delta; applies to the input contraast curve amplitude
        return 0.5*(1 + erf( (np.log10(H_in) -  self.logAmplitude(delta) )/(self.sig*np.sqrt(2))  )  )


    def runDiag(self):
        ""
        ## routine to generate diagnostic plots and report? current stuff here is a bit of place holder / example

        print("Overlight threshold is {} in flare contrast".format(self.ampH))
        print("Total Overlight Probability is {0} in {1}".format(self.probability, self.Dt) )
        print("Danger peaks for Equivalent Durations of {}".format(self.dpeak) )


#        y = np.logspace( np.log10(self.dpeak.value) - 3, np.log10(self.dpeak.value) + 3, 1000 )*u.s
#        eta = self.eta(y)
#        plt.figure(1)
#        plt.plot(y,eta,label="H = {}".format(self.ampH) )
#        plt.vlines(self.dpeak,eta.max().value*1e-3,eta.max().value*2,color='gray',linestyles=':',alpha=0.8)
#        plt.legend()
#        plt.xscale('log')
#        plt.yscale('log')
#        plt.ylabel('$\\eta$ (s$^{-1}$) ')
#        plt.xlabel('Equivalent Duraation (s)')
#        plt.ylim([eta.max().value*1e-3,eta.max().value*2])
#

#        Hs = np.logspace(np.log10(2),4,1000)
#        zed = self.totprob(Hs,self.Dt)
#        plt.figure(2)
#        plt.plot(Hs, zed ,label="$\\Delta t$ = {}".format(self.Dt) )
#        plt.hlines(self.refprob,Hs.min(),Hs.max(),color='gray',linestyles=':',alpha=0.8,label='Policy Reference')
#        plt.vlines(self.Hlimit, zed.min() ,1,color='gray',linestyles=':',alpha=0.8)
#        plt.text(self.Hlimit*1.5,5e-3,"Limiting H = {0:5.1f}".format(self.Hlimit))
#
#        plt.legend()
#        #plt.xscale('log')
#        plt.yscale('log')
#        plt.ylabel('Event Probability')
#        plt.xlabel('Threshold Amplitude, H')
##
