
## Run for generating data underlying function
## mapping threshold to flare duration above that threshold
## and other aspects for flare profile usage

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import root

from flprofile import flare_model, flare_peak



## Example of profile shape based on Mendoza et al. 2022
t = np.logspace(0,1.7,1000) - 3

fig, ax = plt.subplots()
ax.plot(t, flare_model(t,0,1,1) )

ax.hlines(flare_peak,t.min(),t.max(),linestyles='--',color='gray',alpha=0.5)

ax.text(15,0.7,"Flare Peak at {0:5.4f}".format(flare_peak))

ax.set_ylabel("Relative Amplitude")
ax.set_xlabel("Relative Time [Units of Flare FWHM]")

fig.savefig("./fits/plots/profile diagnostics/flareprofile.pdf")



### constraining, and generating flare duration prescriptions
## Dur = amount of time flare spends above threshold
N = 200 #baseline sampling
xin = np.linspace(0.01,flare_peak*0.9975,N)

relflare = lambda x: flare_model(x,0,1,1) - xin

out_l = root(relflare, -0.3*np.ones_like(xin))
out_u = root(relflare, 2*np.ones_like(xin))

#out_u.x

dur = out_u.x - out_l.x
dur = np.append(dur,0)
xin = np.append(xin,flare_peak)

fig1, ax1 = plt.subplots()
ax1.plot(xin,dur)

ax1.set_ylabel("Threshold Duration [Relative Time]")
ax1.set_xlabel("Relative Amplitude Threshold")

fig1.savefig("./fits/plots/profile diagnostics/flaredurations.pdf")


nfile = open("./fits/output tables/flare_prof_dur.dat",'w')
nfile.write("# Flare profile durations spend above given threshold \n")
nfile.write("# threshold,  duration [rel. time] \n")

for j in range(len(dur)):
    nfile.write("{0}, {1} \n".format(xin[j],dur[j]))

nfile.close()
