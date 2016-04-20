from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt

modelfile = '../../RG_ELCmodeling/9291629/trial5/modelU.mag'
datafile =  '../../RG_ELCmodeling/9291629/trial5/KIC_9291629_LC_mag_Q017.txt'

KIC = '9291629'; period = 20.68639; BJD0 = 2454966.882
magdim = 14.17; magbright = 13.93
timemin = 100; timemax = 1600

mtimes, mmags, merrs = np.loadtxt(modelfile, usecols=(0,1,2), unpack=True)
dtimes, dmags, derrs = np.loadtxt(datafile, usecols=(0,1,2), unpack=True)

plt.plot(mtimes, mmags, color='k', ls='None', marker='.')
plt.plot(dtimes, dmags, color='r', ls='None', marker='.')
plt.axis([timemin, timemax, magdim, magbright])
plt.show()