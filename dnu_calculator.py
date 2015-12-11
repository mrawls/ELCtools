from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
'''
Meredith Rawls, 2015
Takes masses, radii, and chi2s from an ELC run and makes delta_nu distributions.
Assumes a fixed temperature for each star.
'''
dir = '../../RG_ELCmodeling/9246715/demcmc001/'
#dir = '../../RG_ELCmodeling/7037405/trial3/'
starparmfile = 'starparm.all'
chifile = 'chi.all'
outfile = 'deltanu_ELCcalc.txt'
nplotbins = 100

print('Reading in giant files, be patient...')
M1s, M2s, R1s, R2s = np.loadtxt(dir+starparmfile, usecols=(0,1,2,3), unpack=True)
chi2s = np.genfromtxt(dir+chifile, usecols=(1,), unpack=True)
print('Finished reading in giant files!')
T1 = 5000.
T2 = 5000.

def calc_dnu(mass, radius, temp):
    density_sun = (1. * u.Msun) /(4./3. * np.pi * np.power((1. * u.Rsun),3))
    dnu_sun = 135.5 # muHz
    density = (mass * u.Msun) /(4./3. * np.pi * np.power((radius * u.Rsun),3))
    dnu = dnu_sun * np.sqrt(density/density_sun)
    return dnu

print('Calculating delta nu for each model, be patient...')
dnu1s = []
dnu2s = []
for M1, M2, R1, R2 in zip(M1s, M2s, R1s, R2s):
    dnu1s.append(calc_dnu(M1, R1, T1))
    dnu2s.append(calc_dnu(M2, R2, T2))
print('Finished calculating delta nus!')

print(np.median(dnu1s), np.std(dnu2s))
print(np.median(dnu2s), np.std(dnu2s))
delta_dnus = np.array(dnu1s) - np.array(dnu2s)
print(np.median(delta_dnus), np.std(delta_dnus))

threshold = 0.5 # muHz, typical mode width
count = 0 # tally for how many delta_dnus are within some threshold
f1 = open(dir+outfile, 'w')
for chi2, dnu1, dnu2, delta_dnu in zip(chi2s, dnu1s, dnu2s, delta_dnus):
    print(chi2, dnu1, dnu2, delta_dnu, file=f1)
    if delta_dnu < threshold:
        count += 1
f1.close()

print('Out of {0} models, {1} ({2}\%) have a difference in dnu less than {3}.'.format \
    (len(chi2s), count, float(count)/float(len(chi2s)), threshold))

# Plot histograms of dnu1s, dnu2s, and the difference
fig = plt.figure()
ax1 = fig.add_subplot(1, 3, 1)
x1 = plt.xlabel('Delta nu 1')
histogram = plt.hist(dnu1s, nplotbins, histtype='stepfilled')
ax2 = fig.add_subplot(1, 3, 2)
x2 = plt.xlabel('Delta nu 2')
histogram = plt.hist(dnu2s, nplotbins, histtype='stepfilled')
ax3 = fig.add_subplot(1, 3, 3)
x3 = plt.xlabel('Delta nu 1 - Delta nu 2')
histogram = plt.hist(delta_dnus, nplotbins, histtype='stepfilled')

plt.show()