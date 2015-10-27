from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.stats import norm
'''
Graph some of the fit parameters from an ELC run, histogram style

INPUT:
Correctly specify the directory where the files listed below are listed
(note that the fitparm.all, starparm.all, chi.all files must be created by
'cat'-ing demcmc_fitparm.1*, demcmc_starparm.1*, chi.1* respectively)
--> For decent CDF plots, set N ~ number of ELC models in the cdferrplot function

OUTPUT:
Prints each fit parameter and derived parameter to the screen
Makes a set of histograms for fit parameters and derived parameters
Option to make individual CDF plots when calling the cdferrplot function
'''

dirstub = '../../RG_ELCmodeling/9246715/demcmc001/'
#dirstub = 'demcmc_chunk1/'

fitparmfile =   dirstub + 'fitparm.all'
starparmfile =  dirstub + 'starparm.all'
gridloopfile =  dirstub + 'gridloop.opt'
parmkeyfile =   dirstub + 'key.ELCparm'
chi2file =      dirstub + 'chi.all'

burnin = 10000
nplotbins = 100

def cdferrplot(var, varname, N=100000, plot=True):
    '''
    Awesome function originally written by Jean!
    Plots a cumulative distribution function, and calculates value +/- 1-sigma errors
    a = 50%, b = 15.75%, c = 84.25%
    Thus, a = value, c-a = upper 1-sigma error, a-b = lower 1-sigma error
    '''
    cdf = plt.hist(var, bins=N, normed=True, cumulative=True, histtype='step', color='k')
    ai = np.where(cdf[0] < 0.5)[0][-1]
    bi = np.where(cdf[0] < 0.1575)[0][-1]
    ci = np.where(cdf[0] < 0.8425)[0][-1]
    a = (cdf[1][ai] + cdf[1][ai+1])/2
    b = (cdf[1][bi] + cdf[1][bi+1])/2
    c = (cdf[1][ci] + cdf[1][ci+1])/2
    if plot == True:
        plt.figure()
        ##fig = plt.figure(3, figsize=(15,10))
        ##for idx in range(1,16):
        ##ax = fig.add_subplot(4, 4, idx)
        ##hist = plt.hist(var, bins=N, normed=True, cumulative=True, histtype='step', color='k')
        plt.ylabel(varname)
        plt.axvline(a, linewidth=2, color='k')
        plt.axvline(b, linestyle='dashed', color='k')
        plt.axvline(c, linestyle='dashed', color='k')
        plt.ylim(0, 1)
        plt.xlim(min(var), max(var))
        plt.show()
    print('{0} = {1} +{2} -{3}'.format(varname, a, c-a, a-b))
    return

# Read in names for FITPARM variables from gridloop file
gridloop = [line.rstrip('\n') for line in open(gridloopfile)]
nfitparms = int(gridloop[10]) # reads the number of fit variables from gridloop file
fitparmnames = []
for i in range(0, nfitparms):
    fitparmnames.append(gridloop[10+i+1].rstrip())
# manually include 2 systemic velocity columns and 4 final columns (t0, tconj, ecc, argper)
# this increases the number of columns in fitparmnames by 6
nfitparms += 6
fitparmnames.append('gamma1'); fitparmnames.append('gamma2'); fitparmnames.append('t0v2')
fitparmnames.append('tconjv2'); fitparmnames.append('ecc'); fitparmnames.append('argper')

# Read in names for STARPARM variables from parmkey file
starparmnames = np.loadtxt(parmkeyfile, comments='#', usecols=(1,), dtype={'names':('starparmnames',),'formats':('|S11',)}, unpack=True)
for idx, entry in enumerate(starparmnames): # remove 'index' from parmkeyfile
    entry = str(entry)
    if ('index' in entry):
        starparmnames = np.delete(starparmnames, idx, axis=0) # remove 'chi^2' from parmkey file
for idx, entry in enumerate(starparmnames):
    entry = str(entry)
    if ('chi^2' in entry):
        starparmnames = np.delete(starparmnames, idx, axis=0)
nstarparms = len(starparmnames)

print('Reading in fitparm.all, starparm.all, and chi.all, please be patient...')

fitparms = np.loadtxt(fitparmfile, usecols=(range(0,nfitparms)), dtype=np.float64, unpack=True)
starparms = np.loadtxt(starparmfile, usecols=(range(0,nstarparms)), dtype=np.float64, unpack=True)
chi2s = np.loadtxt(chi2file, usecols=(1,), dtype=np.float64, unpack=True)

# Ensure the FITPARM, STARPARM, and CHI2 arrays are all the same length
# Omit the first trials corresponding to the burn-in period
newlength = min(fitparms.shape[1], starparms.shape[1], chi2s.shape[0]) # find length of shortest array
fitparms = fitparms[:,burnin:newlength]
starparms = starparms[:,burnin:newlength]
chi2s = chi2s[burnin:newlength]

print('Finished reading everything in, calculating values and generating plots...')

# First plot: histograms of all the FIT PARAMETERS (FITPARM)
fig = plt.figure(1, figsize=(15,10))
windowcols = 4
windowrows = int([np.rint(nfitparms/windowcols) if (np.float(nfitparms)/windowcols)%windowcols == 0 else np.rint(nfitparms/windowcols)+1][0])
for idx, param in enumerate(fitparms):
    paramname = fitparmnames[idx]
    ax = fig.add_subplot(windowrows, windowcols, idx+1)
    small = plt.tick_params(axis='both', which='major', labelsize=10)
    for label in ax.get_xticklabels()[::2]: # hide every other tick label
        label.set_visible(False)
    xval = param
    xmin = np.min(param)
    xmax = np.max(param)
    histogram = plt.hist(xval, nplotbins, histtype='stepfilled')
    ymin, ymax = ax.get_ylim()
    y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    yformat = ax.yaxis.set_major_formatter(y_formatter)
    x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    xformat = ax.xaxis.set_major_formatter(x_formatter)
    label = plt.text(xmin + 0.1*(np.abs(xmax-xmin)), 0.8*ymax, paramname, size=20)
    cdferrplot(param, paramname, newlength, plot=False)
    # set plot=True if you want a CDF plot in a new window for each parameter

# Second plot: histograms of all the DERIVED PARAMETERS (STARPARM)
fig = plt.figure(2, figsize=(15,10))
windowcols = 4
windowrows = 4
#windowrows = int([np.rint(nfitparms/windowcols) if (np.float(nfitparms)/windowcols)%windowcols == 0 else np.rint(nfitparms/windowcols)+1][0])
for idx, param in enumerate(starparms[0:16]):   # we only care about the first 16 starparms
    paramname = str(starparmnames[idx])[2:-3]
    ax = fig.add_subplot(windowrows, windowcols, idx+1)
    nope = ax.set_yticklabels(())
    small = plt.tick_params(axis='both', which='major', labelsize=10)
    for label in ax.get_xticklabels()[::2]: # hide every other tick label
        label.set_visible(False)
    xval = param
    xmin = np.min(param)
    xmax = np.max(param)
    histogram = plt.hist(xval, nplotbins, histtype='stepfilled')
    ymin, ymax = ax.get_ylim()
    y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    yformat = ax.yaxis.set_major_formatter(y_formatter)
    x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    xformat = ax.xaxis.set_major_formatter(x_formatter)
    label = plt.text(xmin + 0.1*(np.abs(xmax-xmin)), 0.8*ymax, paramname, size=20)
    cdferrplot(param, paramname, newlength, plot=True)
    # set plot=True if you want a CDF plot in a new window for each parameter
    
plt.show()