from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import IndexLocator, FormatStrFormatter
import matplotlib.mlab as mlab
from scipy.stats import norm
import gaussfitter as gf
'''
Graph some of the fit parameters from an ELC run, histogram style
This program is very rough and needs to be improved, but it's a start
'''
fitparmfile =   '../../RG_ELCmodeling/9246715/demcmc001/fitparm.all'
starparmfile =  '../../RG_ELCmodeling/9246715/demcmc001/starparm.all'

def cdferrplot(var,varname):
    '''
    Amazing function written by Jean
    Something something cumulative distribution gaussian errorbar magic something
    '''
    plt.figure()
    cdf = plt.hist(var,bins=N,normed=True,cumulative=True,histtype='step',color='k')
    ai = np.where(cdf[0]<.5)[0][-1]
    bi = np.where(cdf[0]<.1575)[0][-1]
    ci = np.where(cdf[0]<.8425)[0][-1]
    a = (cdf[1][ai]+cdf[1][ai+1])/2
    b = (cdf[1][bi]+cdf[1][bi+1])/2
    c = (cdf[1][ci]+cdf[1][ci+1])/2
    plt.ylabel(varname)
    plt.axvline(a,linewidth=2,color='k')
    plt.axvline(b,linestyle='dashed',color='k')
    plt.axvline(c,linestyle='dashed',color='k')
    plt.ylim(0,1)
    plt.xlim(min(var),max(var))
    plt.show()
    print('%s =  %.3f +%.3f -%.3f' % (varname, a, c-a, a-b))
    return

# Reading in all the things is SLOW so I just read some in manually for now
labels = ['M1', 'M2', 'R1', 'R2', 'T1', 'T2/T1', 'R1/a', 'R2/a', 'q1 star1', 'q2 star1', 'q1 star2', 'q2 star2']
M1s, M2s, R1s, R2s = np.loadtxt(starparmfile, usecols=(0,1,2,3), unpack=True)
T1s, temprats, R1as, R2as, q1star1s, q2star1s, q1star2s, q2star2s = np.loadtxt(
    fitparmfile, usecols=(5,6,7,8,10,11,12,13), unpack=True)

allpars = [M1s[4000:], M2s[4000:], R1s[4000:], R2s[4000:], T1s[4000:], temprats[4000:], 
    R1as[4000:], R2as[4000:], q1star1s[4000:], q2star1s[4000:], q1star2s[4000:], q2star2s[4000:]]

print('Finished reading stuff in, generating a plot...')

fig = plt.figure(1, figsize=(15,10))
ymin = 0
ymax = 100000
for idx, param in enumerate(allpars):
    ax = fig.add_subplot(3, 4, idx)
    ax.set_yticklabels(())
    ax.set_ylim([ymin,ymax])
    #ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    for label in ax.get_xticklabels()[::2]: # hide every other tick label
        label.set_visible(False)
    xval = param
    xmin = np.min(param)
    xmax = np.max(param)
    n, bins, patches = plt.hist(xval, 100)
    bincenters = 0.5*(bins[1:] + bins[:-1])
    (mu, sigma) = norm.fit(xval)
    print(mu, sigma)
    #y = mlab.normpdf(bincenters, mu, sigma)
    #plt.plot(bincenters, y, 'r--')
    #plt.xlabel(labels[idx])
    plt.text(xmin + 0.1*(np.abs(xmax-xmin)), 0.8*ymax, labels[idx], size=20)
plt.show()

print('Run this interactively and type cdferrplot(var, varname) for now.')