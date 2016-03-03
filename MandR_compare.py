from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
'''
Plots a comparison of M,R from ELC vs M,R from JKTEBOP.
'''
red = '#e34a33' # red for RG
yel = '#fdbb84' # yellow for MS

infile = '../../RG_ELCmodeling/RGEB_MandR_results.txt'

KICs = str(np.loadtxt(infile, usecols=(0,), comments='#', unpack=True))

ELC_MRGs = np.loadtxt(infile, usecols=(1,2,3), comments='#', unpack=True)
ELC_RRGs = np.loadtxt(infile, usecols=(4,5,6), comments='#', unpack=True)
ELC_MMSs = np.loadtxt(infile, usecols=(7,8,9), comments='#', unpack=True)
ELC_RMSs = np.loadtxt(infile, usecols=(10,11,12), comments='#', unpack=True)

BOP_MRGs = np.loadtxt(infile, usecols=(13,14,15), comments='#', unpack=True)
BOP_RRGs = np.loadtxt(infile, usecols=(16,17,18), comments='#', unpack=True)
BOP_MMSs = np.loadtxt(infile, usecols=(19,20,21), comments='#', unpack=True)
BOP_RMSs = np.loadtxt(infile, usecols=(22,23,24), comments='#', unpack=True)

fig = plt.figure()

fig.add_subplot(2,2,1)
plt.axis([0.6,1.6,0.6,1.6])
plt.plot([0,20], [0,20], color='k', ls=':')
plt.xlabel('ELC RG Mass')
plt.ylabel('JKTEBOP RG Mass')
plt.errorbar(ELC_MRGs[0], BOP_MRGs[0], 
    xerr=[ELC_MRGs[2], ELC_MRGs[1]], yerr=[BOP_MRGs[2], BOP_MRGs[1]], 
    ls='None', marker='o', color=red, markersize=8)

fig.add_subplot(2,2,2)
plt.axis([6,15,6,15])
plt.plot([0,20], [0,20], color='k', ls=':')
plt.xlabel('ELC RG Radius')
plt.ylabel('JKTEBOP RG Radius')
plt.errorbar(ELC_RRGs[0], BOP_RRGs[0], 
    xerr=[ELC_RRGs[2], ELC_RRGs[1]], yerr=[BOP_RRGs[2], BOP_RRGs[1]], 
    ls='None', marker='o', color=red, markersize=8)

fig.add_subplot(2,2,3)
plt.axis([0.6,1.6,0.6,1.6])
plt.plot([0,20], [0,20], color='k', ls=':')
plt.xlabel('ELC MS Mass')
plt.ylabel('JKTEBOP MS Mass')
plt.errorbar(ELC_MMSs[0], BOP_MMSs[0], 
    xerr=[ELC_MMSs[2], ELC_MMSs[1]], yerr=[BOP_MMSs[2], BOP_MMSs[1]], 
    ls='None', marker='o', color=yel, markersize=8)

fig.add_subplot(2,2,4)
plt.axis([0.7,2,0.7,2])
plt.plot([0,20], [0,20], color='k', ls=':')
plt.xlabel('ELC MS Radius')
plt.ylabel('JKTEBOP MS Radius')
plt.errorbar(ELC_RMSs[0], BOP_RMSs[0], 
    xerr=[ELC_RMSs[2], ELC_RMSs[1]], yerr=[BOP_RMSs[2], BOP_RMSs[1]], 
    ls='None', marker='o', color=yel, markersize=8)


plt.show()