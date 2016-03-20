from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
'''
Plots a comparison of M,R from ELC vs M,R from asteroseismology.
'''
red = '#e34a33' # red for RG

infile1 = '../../RG_ELCmodeling/RGEB_MandR_results.txt'
infile2 = '../../seismic_results.txt'
infile3 = '../../RGEB_ELCdensities.txt'

KIC1s = np.loadtxt(infile1, usecols=(0,), comments='#', unpack=True)
KIC2s = np.loadtxt(infile2, usecols=(0,), comments='#', unpack=True)
KIC3s = np.loadtxt(infile3, usecols=(0,), comments='#', unpack=True)
for idx, KIC in enumerate(KIC1s):
    if KIC != KIC2s[idx] or KIC != KIC3s[idx]:
        print('WARNING: input systems differ between files')

# load dynamic values here
ELC_MRGs = np.loadtxt(infile1, usecols=(1,2,3), comments='#', unpack=True)
ELC_RRGs = np.loadtxt(infile1, usecols=(4,5,6), comments='#', unpack=True)
ELCdensities = np.loadtxt(infile3, usecols=(1,2), comments='#', unpack=True)
ELCgravities = np.loadtxt(infile3, usecols=(3,4), comments='#', unpack=True)

# load seismic values here
SEI_MRGs = np.loadtxt(infile2, usecols=(1,2), comments='#', unpack=True)
SEI_RRGs = np.loadtxt(infile2, usecols=(3,4), comments='#', unpack=True)
SEIdensities = np.loadtxt(infile2, usecols=(7,8), comments='#', unpack=True)
SEIgravities = np.loadtxt(infile2, usecols=(5,6), comments='#', unpack=True)

#print(ELCdensities[0])
#print(SEIdensities)

fig = plt.figure()

fig.add_subplot(2,2,1)
plt.axis([0.6,2.5,0.6,2.5])
plt.plot([0,20], [0,20], color='k', ls=':')
plt.xlabel('Dynamic RG Mass')
plt.ylabel('Seismic RG Mass')
plt.errorbar(ELC_MRGs[0], SEI_MRGs[0], 
    xerr=[ELC_MRGs[2], ELC_MRGs[1]], yerr=SEI_MRGs[1], 
    ls='None', marker='o', color=red, markersize=8)
for idx, KIC in enumerate(KIC1s):
    plt.annotate(str(int(KIC)), xy=(ELC_MRGs[0][idx], SEI_MRGs[0][idx]), 
                 xytext=(ELC_MRGs[0][idx]+0.05, SEI_MRGs[0][idx]+0.05), size=10)

fig.add_subplot(2,2,2)
plt.axis([6.1,16,6.1,16])
plt.plot([0,20], [0,20], color='k', ls=':')
plt.xlabel('Dynamic RG Radius')
plt.ylabel('Seismic RG Radius')
plt.errorbar(ELC_RRGs[0], SEI_RRGs[0], 
    xerr=[ELC_RRGs[2], ELC_RRGs[1]], yerr=SEI_RRGs[1], 
    ls='None', marker='o', color=red, markersize=8)
for idx, KIC in enumerate(KIC1s):
    plt.annotate(str(int(KIC)), xy=(ELC_RRGs[0][idx], SEI_RRGs[0][idx]), 
                 xytext=(ELC_RRGs[0][idx]+0.2, SEI_RRGs[0][idx]+0.2), size=10)

fig.add_subplot(2,2,3)
plt.axis([0,5,0,5])
plt.plot([0,20], [0,20], color='k', ls=':')
plt.xlabel('Dynamic RG Density')
plt.ylabel('Seismic RG Density')
plt.errorbar(ELCdensities[0], SEIdensities[0], 
    xerr=ELCdensities[1], yerr=SEIdensities[1], 
    ls='None', marker='o', color=red, markersize=8)
for idx, KIC in enumerate(KIC1s):
    plt.annotate(str(int(KIC)), xy=(ELCdensities[0][idx], SEIdensities[0][idx]), 
                 xytext=(ELCdensities[0][idx]+0.1, SEIdensities[0][idx]+0.1), size=10)

fig.add_subplot(2,2,4)
plt.axis([2.01,3.5,2.01,3.5])
plt.plot([0,20], [0,20], color='k', ls=':')
plt.xlabel('Dynamic RG $\log g$')
plt.ylabel('Seismic RG $\log g$')
plt.errorbar(ELCgravities[0], SEIgravities[0], 
    xerr=ELCgravities[1], yerr=SEIgravities[1], 
    ls='None', marker='o', color=red, markersize=8)
for idx, KIC in enumerate(KIC1s):
    plt.annotate(str(int(KIC)), xy=(ELCgravities[0][idx], SEIgravities[0][idx]), 
                 xytext=(ELCgravities[0][idx]+0.05, SEIgravities[0][idx]+0.05), size=10)


plt.show()