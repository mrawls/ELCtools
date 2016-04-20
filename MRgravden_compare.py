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

# print the values being plotted, just to double check
#print(ELC_MRGs[0], SEI_MRGs[0])
#print(ELC_RRGs[0], SEI_RRGs[0])
#print(ELCdensities[0], SEIdensities[0])
#print(ELCgravities[0], SEIgravities[0])

fig = plt.figure()

# masses
ax1 = plt.subplot2grid((13,13),(0,1), rowspan=5, colspan=5)
ax1.set_xticklabels([])
plt.subplots_adjust(wspace = 0, hspace=0)
#ax1 = fig.add_subplot(2,2,1)
plt.axis([0.6,2.5,0.6,2.5])
plt.plot([0,20], [0,20], color='k', ls=':')
plt.ylabel(r'Seismic Mass ($M_{\odot}$)')
plt.errorbar(ELC_MRGs[0], SEI_MRGs[0], 
    xerr=[ELC_MRGs[2], ELC_MRGs[1]], yerr=SEI_MRGs[1], 
    ls='None', marker='o', color=red, markersize=8)
# systems:   1000,  578,  703,  924, 997
offsetx1s = [ 0.00,-0.20, 0.00,-0.42, 0.05]
offsety1s = [-0.15, 0.10,-0.15, 0.00, 0.02]
for idx, KIC in enumerate(KIC1s):
    plt.annotate(str(int(KIC)), xy=(ELC_MRGs[0][idx], SEI_MRGs[0][idx]), size=14,
                 xytext=(ELC_MRGs[0][idx]+offsetx1s[idx], SEI_MRGs[0][idx]+offsety1s[idx]))
# residuals
ax1r = plt.subplot2grid((13,13),(5,1), colspan=5)
ax1r.set_xlim([0.6,2.5])
ax1r.set_ylim([-0.2,0.5])
ax1r.set_yticks([-0.2,0,0.2,0.4])
plt.xlabel(r'Dynamic Mass ($M_{\odot}$)')
plt.ylabel('$\%$')
residuals = (SEI_MRGs[0] - ELC_MRGs[0])/ELC_MRGs[0]
plt.errorbar(ELC_MRGs[0], residuals, yerr=SEI_MRGs[1]/SEI_MRGs[0], 
             ls='None', marker='o', color=red, markersize=8)
plt.axhline(y=0, ls=':', color='k')

# radii
ax2 = plt.subplot2grid((13,13),(0,8), rowspan=5, colspan=5)
ax2.set_xticklabels([])
plt.axis([6.1,16,6.1,16])
plt.plot([0,20], [0,20], color='k', ls=':')
plt.ylabel(r'Seismic Radius ($R_{\odot}$)')
plt.errorbar(ELC_RRGs[0], SEI_RRGs[0], 
    xerr=[ELC_RRGs[2], ELC_RRGs[1]], yerr=SEI_RRGs[1], 
    ls='None', marker='o', color=red, markersize=8)
# systems:   1000,  578,  703,  924, 997
offsetx2s = [-2.20,-2.20,-2.20, 0.20,-1.40]
offsety2s = [ 0.20, 0.20, 0.20,-0.60, 0.40]
for idx, KIC in enumerate(KIC1s):
    plt.annotate(str(int(KIC)), xy=(ELC_RRGs[0][idx], SEI_RRGs[0][idx]), size=14,
                 xytext=(ELC_RRGs[0][idx]+offsetx2s[idx], SEI_RRGs[0][idx]+offsety2s[idx]))
# residuals
ax2r = plt.subplot2grid((13,13),(5,8), colspan=5)
ax2r.set_xlim([6.1,16])
ax2r.set_ylim([-0.1,0.2])
ax2r.set_yticks([-0.1,0,0.1,0.2])
plt.xlabel(r'Dynamic Radius ($R_{\odot}$)')
plt.ylabel('$\%$')
residuals = (SEI_RRGs[0] - ELC_RRGs[0])/ELC_RRGs[0]
plt.errorbar(ELC_RRGs[0], residuals, yerr=SEI_RRGs[1]/SEI_RRGs[0], 
             ls='None', marker='o', color=red, markersize=8)
plt.axhline(y=0, ls=':', color='k')

# densities
ax3 = plt.subplot2grid((13,13),(7,1), rowspan=5, colspan=5)
ax3.set_xticklabels([])
plt.axis([0.01,4,0.01,4])
plt.plot([0,20], [0,20], color='k', ls=':')
plt.ylabel(r'Seismic Density ($\bar{\rho}_{\odot} \times 10^{-3}$)')
plt.errorbar(ELCdensities[0], SEIdensities[0], 
    xerr=ELCdensities[1], yerr=SEIdensities[1], 
    ls='None', marker='o', color=red, markersize=8)
# systems:   1000,  578,  703,  924, 997
offsetx3s = [-0.2,  0.1,  0.1, -0.9, 0.1]
offsety3s = [-0.35, 0.1, -0.1, -0.2, 0.1]
for idx, KIC in enumerate(KIC1s):
    plt.annotate(str(int(KIC)), xy=(ELCdensities[0][idx], SEIdensities[0][idx]), size=14,
                 xytext=(ELCdensities[0][idx]+offsetx3s[idx], SEIdensities[0][idx]+offsety3s[idx]))
# residuals
ax3r = plt.subplot2grid((13,13),(12,1), colspan=5)
ax3r.set_xlim([0.01,4])
ax3r.set_ylim([-0.2,0.1])
ax3r.set_yticks([-0.2,-0.1,0])
plt.xlabel(r'Dynamic Density ($\bar{\rho}_{\odot} \times 10^{-3}$)')
plt.ylabel('$\%$')
residuals = (SEIdensities[0] - ELCdensities[0])/ELCdensities[0]
plt.errorbar(ELCdensities[0], residuals, yerr=SEIdensities[1]/SEIdensities[0], 
             ls='None', marker='o', color=red, markersize=8)
plt.axhline(y=0, ls=':', color='k')

# gravities
ax4 = plt.subplot2grid((13,13),(7,8), rowspan=5, colspan=5)
ax4.set_xticklabels([])
plt.axis([2.01,3.4,2.01,3.4])
plt.plot([0,20], [0,20], color='k', ls=':')
plt.ylabel(r'Seismic $\log g$ (cgs)')
plt.errorbar(ELCgravities[0], SEIgravities[0], 
    xerr=ELCgravities[1], yerr=SEIgravities[1], 
    ls='None', marker='o', color=red, markersize=8)
# systems:   1000,  578,  703,  924, 997
offsetx4s = [ 0.01, 0.01, 0.01, 0.01, 0.01]
offsety4s = [-0.10, 0.03, 0.03,-0.08,-0.08]
for idx, KIC in enumerate(KIC1s):
    plt.annotate(str(int(KIC)), xy=(ELCgravities[0][idx], SEIgravities[0][idx]), size=14,
                 xytext=(ELCgravities[0][idx]+offsetx4s[idx], SEIgravities[0][idx]+offsety4s[idx]))
# what the heck, let's throw spectroscopic log g on this plot too
specloggs =    [2.65, 2.65, 2.46, 3.33, 3.14]
specloggerrs = [0.10, 0.21, 0.17, 0.37, 0.12]
plt.errorbar(ELCgravities[0], specloggs, xerr=ELCgravities[1], yerr=specloggerrs, 
             ls='None', marker='o', color='0.75', markersize=8)
fig.text(0.91, 0.395, '(Spectroscopic $\log g$)', size=20, rotation='vertical', color='k')

# residuals
ax4r = plt.subplot2grid((13,13),(12,8), colspan=5)
ax4r.set_xlim([2.01,3.4])
ax4r.set_ylim([-0.05,0.08])
ax4r.set_yticks([-0.04,0,0.04,0.08])
plt.xlabel(r'Dynamic $\log g$ (cgs)')
plt.ylabel('$\Delta$')
residuals = SEIgravities[0] - ELCgravities[0]
plt.errorbar(ELCgravities[0], residuals, yerr=SEIgravities[1]/SEIgravities[0], 
             ls='None', marker='o', color=red, markersize=8)
plt.axhline(y=0, ls=':', color='k')

plt.show()