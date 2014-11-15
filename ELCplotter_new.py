from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import IndexLocator, FormatStrFormatter
'''
Meredith Rawls, Sep 2014
Updated, somewhat friendlier than 'ELCplotter.py'
Plotting routine for geneticELC output.
Save this program in a directory where you've run geneticELC or markovELC, and run it.
It will make a plot that has both light curve data w/fit and RV data w/fit.
There are also residuals in the plots!
'''
# Colors for plots. Selected with help from colorbrewer.
red = '#e34a33' # red, star 1
yel = '#fdbb84' # yellow, star 2

# Read in everything
f1 = open('modelU.mag')
f2 = open('ELCdataU.fold')
f3 = open('star1.RV')
f4 = open('star2.RV')
f5 = open('ELCdataRV1.fold')
f6 = open('ELCdataRV2.fold')

phase_mod,mag_mod = np.loadtxt(f1, comments='#', dtype=np.float64, usecols=(0,1), unpack=True)
phase_dat,mag_dat = np.loadtxt(f2, comments='#', dtype=np.float64, usecols=(0,1), unpack=True)
phase_rv1,rv1 = np.loadtxt(f3, comments='#', dtype=np.float64, usecols=(0,1), unpack=True)
phase_rv2,rv2 = np.loadtxt(f4, comments='#', dtype=np.float64, usecols=(0,1), unpack=True)
phase_rv1dat,rv1dat,rv1err = np.loadtxt(f5, comments='#', dtype=np.float64, usecols=(0,1,2), unpack=True)
phase_rv2dat,rv2dat,rv2err = np.loadtxt(f6, comments='#', dtype=np.float64, usecols=(0,1,2), unpack=True)

f1.close(); f2.close(); f3.close(); f4.close(); f5.close(); f6.close()

# OPTIONAL ADJUSTMENT B/C FINAL GENETICELC RV MODEL OUTPUT IS SHIFTED BY GAMMA
#gamma = input("Enter gamma adjustment (0 for none): ")
#rv1 = rv1 + gamma
#rv2 = rv2 + gamma

print ("Done reading data!")

if np.abs(np.median(mag_mod) - np.median(mag_dat)) > 1:
	print('Adjusting magnitude of model light curve...')
	mag_mod = mag_mod + (np.median(mag_dat) - np.median(mag_mod))

# Interpolate model onto data phase grid, for residuals
newmag_model = np.interp(phase_dat, phase_mod, mag_mod)
newrv1 = np.interp(phase_rv1dat, phase_rv1, rv1)
newrv2 = np.interp(phase_rv2dat, phase_rv2, rv2)

lcresid = mag_dat - newmag_model
rv1resid = rv1dat - newrv1
rv2resid = rv2dat - newrv2

print ("Done interpolating!")

# Make plots
# First, define some handy global parameters for the plots
phasemin = 0
phasemax = 1
magdim = 9.52			# remember magnitudes are backwards, dangit
magbright = 9.251
rvmin = -55
rvmax = 55
primary_phasemin = 0.48
primary_phasemax = 0.52
secondary_phasemin = 0.194
secondary_phasemax = 0.234
magresid_min = 0.010	# remember magnitudes are backwards, dangit
magresid_max = -0.010
rvresid_min = -6
rvresid_max = 6

# Light curve
ax1 = plt.subplot2grid((12,1),(4,0), rowspan=3)
plt.axis([phasemin, phasemax, magdim, magbright])
plt.tick_params(axis='both', which='major')
plt.plot(phase_dat, mag_dat, color=red, marker='.', ls='None', ms=1, mew=0) #lc data
plt.plot(phase_mod, mag_mod, 'k', lw=1.5) #lc model
plt.plot(phase_dat+1, mag_dat, color=red, marker='.', ls='None', ms=1, mew=0) #lc data 1-2
plt.plot(phase_mod+1, mag_mod, 'k', lw=1.5, label='ELC Model') #lc model 1-2
ax1.set_ylabel('Magnitude', size=18)
ax1.set_xticklabels([])

# Radial velocities
ax2 = plt.subplot2grid((12,1),(1,0), rowspan=3)
plt.subplots_adjust(wspace = 0.0001, hspace=0.0001)
plt.axis([phasemin, phasemax, rvmin, rvmax])
plt.errorbar(phase_rv1dat, rv1dat, yerr=rv1err, marker='o', color=red, mec=red, ls='None') #rv1 data
plt.errorbar(phase_rv2dat, rv2dat, yerr=rv2err, marker='o', color=yel, mec=yel, ls='None') #rv2 data
plt.plot(phase_rv1, rv1, color='k', lw=1.5) #rv1 model
plt.plot(phase_rv2, rv2, color='k', lw=1.5) #rv2 model
plt.errorbar(phase_rv1dat+1, rv1dat, yerr=rv1err, marker='o', color=red, mec=red, ls='None') #rv1 data 1-2
plt.errorbar(phase_rv2dat+1, rv2dat, yerr=rv2err, marker='o', color=yel, mec=yel, ls='None') #rv2 data 1-2
plt.plot(phase_rv1+1, rv1, color='k', lw=1.5) #rv1 model 1-2
plt.plot(phase_rv2+1, rv2, color='k', lw=1.5) #rv2 model 1-2
ax2.set_ylabel('Radial Velocity (km s$^{-1}$)', size=18)
ax2.set_xticklabels([])

# Light curve residuals
axr1 = plt.subplot2grid((12,1),(7,0))
axr1.axis([phasemin, phasemax, magresid_min, magresid_max])
axr1.set_yticks([-0.006, 0, 0.006])
axr1.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
plt.axhline(y=0, xmin=phasemin, xmax=phasemax, color='0.75', ls=':')
plt.plot(phase_dat, lcresid, color=red, marker='.', ls='None', ms=1, mew=0) #lc residual
plt.plot(phase_dat+1, lcresid, color=red, marker='.', ls='None', ms=1, mew=0) #lc residual 1-2
#axr1.set_xticklabels([])

# Radial velocity residuals
axr2 = plt.subplot2grid((12,1),(0,0))
axr2.axis([phasemin, phasemax, rvresid_min, rvresid_max])
axr2.set_yticks([-4,0,4])
plt.axhline(y=0, xmin=phasemin, xmax=phasemax, color='0.75')
plt.errorbar(phase_rv1dat, rv1resid, yerr=rv1err, marker='o', color=red, mec=red, ls='None') #rv1 residual
plt.errorbar(phase_rv2dat, rv2resid, yerr=rv2err, marker='o', color=yel, mec=yel, ls='None') #rv2 residual
plt.errorbar(phase_rv1dat+1, rv1resid, yerr=rv1err, marker='o', color=red, mec=red, ls='None') #rv1 residual 1-2
plt.errorbar(phase_rv2dat+1, rv2resid, yerr=rv2err, marker='o', color=yel, mec=yel, ls='None') #rv2 residual 1-2
#plt.xlabel('Orbital Phase (conjunction at $\phi = 0.5$)', size=20) # EXTRA LABEL
axr2.set_xticklabels([])

# Zoom-in of shallower (secondary) eclipse
ax3 = plt.subplot2grid((12,2),(9,0), rowspan=2)
plt.axis([secondary_phasemin, secondary_phasemax, magdim, magbright])
ax3.xaxis.set_major_locator(IndexLocator(0.01, 0.20))
ax3.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.plot(phase_dat, mag_dat, color=yel, marker='.', ls='None', ms=2, mew=0) #lc data
plt.plot(phase_mod, mag_mod, color='k', lw=1.5) #lc model
ax3.set_ylabel('Magnitude')
ax3.set_xticklabels([])

# Zoom-in of deeper (primary) eclipse
ax4 = plt.subplot2grid((12,2),(9,1), rowspan=2)
ax4.set_yticklabels([])
plt.axis([primary_phasemin, primary_phasemax, magdim, magbright])
ax4.xaxis.set_major_locator(IndexLocator(0.01, 0.49))
ax4.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.plot(phase_dat, mag_dat, color=red, marker='.', ls='None', ms=2, mew=0) #lc data
plt.plot(phase_mod, mag_mod, color='k', lw=1.5) #lc model
ax4.set_xticklabels([])

# Zoom plot residuals, shallower (secondary) eclipse
axr3 = plt.subplot2grid((12,2),(11,0))
plt.axis([secondary_phasemin, secondary_phasemax, magresid_min, magresid_max])
axr3.xaxis.set_major_locator(IndexLocator(0.01, 0.20))
axr3.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
axr3.set_yticks([-0.006, 0, 0.006])
axr3.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
plt.axhline(y=0, xmin=0, xmax=2, color='0.75')
plt.plot(phase_dat, lcresid, color=red, marker='.', ls='None', ms=2, mew=0) #lc residual

# Zoom plot residuals, deeper (primary) eclipse
axr4 = plt.subplot2grid((12,2),(11,1))
plt.axis([primary_phasemin, primary_phasemax, magresid_min, magresid_max])
axr4.xaxis.set_major_locator(IndexLocator(0.01, 0.49))
axr4.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
axr4.set_yticks([-0.006, 0, 0.006])
axr4.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
plt.axhline(y=0, xmin=0, xmax=2, color='0.75')
plt.plot(phase_dat, lcresid, color=red, marker='.', ls='None', ms=2, mew=0) #lc residual
axr4.set_yticklabels([])

# Labels using overall figure as a reference
plt.figtext(0.5, 0.04, 'Orbital Phase (conjunction at $\phi = 0.5$)', ha='center', va='center', size=25)
plt.figtext(0.135, 0.18, 'Secondary')
plt.figtext(0.525, 0.18, 'Primary')
plt.figtext(0.06, 0.86, '$\Delta$')
plt.figtext(0.04, 0.395, '$\Delta$')
plt.figtext(0.04, 0.125, '$\Delta$')
ax1.legend(loc='lower right', frameon=False, prop={'size':20})

print ("Done preparing plot!")

plt.show()