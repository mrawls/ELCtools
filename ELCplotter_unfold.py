from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import IndexLocator, FormatStrFormatter
'''
Meredith Rawls, Dec 2014
Updated, somewhat friendlier than 'ELCplotter.py'
Plotting routine for geneticELC / markovELC output.
Save this program in a directory where you've run geneticELC or markovELC, and run it.
It will make a plot that has both light curve data w/fit and RV data w/fit.
There are also residuals in the plots!

***IMPORTANT***
This version assumes the files above are NOT yet folded in phase, and are in time.
This would happen if you are using ELCgap.inp, or anytime when ELC.inp has itime = 2.
So we need to fold them.
'''
# Colors for plots. Selected with help from colorbrewer.
red = '#e34a33' # red, star 1
yel = '#fdbb84' # yellow, star 2

# Columns in fitparm file that correspond to T0 and Period
tconj_col = 0
porb_col = 15

# Read in everything
f1 =        '../../RG_ELCmodeling/9246715/demcmc001/modelU.mag'
f2 =        '../../RG_ELCmodeling/9246715/demcmc001/ELCdataU.fold'
f3 =        '../../RG_ELCmodeling/9246715/demcmc001/star1.RV'
f4 =        '../../RG_ELCmodeling/9246715/demcmc001/star2.RV'
f5 =        '../../RG_ELCmodeling/9246715/demcmc001/ELCdataRV1.fold'
f6 =        '../../RG_ELCmodeling/9246715/demcmc001/ELCdataRV2.fold'
fitparm =   '../../RG_ELCmodeling/9246715/demcmc001/fitparm.all'
#ELCoutfile = '../../RG_ELCmodeling/9246715/newtrial3/ELC.out'

phase_mod,mag_mod = np.loadtxt(f1, comments='#', dtype=np.float64, usecols=(0,1), unpack=True)
phase_dat,mag_dat = np.loadtxt(f2, comments='#', dtype=np.float64, usecols=(0,1), unpack=True)
phase_rv1,rv1 = np.loadtxt(f3, comments='#', dtype=np.float64, usecols=(0,1), unpack=True)
phase_rv2,rv2 = np.loadtxt(f4, comments='#', dtype=np.float64, usecols=(0,1), unpack=True)
phase_rv1dat,rv1dat,rv1err = np.loadtxt(f5, comments='#', dtype=np.float64, usecols=(0,1,2), unpack=True)
phase_rv2dat,rv2dat,rv2err = np.loadtxt(f6, comments='#', dtype=np.float64, usecols=(0,1,2), unpack=True)

# FUNCTION TO FOLD STUFF so phases are actually phases ... and then sort all the arrays.
def phasecalc(times, period=100, BJD0=2454833):
	phases = []
	cycles = []
	for i in range(0, len(times)):
		fracP = (times[i] - BJD0) / period
		if fracP < 0:
			phases.append(fracP % 1)
			cycles.append(int(fracP))
		else:
			phases.append(fracP % 1)
			cycles.append(int(fracP) + 1)
		#print(fracP, phases[i])
	return np.array(phases)

# GET PERIOD AND T0 from files
#with open(ELCoutfile) as f:
#	for i, row in enumerate(f):
#		if i == 27: # 28th row
#			columns = row.split()
#			period = float(columns[0]) # 1st column
#		#if i == 38: # 39th row, i.e. T0	# this one has a funny zeropoint (ok if circular)
#		if i == 133: # 134th row, i.e. Tconj # this one puts primary eclipse at phase 0
#			columns = row.split()
#			Tconj = float(columns[0]) #1st column

periods, tconjs = np.loadtxt(fitparm, usecols=(porb_col, tconj_col), unpack=True)
period = np.median(periods)
Tconj = np.median(tconjs)

print(period, Tconj)
Tconj = Tconj + 0.5*period

phase_mod = phasecalc(phase_mod, period=period, BJD0=Tconj)
phase_dat = phasecalc(phase_dat, period=period, BJD0=Tconj)
phase_rv1 = phasecalc(phase_rv1, period=period, BJD0=Tconj)
phase_rv2 = phasecalc(phase_rv2, period=period, BJD0=Tconj)
phase_rv1dat = phasecalc(phase_rv1dat, period=period, BJD0=Tconj)
phase_rv2dat = phasecalc(phase_rv2dat, period=period, BJD0=Tconj)

p1 = phase_mod.argsort()
p2 = phase_dat.argsort()
p3 = phase_rv1.argsort()
p4 = phase_rv2.argsort()
p5 = phase_rv1dat.argsort()
p6 = phase_rv2dat.argsort()

phase_mod = phase_mod[p1]
phase_dat = phase_dat[p2]
phase_rv1 = phase_rv1[p3]
phase_rv2 = phase_rv2[p4]
phase_rv1dat = phase_rv1dat[p5]
phase_rv2dat = phase_rv2dat[p6]

mag_mod = mag_mod[p1]
mag_dat = mag_dat[p2]
rv1 = rv1[p3]
rv2 = rv2[p4]
rv1dat = rv1dat[p5]
rv2dat = rv2dat[p6]


# OPTIONAL ADJUSTMENT B/C FINAL GENETICELC RV MODEL OUTPUT IS SHIFTED BY GAMMA
#gamma = input("Enter gamma adjustment (0 for none): ")
#rv1 = rv1 + gamma
#rv2 = rv2 + gamma

print ("Done reading (and folding) data!")

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
rvmin = -60
rvmax = 50
primary_phasemin = 0.48 #0.09 #0.48
primary_phasemax = 0.52 #0.14 #0.52
secondary_phasemin = 0.194 #0.80 #0.194
secondary_phasemax = 0.234 #0.85 #0.234
magresid_min = 0.006	# remember magnitudes are backwards, dangit
magresid_max = -0.006
rvresid_min = -2.3
rvresid_max = 2.3

# Light curve
ax1 = plt.subplot2grid((12,1),(4,0), rowspan=3)
plt.axis([phasemin, phasemax, magdim, magbright])
plt.tick_params(axis='both', which='major')
plt.plot(phase_dat, mag_dat, color=red, marker='.', ls='None', ms=6, mew=0) #lc data
plt.plot(phase_mod, mag_mod, 'k', lw=1.5, label='ELC Model') #lc model
ax1.set_ylabel('Magnitude', size=18)
ax1.set_xticklabels([])

# Radial velocities
ax2 = plt.subplot2grid((12,1),(1,0), rowspan=3)
plt.subplots_adjust(wspace = 0.0001, hspace=0.0001)
plt.axis([phasemin, phasemax, rvmin, rvmax])
plt.errorbar(phase_rv1dat, rv1dat, yerr=rv1err, marker='o', color=red, ms=9, mec='None', ls='None') #rv1 data
plt.errorbar(phase_rv2dat, rv2dat, yerr=rv2err, marker='o', color=yel, ms=9, mec='None', ls='None') #rv2 data
plt.plot(phase_rv1, rv1, color='k', lw=1.5) #rv1 model
plt.plot(phase_rv2, rv2, color='k', lw=1.5) #rv2 model
ax2.set_ylabel('Radial Velocity (km s$^{-1}$)', size=18)
ax2.set_xticklabels([])

# Light curve residuals
axr1 = plt.subplot2grid((12,1),(7,0))
axr1.axis([phasemin, phasemax, magresid_min, magresid_max])
axr1.set_yticks([-0.004, 0, 0.004])
plt.axhline(y=0, xmin=phasemin, xmax=phasemax, color='0.75', ls=':')
plt.plot(phase_dat, lcresid, color=red, marker='.', ls='None', ms=4, mew=0) #lc residual

# Radial velocity residuals
axr2 = plt.subplot2grid((12,1),(0,0))
axr2.axis([phasemin, phasemax, rvresid_min, rvresid_max])
axr2.set_yticks([-2,0,2])
plt.axhline(y=0, xmin=phasemin, xmax=phasemax, color='0.75', ls=':')
plt.errorbar(phase_rv1dat, rv1resid, yerr=rv1err, marker='o', color=red, ms=9, mec='None', ls='None') #rv1 residual
plt.errorbar(phase_rv2dat, rv2resid, yerr=rv2err, marker='o', color=yel, ms=9, mec='None', ls='None') #rv2 residual
#plt.xlabel('Orbital Phase (conjunction at $\phi = 0.5$)', size=20) # EXTRA LABEL
axr2.set_xticklabels([])

# Zoom-in of shallower (secondary) eclipse
ax3 = plt.subplot2grid((12,2),(9,0), rowspan=2)
plt.axis([secondary_phasemin, secondary_phasemax, magdim, magbright])
ax3.set_xticks([0.20, 0.21, 0.22, 0.23])
plt.plot(phase_dat, mag_dat, color=yel, marker='.', ls='None', ms=6, mew=0) #lc data
plt.plot(phase_mod, mag_mod, color='k', lw=1.5) #lc model
ax3.set_ylabel('Magnitude')
ax3.set_xticklabels([])
#ax3.set_yticklabels([])

# Zoom-in of deeper (primary) eclipse
ax4 = plt.subplot2grid((12,2),(9,1), rowspan=2)
plt.axis([primary_phasemin, primary_phasemax, magdim, magbright])
ax4.set_xticks([0.49, 0.50, 0.51, 0.52])
plt.plot(phase_dat, mag_dat, color=red, marker='.', ls='None', ms=6, mew=0) #lc data
plt.plot(phase_mod, mag_mod, color='k', lw=1.5) #lc model
ax4.set_xticklabels([])
ax4.set_yticklabels([])

# Zoom plot residuals, shallower (secondary) eclipse
axr3 = plt.subplot2grid((12,2),(11,0))
plt.axis([secondary_phasemin, secondary_phasemax, magresid_min, magresid_max])
axr3.set_yticks([-0.004, 0, 0.004])
axr3.set_xticks([0.20, 0.21, 0.22, 0.23])
plt.axhline(y=0, xmin=0, xmax=2, color='0.75', ls=':')
plt.plot(phase_dat, lcresid, color=red, marker='.', ls='None', ms=4, mew=0) #lc residual
#axr3.set_yticklabels([])

# Zoom plot residuals, deeper (primary) eclipse
axr4 = plt.subplot2grid((12,2),(11,1))
plt.axis([primary_phasemin, primary_phasemax, magresid_min, magresid_max])
axr4.set_yticks([-0.004, 0, 0.004])
axr4.set_xticks([0.49, 0.50, 0.51, 0.52])
plt.axhline(y=0, xmin=0, xmax=2, color='0.75', ls=':')
plt.plot(phase_dat, lcresid, color=red, marker='.', ls='None', ms=4, mew=0) #lc residual
axr4.set_yticklabels([])

# Labels using overall figure as a reference
plt.figtext(0.5, 0.04, 'Orbital Phase (conjunction at $\phi = 0.5$)', ha='center', va='center', size=25)
plt.figtext(0.135, 0.18, 'Secondary')
plt.figtext(0.535, 0.18, 'Primary')
plt.figtext(0.06, 0.86, '$\Delta$')
plt.figtext(0.04, 0.395, '$\Delta$')
plt.figtext(0.04, 0.125, '$\Delta$')
ax1.legend(loc='lower right', frameon=False, prop={'size':20})

print ("Done preparing plot!")

plt.show()
#outfile = 'testplot1.png'
#plt.savefig(outfile)
#print ("Plot saved to %s!" % outfile)