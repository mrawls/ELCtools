from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter, MaxNLocator
import cubehelix # from https://github.com/jradavenport/cubehelix
'''
Meredith Rawls
September 2014

Plots a grid of 2D histograms for a set of fit parameters.
Inspired by http://oceanpython.org/2013/02/25/2d-histogram/

Designed to work with the generation.*** files output by markov/geneticELC,
and/or the markovchainXX.*** files output by markovELC.

	Before running this, you'll want to do something like:
	cat generation.1* > generation.all
	
	Or a set of:
	cat markov01.1* > markov01.txt
	(for each parameter, 01--npar)

Also plots a RV curve (data + fit) and light curve (data + fit) in the upper right.
'''

# Number of fit parameters and their names
filepath = '../../RG_ELCmodeling/9246715/trial4_fromseismo_super/'
folded = False # if you used itime=2 in ELC, set this to False. otherwise, True.
npar = 10 #please make it even for sanity ... actually maybe odd is ok?
parlabels = ['$P_{orb}$', '$e \cos \omega$','$e \sin \omega$','$i$','$M_1$','$T_1$','$T_2/T_1$','$R_1/a$','$R_2/a$','$K_1$']
pars = ['period', 'oc', 'os', 'inclination', 'pm', 'T1', 'temprat', 'q1', 'q2', 'pk'] #ELC names
parindices = [2, 4, 5, 6, 7, 8, 9, 10, 11, 12] #ELC index
nbins = 50 # number of histogram bins (same in both dimensions)

# Files that have stuff in them
gridloop = filepath + 'gridloop.opt'
parfile = filepath + 'generation.all'
lcmodel = filepath + 'modelU.mag'
lcdata = filepath + 'ELCdataU.fold'
rv1model = filepath + 'star1.RV'
rv1data = filepath + 'ELCdataRV1.fold'
rv2model = filepath + 'star2.RV'
rv2data = filepath + 'ELCdataRV2.fold'

# Color map for 2D histogram "contours"
cx = cubehelix.cmap(reverse=True, maxLight=0.7, start=0., rot=0.5)

# Plot colors and LC/RV axis limits
red = '#e34a33' # red, star 1
yel = '#fdbb84' # yellow, star 2
phasemin = 0
phasemax = 1
magdim = 9.54
magbright = 9.21
rvmin = -59
rvmax = 59

# Function to fold stuff if it isn't folded already
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

# Get axis limits for each histogram box from the gridloop.opt file
parmin = []
parmax = []
alllines = []
f1 = open(gridloop)
for idx, line in enumerate(f1):
	line = line.rstrip('\n')
	if idx == 10:
		nfit = int(line)
	for jdx, parname in enumerate(pars):
		if line == pars[jdx]:
			f2 = open(gridloop)
			for newline in f2:
				alllines.append( newline.rstrip('\n') )
			values = alllines[idx+nfit].split()
			parmin.append(values[0])
			parmax.append(values[1])
f1.close()
f2.close()

# Make sure everything is consistent with gridloop.opt
#for idx, par in enumerate(pars):
#	print(par, parmin[idx], parmax[idx], parlabels[idx])

# Read in data that ELC creates for all of the parameters
# For either option, you need to manually concatenate first.
##
## OPTION 1: generation.*** files ##
##
pardata = []
for idx, par in enumerate(pars):
	with open(parfile) as f1:
		pardata.append(np.loadtxt(f1, dtype=np.float64, usecols=(parindices[idx],), unpack=True))
##
## OPTION 2: markovchainXX.*** files ## ... probably not useful
##
#pardata = []
#for idx, par in enumerate(pars):
#	if parindices[idx] < 10: parindex = '0' + str(parindices[idx])
#	else: parindex = str(parindices[idx])
#	parfile = KIC + '/markovchain' + parindex + '.txt'
#	with open(parfile) as f1:
#		pardata.append(np.loadtxt(f1, dtype=np.float64, unpack=True))

# Loop over each plot square
fig1 = plt.figure()
for row in range(0, npar):
	for col in range(0, row+1):
		if row == col: #1D histogram
			ax = fig1.add_subplot(npar,npar,col+row*npar+1)
#			ax.set_xlim([float(parmin[col]),float(parmax[col])])
			ax.xaxis.set_major_locator(MaxNLocator(3)) # no more than 3 ticks per axis
			ax.yaxis.set_major_locator(MaxNLocator(3))
			if row != npar-1:
				ax.set_xticklabels(())
			if col != 0 or row == 0:
				ax.set_yticklabels(())
			plt.subplots_adjust(wspace=0.0, hspace=0.0)
			# x value is the one with index=col
			plt.hist(pardata[col], nbins/4, facecolor=red, edgecolor=None)#, histtype='stepfilled')
		else: #2D histogram
			ax = fig1.add_subplot(npar,npar,col+row*npar+1)
#			ax.set_xlim([float(parmin[col]),float(parmax[col])])
#			ax.set_ylim([float(parmin[row]),float(parmax[row])])
			ax.xaxis.set_major_locator(MaxNLocator(3))
			ax.yaxis.set_major_locator(MaxNLocator(3))
			if row != npar-1:
				ax.set_xticklabels(())
			if col != 0:
				ax.set_yticklabels(())
			# x value is the one with index=col; y value is the one with index=row
			while len(pardata[col]) > len(pardata[row]): # pad pardata[row] with zeros
				pardata[row] = np.append(pardata[row], 0)
			while len(pardata[col]) < len(pardata[row]): # pad pardata[col] with zeros
				pardata[col] = np.append(pardata[col], 0)
			H, xedges, yedges = np.histogram2d(pardata[col], pardata[row], bins=(nbins, nbins))	
			H = np.rot90(H) # H needs to be rotated and flipped, because reasons
			H = np.flipud(H)
			Hmasked = np.ma.masked_where(H==0,H) # mask anything with a value of zero
			plt.pcolormesh(xedges, yedges, Hmasked, cmap=cx)
		if row == npar-1:
			plt.xlabel(parlabels[col], labelpad=30) #x-axis labels
		if col == 0 and row != 0:
			plt.ylabel(parlabels[row], rotation='horizontal', labelpad=30) #y-axis labels

# Open stuff necessary for LC + RV plot
with open(rv1data) as f1:
	phase_rv1dat, rv1dat = np.loadtxt(f1, comments='#', dtype=np.float64, usecols=(0,1), unpack=True)
with open(rv2data) as f1:
	phase_rv2dat, rv2dat = np.loadtxt(f1, comments='#', dtype=np.float64, usecols=(0,1), unpack=True)
with open(rv1model) as f1:
	phase_rv1mod, rv1mod = np.loadtxt(f1, comments='#', dtype=np.float64, usecols=(0,1), unpack=True)
with open(rv2model) as f1:
	phase_rv2mod, rv2mod = np.loadtxt(f1, comments='#', dtype=np.float64, usecols=(0,1), unpack=True)
with open(lcmodel) as f1:
	phase_lcmod, mag_lcmod = np.loadtxt(f1, comments='#', dtype=np.float64, usecols=(0,1), unpack=True)
with open(lcdata) as f1:
	phase_lcdat, mag_lcdat = np.loadtxt(f1, comments='#', dtype=np.float64, usecols=(0,1), unpack=True)

if folded == False:
# FOLD STUFF so phases are actually phases ... and then sort all the arrays.
	with open(filepath+'ELC.out') as f:
		for i, row in enumerate(f):
			if i == 27: # 28th row
				columns = row.split()
				period = float(columns[0]) # 1st column
			#if i == 38: # 39th row, i.e. T0	# this one has a funny zeropoint (ok if circular)
			if i == 133: # 134th row, i.e. Tconj # this one puts primary eclipse at phase 0
				columns = row.split()
				Tconj = float(columns[0]) #1st column
	Tconj = Tconj + 0.5*period
	phase_lcmod = phasecalc(phase_lcmod, period=period, BJD0=Tconj)
	phase_lcdat = phasecalc(phase_lcdat, period=period, BJD0=Tconj)
	phase_rv1mod = phasecalc(phase_rv1mod, period=period, BJD0=Tconj)
	phase_rv2mod = phasecalc(phase_rv2mod, period=period, BJD0=Tconj)
	phase_rv1dat = phasecalc(phase_rv1dat, period=period, BJD0=Tconj)
	phase_rv2dat = phasecalc(phase_rv2dat, period=period, BJD0=Tconj)
	p1 = phase_lcmod.argsort()
	p2 = phase_lcdat.argsort()
	p3 = phase_rv1mod.argsort()
	p4 = phase_rv2mod.argsort()
	p5 = phase_rv1dat.argsort()
	p6 = phase_rv2dat.argsort()
	phase_lcmod = phase_lcmod[p1]
	phase_lcdat = phase_lcdat[p2]
	phase_rv1mod = phase_rv1mod[p3]
	phase_rv2mod = phase_rv2mod[p4]
	phase_rv1dat = phase_rv1dat[p5]
	phase_rv2dat = phase_rv2dat[p6]
	mag_lcmod = mag_lcmod[p1]
	mag_lcdat = mag_lcdat[p2]
	rv1mod = rv1mod[p3]
	rv2mod = rv2mod[p4]
	rv1dat = rv1dat[p5]
	rv2dat = rv2dat[p6]

# Upper right corner: RV plot
fig2 = plt.subplot2grid((npar*2,npar*2), (0,npar+1), rowspan=int((npar-1)/2), colspan=npar-1)
fig2.axis([phasemin,phasemax,rvmin,rvmax])
fig2.set_xticklabels(())
fig2.spines['top'].set_visible(False)
fig2.spines['right'].set_visible(False)
fig2.xaxis.set_ticks_position('bottom')
fig2.yaxis.set_ticks_position('left')
plt.plot(phase_rv1dat, rv1dat, marker='o', color=red, mec=red, ms=6, ls='None')
plt.plot(phase_rv2dat, rv2dat, marker='o', color=yel, mec=yel, ms=6, ls='None')
plt.plot(phase_rv1mod, rv1mod, color='k', lw=1.5)
plt.plot(phase_rv2mod, rv2mod, color='k', lw=1.5, label='ELC model')
fig2.legend(loc='upper right', frameon=False, prop={'size':20})
plt.ylabel('RV (km s$^{-1}$)')

# Upper right corner: Light curve plot
fig3 = plt.subplot2grid((npar*2,npar*2), (int((npar-1)/2),npar+1), rowspan=(npar-1)-int((npar-1)/2), colspan=npar-1)
fig3.axis([phasemin,phasemax,magdim,magbright])
fig3.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
fig3.set_yticks([9.3,9.4,9.5])
fig3.spines['right'].set_visible(False)
fig3.yaxis.set_ticks_position('left')
plt.plot(phase_lcdat, mag_lcdat, color=red, marker='.', ls='None', ms=2, mew=0) # data
plt.plot(phase_lcmod, mag_lcmod, color='k', lw=1.5) # model
plt.ylabel('Kepler Mag')
plt.xlabel('Orbital Phase')

plt.show()

# ORIGINAL EXAMPLE is below
# This is from http://oceanpython.org/2013/02/25/2d-histogram/
# Plot 2D histogram
#fig1 = plt.figure()
#plt.pcolormesh(xedges, yedges, Hmasked, cmap=cx)
#plt.xlabel('x')
#plt.ylabel('y')
#cbar = plt.colorbar()
#cbar.ax.set_ylabel('Counts')
#plt.show()