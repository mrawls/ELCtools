from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter, MaxNLocator
'''
This handy program turns your ELC output file party into something more useful.

--> Makes a plot of chi^2 vs. fit parameters from a markovELC / geneticELC run
    (Only parameters in the generation files with the lowest 10k chi2 values are plotted)
--> Spits out the fit parameters with a measure of uncertainty (+ and -) for all
    the parameters in both the generation and ELCparm files

Required file: 'generation.all', created by 'cat generation.1* > generation.all'
Required file: 'ELCparm.all', created by 'cat ELCparm.1* > ELCparm.all'
Required file: 'key.ELCparm' (this is automatically generated by ELC)
Required file: 'gridloop.opt' (you need this to run ELC to start with)

The fit parameters in gridloop.opt are reported in the generation files.
Other parameters of interest which follow from these are reported in the ELCparm files.

The plot will be a 4 x 5 grid. If you are fitting more than 20 parameters in gridloop.opt,
the last ones will be omitted.
'''
# Important filename definitions
gridloopfile = 	'../../RG_ELCmodeling/9246715/chunk4/gridloop.opt'
generationfile = 	'../../RG_ELCmodeling/9246715/chunk4/generation.all'
parmfile = 		'../../RG_ELCmodeling/9246715/chunk4/ELCparm.all'
parmkeyfile = 		'../../RG_ELCmodeling/9246715/chunk4/key.ELCparm'
outfile = 		'../../RG_ELCmodeling/9246715/chunk4/chiplotout.txt'
out = open(outfile, 'w')
gridloop = [line.rstrip('\n') for line in open(gridloopfile)]
nvars = int(gridloop[10]) # reads the number of fit variables from gridloop file

print('Working, please be patient...')
print('Results will be written to {0}'.format(outfile))

# Read in names and limits for fit variables from gridloop file
varnames = []; varlower = []; varupper = []; varlist = []
for i in range(0, nvars):
	varnames.append(gridloop[10+i+1].rstrip())
	varlimits = gridloop[10+nvars+i+1]
	values = varlimits.split()
	varlower.append(float(values[0]))
	varupper.append(float(values[1]))

varnames.append('gamma1'); varnames.append('gamma2') # manually include systemic velocity columns
varlower.append(-100); varlower.append(-100)
varupper.append(100); varupper.append(100)

# Read in chi^2 and parameter values from generation.all file
varlist_gen = np.loadtxt(generationfile, usecols=(range(1,nvars+4)), dtype=np.float64, unpack=True)
chi2_gen = varlist_gen[0]
varlist_gen = np.delete(varlist_gen, 0, 0)

# Read in chi^2 and parameter values from ELCparm.all file
parmkeys = np.loadtxt(parmkeyfile, comments='#', usecols=(1,), dtype={'names':('parmkeys',),'formats':('|S11',)}, unpack=True)
parmkeys = np.delete(parmkeys, 0, 0) # get rid of the index column
nparms = len(parmkeys)
varlist_par = np.loadtxt(parmfile, usecols=(range(1,nparms+2)), dtype=np.float64, unpack=True)
chi2_par = varlist_par[0]
varlist_par = np.delete(varlist_par, 0, 0)
chi2s = chi2_par # higher precision than chi2_gen, but otherwise identical

# Sort parameter arrays by chi2, and only keep the 10,000 values with the lowest chi2 values
sorted_varlist_gen = []; sorted_varlist_par = []
for array in varlist_gen:
	sorted_varlist_gen.append(array[np.argsort(chi2s)][:10000])
for array in varlist_par:
	sorted_varlist_par.append(array[np.argsort(chi2s)][:10000])
sorted_chi2s = chi2s[np.argsort(chi2s)][:10000]
deltachi = np.min(sorted_chi2s) # in a perfect world, deltachi (error bar threshold) = 1.0

# Loop over generation file things
for i in range (1, nvars+3):
	xvalues = sorted_varlist_gen[i-1]

##############################
### PLOTTING STUFF IS HERE ###
	if i <= 5:
		ax = plt.subplot2grid((4,5),(0,i-1))
	elif i > 5 and i <= 10:
		ax = plt.subplot2grid((4,5),(1,i-5-1))
	elif i > 10 and i <= 15:
		ax = plt.subplot2grid((4,5),(2,i-10-1))
	elif i > 15 and i <= 20:
		ax = plt.subplot2grid((4,5),(3,i-15-1))
	#elif i > 20 and i <= 25:
	#	ax = plt.subplot2grid((5,5),(4,i-20-1))
	ax.xaxis.set_major_locator(MaxNLocator(3)) # no more than 2 ticks per axis
	#ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
	ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
	plt.subplots_adjust(wspace=0, hspace=0.4)
	if i % 5 != 1: ax.set_yticklabels(())
	plt.plot(xvalues, sorted_chi2s, marker='o', color='k', mec=None, ms=2, ls='None', markevery=10)
	plt.xlabel(varnames[i-1])
	#plt.text(xmin*1.1, ymin*1.1, varnames[i-1])
	#xmin = np.min(xfilter)
	#xmax = np.max(xfilter)
	#plt.axis([xmin, xmax, ymin, ymax])
	#plt.ylim((plotymin, plotymax))
	#plt.axis([np.min(errorxs), np.max(errorxs), plotymin, plotymax])
### PLOTTING STUFF IS HERE ###
##############################
	
# Print out fit variables from generation file with uncertainties
	bestx = xvalues[0]
	errorx = []
	for idx, (value, chi2) in enumerate(zip(xvalues, chi2s)):
		if (chi2 - chi2s[0]) <= deltachi:
			errorx.append(value)
	errorxplus = np.max(errorx)
	errorxminus = np.min(errorx)
	print(varnames[i-1], '\t =', bestx, '+', errorxplus-bestx, '\ -', bestx-errorxminus, file=out)

# Loop over ELCparm file things	
# Print out fit variables from ELCparm file with uncertainties
for i in range(1, nparms+1):
	xvalues = sorted_varlist_par[i-1]
	bestx = xvalues[0]
	errorx = []
	for idx, (value, chi2) in enumerate(zip(xvalues, chi2s)):
		if (chi2 - chi2s[0]) <= deltachi:
			errorx.append(value)
	errorxplus = np.max(errorx)
	errorxminus = np.min(errorx)
	print(parmkeys[i-1][0], '\t =', bestx, '+', errorxplus-bestx, '\ -', bestx-errorxminus, file=out)

out.close()
print('Chi2 ranges from {0} to {1}'.format(np.min(chi2s), np.max(chi2s)))
print('There are {0} parameters explicitly being fit in gridloop.'.format(nvars))
print('Error bars assume a delta-chi2 threshold of {0}.'.format(deltachi))
print('Here comes a plot...')
plt.show()