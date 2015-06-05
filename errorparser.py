from __future__ import print_function
import numpy as np
'''
Takes a set of chiplotter.py-generated text files (from several ELC runs)
and reports "formal" parameter values with error bars.

Assumes 0th file is the "main ELC run" with the correct parameter values
Assumes 1st file is the secondary "error ELC run" to be used for error bar values
Assumes 2nd - Nth files are for each light curve chunk or segment
'''
infilelist = '../../RG_ELCmodeling/9246715/chi2outfilelist.txt'
infiles = [line.rstrip('\n') for line in open(infilelist)]
parnames = np.loadtxt(infiles[0], usecols=(0,), dtype={'names':('parnames',),'formats':('|S6',)}, unpack=True)

# Read stuff in from the chiplotter outfiles
valuelist = []; pluerrlist = []; minerrlist = []
for infile in infiles:
	value, pluerr, minerr = np.loadtxt(infile, usecols=(2,4,7), unpack=True)
	valuelist.append(value)
	pluerrlist.append(pluerr)
	minerrlist.append(minerr)

print('#here are global errors to choose from - values are from main run')
print('#diff = difference in parameter value between main run and error run')
print('#err1 = error value calculated for main run (larger of + or - error)')
print('#err2 = error value calculated for error run (larger of + or - error)')
prnit('#')

for idx, par in enumerate(parnames):
	if np.abs(valuelist[1][idx]-valuelist[0][idx]) != np.max([pluerrlist[1][idx], minerrlist[1][idx]]):
		# only prints values that are meaningful, i.e., with nonzero errors
		print(par[0], valuelist[0][idx], '\t diff', np.abs(valuelist[1][idx]-valuelist[0][idx]), 
			'\t err1', np.max([pluerrlist[1][idx], minerrlist[1][idx]]), 
			'\t err2', np.max([pluerrlist[0][idx], minerrlist[0][idx]]))

print(' ')
print('#errors based on chunk runs are next - values are rms')
print('#plus and minus are medians +/- errors added in quadrature')
print('#rmse is the rms error (based on the spread of the parameter values)')
print('#')

for idx, par in enumerate(parnames):
	parsave = []; plusave = []; minsave = []
	for chunk in range(2,len(valuelist)):
		parsave.append(valuelist[chunk][idx])
		plusave.append(pluerrlist[chunk][idx])
		minsave.append(minerrlist[chunk][idx])
	# adding errors in quadrature is fun
	# I could probably do this in one line if I cared enough
	pluerrquad = 0; minerrquad = 0
	for plu, min in zip(plusave, minsave):
		pluerrquad += np.power(plu,2)
		minerrquad += np.power(min,2)
	pluerrquad = np.sqrt(pluerrquad)
	minerrquad = np.sqrt(minerrquad)
	rms = np.sqrt(np.mean(np.power(parsave,2)))
	rmserr = np.sqrt(np.mean(valuelist[0][idx] - parsave)**2)
	if pluerrquad != 0 and minerrquad != 0:
	# only prints values that are meaningful, i.e., with nonzero errors
		print(par[0], rms, #np.median(parsave),
			'\t plus', np.median(pluerrquad),
			'\t minus', np.median(minerrquad),
			'\t rmse', rmserr)
