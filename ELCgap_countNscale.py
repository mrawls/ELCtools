from __future__ import print_function
import numpy as np
'''
A not-overly-elegant program in two parts.
There is no good reason for these two parts to be in the same program... sorry.

1.
Read light curve file with time as first column
Read ELCgap.inp file with gap start, end as only two columns
Figure out how many points from the full light curve are actually being used
(this is slow, and there's an option to set it manually)

2.
Write new LC & RV files with scaled errors
(for use with determining ELC parameter error bars)
'''
## FIRST PART ##

# Read file with light curve. 1st column time in days, 2nd column magnitude
lcin = 'KIC_9246715_201408_Patrick.txt'
f1 = open(lcin)
times, mags, errors = np.loadtxt(f1, dtype=np.float64, usecols=(0,1,2), unpack=True)
f1.close()

# Read file with gap info. 1st column gap start, 2nd column gap end (both in days)
filename = 'ELCgap.inp'
f1 = open(filename)
gapstarts, gapends = np.loadtxt(f1, dtype=np.float64, usecols=(0,1), unpack=True)
f1.close()

# Identify how many points are NOT in gaps
gaptimes = []
for time in times:
	for gapstart, gapend in zip(gapstarts, gapends):
		if time > gapstart and time < gapend: # the point is in a gap
			gaptimes.append(time)

totalcount = 0
count = 0

###
### OPTION: manually set these so you don't recount every time you run this program
totalcount = 65896
count = 25523
###
###

if totalcount == 0:
	print('Counting 60,000+ lines is slow, be patient...')
	for time in times:
		totalcount = totalcount + 1
		if time not in gaptimes:
			count = count + 1
else:
	print('The number of points has been set manually in this program:')

print('There are {0} points in the total light curve'.format(totalcount))
print('There are {0} points in the light curve with gaps removed'.format(count))

## SECOND PART ##

# Make a new light curves & RV data files with scaled error bars
npar = 19
rvcount = 23
rv1in = 'KIC_9246715_RV1new.txt'
rv2in = 'KIC_9246715_RV2new.txt'
lcout = 'KIC_9246715_Patrick_scaled2.txt'
rv1out = 'KIC_9246715_RV1_scaled2.txt'
rv2out = 'KIC_9246715_RV2_scaled2.txt'

# these chi2s are for the GENETIC RUN !!!
#lc_chi2 = 173986.67131884
#rv1_chi2 = 4474.7432846
#rv2_chi2 = 4498.8891258

# these chi2s are for the MCMC RUN !!!
lc_chi2 = 168568.96
rv1_chi2 = 2219.8385
rv2_chi2 = 2323.6059

newerrors = []
newrv1errs = []
newrv2errs = []

for error in errors:
	newerrors.append( error / np.sqrt(lc_chi2/(totalcount-npar)))

f1 = open(rv1in)
f2 = open(rv2in)
rv1times, rv1s, rv1errs = np.loadtxt(f1, dtype=np.float64, usecols=(0,1,2), unpack=True)
rv2times, rv2s, rv2errs = np.loadtxt(f2, dtype=np.float64, usecols=(0,1,2), unpack=True)
f1.close()
f2.close()

for rv1err, rv2err in zip(rv1errs, rv2errs):
	newrv1errs.append(rv1err / np.sqrt(rv1_chi2/(rvcount-npar)))
	newrv2errs.append(rv2err / np.sqrt(rv2_chi2/(rvcount-npar)))

f1 = open(lcout, 'w')
for time, mag, newerror in zip(times, mags, newerrors):
	print(time, mag, newerror, file=f1)
f1.close()

f1 = open(rv1out, 'w')
f2 = open(rv2out, 'w')
for rv1time, rv1, rv1err in zip(rv1times, rv1s, rv1errs):
	print(rv1time, rv1, rv1err, file=f1)
for rv2time, rv2, rv2err in zip(rv2times, rv2s, rv2errs):
	print(rv2time, rv2, rv2err, file=f2)
f1.close()
f2.close()

print('New files with scaled errors written to {0}, {1}, {2}'.format(lcout,rv1out,rv2out))