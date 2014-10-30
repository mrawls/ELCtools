from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
'''
Meredith Rawls
October 2014

Input a light curve text file, find the eclipses and their durations, and create a 
text file called ELCgap.inp. This is intended for use with Jerry Orosz's ELC code.
It specifies regions BETWEEN eclipses so ELC will only model portions of the
light curve during eclipse AND one eclipse duration before and after.

For example:
	Eclipse from 0 to 10
	No eclipse from 10 to 55
	Eclipse from 55 to 65
	No eclipse from 65 to 110
	Eclipse from 110 to 120
	... etc.
Then ELCgap.inp should read
	20 45
	75 100
	... etc.
	
Credit for the 'contains' function used:
http://stackoverflow.com/questions/3847386/testing-if-a-list-contains-another-list-with-python
'''
# Read file with light curve. 1st column time in days, 2nd column magnitude
filename = 'ELC_lcall.txt'
f1 = open(filename)
times, mag = np.loadtxt(f1, dtype=np.float64, usecols=(0,1), unpack=True)
f1.close()

times = times.tolist() # so we can use the .index() function later

# Knobs you can turn to tweak the eclipse finding algorithm
nsig = 1.5	# if the derivative differs by more than this from the median, flag the time
tsep = 30 	# if two times are separated by more than this, they're separate events

# Function to see if a list is a subset of another list
# Example: a = [1,2,3,4,5] and b = [2,3,4]. Running contains(b,a) will return tuple (1,4).
# You then know that the list a[1:4] is equivalent to the list b.
def contains(small, big):
    for i in xrange(len(big)-len(small)+1):
        for j in xrange(len(small)):
            if big[i+j] != small[j]:
                break
        else:
            return i, i+len(small)
    return False

# Take the derivative of the light curve
deriv1 = np.diff(mag)
sigma = np.std(deriv1)

# Identify time points that are probably in an eclipse, based on the sigma threshold
etimes = []
for idx, (time, value) in enumerate(zip(times[0:-1], deriv1)):
	if np.abs(value) > nsig*sigma and mag[idx] > np.median(mag):
		etimes.append(time)

# Eliminate false positive eclipses
print('This part is slow, please be patient...')
# --> do a rolling window through etimes to make small lists of length 5
# --> check if these small lists are subsets of times using contains
# --> if contains returns something != False, this is a for-reals eclipse
# --> while the contains function doesn't return False, build up a big array
# --> this array, called ereal, is a new and improved (and shorter) version of etimes
ereal = []
ereal.append(etimes[0])
for idx, etime in enumerate(etimes):
	if idx < len(etimes)-4:
		echunk = [etime, etimes[idx+1], etimes[idx+2], etimes[idx+3], etimes[idx+4]]
		if contains(echunk, times) != False: 
			temp0 = times[contains(echunk, times)[0]]
			temp1 = times[contains(echunk, times)[1]]
			if temp0 not in ereal and temp0 >= ereal[-1]: ereal.append(temp0)
			if temp1 not in ereal and temp1 >= ereal[-1]: ereal.append(temp1)
#ereal.append(etimes[-1]) # don't think this is necessary
print('...OK, the slow part is done!')

# Define eclipse start and end times
eclipsecount = 0
eclipsestart = []
eclipseend = []
eclipsestart.append(ereal[0])
for idx, etime in enumerate(ereal):
	if etime - ereal[idx-1] > tsep:
		if etime not in eclipseend:
			eclipsestart.append(etime)
		if ereal[idx-1] not in eclipsestart:
			eclipseend.append(ereal[idx-1])
		eclipsecount += 1
eclipseend.append(ereal[-1])

# Define appropriate time windows that include the entire eclipse + extra on each side
windowstart = []
windowend = []
duration = []
print('We found {0} eclipses'.format(eclipsecount+1))
print('Eclipse durations (some may be weird due to gaps):')
for start, end in zip(eclipsestart, eclipseend):
	duration.append(end - start)
	print(end - start)
med_duration = np.median(duration)
print('Median eclipse duration: ', med_duration)
for start, end in zip(eclipsestart, eclipseend):
	mid = (start + end) / 2
	windowstart.append(mid - 2*med_duration)
	windowend.append(mid + 2*med_duration)

# Sanity check: these should all be the same length
#print(' ')
#print(len(eclipsestart), len(eclipseend))
#print(len(windowstart), len(windowend))
	
# Write data to the ELCgap.inp file
# (need to do this still)
f2 = open('ELCgap.inp', 'w')
print(times[0], windowstart[0], file=f2)
for idx, (start, end) in enumerate(zip(windowstart, windowend)):
	if idx < len(windowstart)-1:
		print(end, windowstart[idx+1], file=f2)
print(windowend[-1], times[-1], file=f2)
f2.close()
print('ELCgap.inp file created')

# Plot the light curve, its derivative, eclipse points, and eclipse window starts/ends
plt.gca().invert_yaxis() # because magnitudes
plt.xlabel('Time (days)')
plt.ylabel('Magnitude')
plt.plot(times, mag, color='0.75', marker='.', linestyle='None', label='Light curve')
plt.plot(times[0:-1], np.abs(deriv1)+np.median(mag)-0.020, color='g', marker='.', linestyle='None', label='LC derivative')
plt.plot(windowstart, np.array([np.median(mag)-0.010]*len(eclipsestart)), 'or', label='Eclipse window start')
plt.plot(windowend, np.array([np.median(mag)-0.011]*len(eclipseend)), 'ob', label='Eclipse window end')
plt.plot(ereal, np.array([np.median(mag)-0.022]*len(ereal)), color='k', marker='.', linestyle='None', label='Points in eclipse')
plt.legend(numpoints=1, loc=4)
plt.show()