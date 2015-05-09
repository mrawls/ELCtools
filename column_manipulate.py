from __future__ import print_function
import numpy as np
'''
Python script to adjust error column of an ELC input file (light curve or RV data).
This is useful when you have run ELC and want to run it again to get estimates
for error bars on your fit parameters.

Jerry recommends scaling the error values so that the fit yields a reduced chi^2 = 1, 
and then estimating 1-sigma errors by looking at how far you have to change the 
parameter for delta-chi^2 = 1. To do this, he suggests dividing the error column of 
each input data set by the square root of the corresponding minimum chi^2 value. 
For example: if the LC min chi^2 is 500, there are 3000 data points, 16 fit parameters, 
and all the error values are 0.001 mag, I would compute new 
error values = 0.001 / sqrt(500/(3000-16)) = 0.0006. Then I run a new model with these 
“fake” error values in order to get error bars. 

Prompts for user to type:
		minimum chi^2, number of points in the file, and number of fit parameters.
Input: 	text file with three columns (time, data, error)
Output: 	new text file with three columns (time, data, newerror)
'''
#infile = 'KIC_9246715_LC_201401.txt'
#infile = 'kplr009246715_RV1better.txt'
infile = 'kplr009246715_RV2better.txt'
time, mag, err = np.loadtxt(infile, comments='#', dtype=np.float64, usecols=(0,1,2), unpack=True)

minchi2 = input("Minimum chi^2 value? ")
npts = input("Number of points in the data file? ")
npar = input("Number of fit parameters? ")

#outfile = 'out1.txt'
#outfile = 'out2.txt'
outfile = 'out3.txt'
f = open(outfile, 'w')

for i in range (0, npts):
	newerr = err[i] / np.sqrt(minchi2/(npts-npar))
	print("%.7f \t %.7f \t %.3f \n" % (time[i], mag[i], newerr), file=f)

f.close()