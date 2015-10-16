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
'fake' error values in order to get error bars. 

Prompts for user to type:   number of fit parameters
Input:                      three text files with three columns each (time, data, error)
                            chi2.all (concatenation of chi.1*; same as in chiplotter.py)
Output:                     three new text files with three columns (time, data, newerror)
'''
chi2file =  '../../RG_ELCmodeling/9246715/demcmc001/chi.all'

infileLC =  '../../RG_ELCmodeling/9246715/demcmc001/KIC_9246715_201408_Patrick_eclipse.txt'
infileRV1 = '../../RG_ELCmodeling/9246715/demcmc001/9246715_rv1_final.txt'
infileRV2 = '../../RG_ELCmodeling/9246715/demcmc001/9246715_rv2_final.txt'

outfileLC =  '../../RG_ELCmodeling/9246715/LC_erroradjust_demcmc.txt'
outfileRV1 = '../../RG_ELCmodeling/9246715/RV1_erroradjust_demcmc.txt'
outfileRV2 = '../../RG_ELCmodeling/9246715/RV2_erroradjust_demcmc.txt'

npar = input("Number of fit parameters? ")

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

timesLC, dataLC, errsLC = np.loadtxt(infileLC, comments='#', dtype=np.float64, usecols=(0,1,2), unpack=True)
timesRV1, dataRV1, errsRV1 = np.loadtxt(infileRV1, comments='#', dtype=np.float64, usecols=(0,1,2), unpack=True)
timesRV2, dataRV2, errsRV2 = np.loadtxt(infileRV2, comments='#', dtype=np.float64, usecols=(0,1,2), unpack=True)

nptsLC = file_len(infileLC)
nptsRV1 = file_len(infileRV1)
nptsRV2 = file_len(infileRV2)

chi2s, chi2sLC, chi2sRV1, chi2sRV2 = np.loadtxt(chi2file, comments='#', dtype=np.float64, usecols=(1,2,4,5), unpack=True)
minchi2LC = np.min(chi2sLC)
minchi2RV1 = np.min(chi2sRV1)
minchi2RV2 = np.min(chi2sRV2)

LC = open(outfileLC, 'w')
RV1 = open(outfileRV1, 'w')
RV2 = open(outfileRV2, 'w')

for (time, data, err) in zip(timesLC, dataLC, errsLC):
    newerr = err / np.sqrt(minchi2LC / (nptsLC - npar))
    print("%.9f \t %.9f \t %.7f" % (time, data, newerr), file=LC)

for (time, data, err) in zip(timesRV1, dataRV1, errsRV1):
    newerr = err / np.sqrt(minchi2RV1 / (nptsRV1 - npar))
    print("%.9f \t %.9f \t %.7f" % (time, data, newerr), file=RV1)
    
for (time, data, err) in zip(timesRV2, dataRV2, errsRV2):
    newerr = err / np.sqrt(minchi2RV2 / (nptsRV2 - npar))
    print("%.9f \t %.9f \t %.7f" % (time, data, newerr), file=RV2)

LC.close()
RV1.close()
RV2.close()

print('New files written: {0}, {1}, {2}'.format(outfileLC, outfileRV1, outfileRV2))