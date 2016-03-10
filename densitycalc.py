import numpy as np
import astropy.units as u
'''
Given masses and radii of some star, calculate density in solar units x 1000.
Everything is hard-wired, and uncertainties are propagated!
'''

def densitycalc(M=1.0*u.Msun, R=1.0*u.Rsun, Merr=0.*u.Msun, Rerr=0.*u.Rsun):
    density = 3.*M / (4.*np.pi*R**3.)
    err = np.sqrt( (3./(4.*np.pi*R**3.)*Merr)**2. + (-9.*M/(4.*np.pi*R**4.)*Rerr)**2. )
    return density, err

solardensity = densitycalc()[0]

KICs =  [ 929,   395,   100,   578,   703, 979]
Ms =    [1.12, 1.103, 0.857, 1.062, 1.267, 1.39]
Rs =    [7.99, 8.238, 12.73, 11.01, 13.72, 8.84]
Merrs = [0.07, 0.002, 0.060, 0.011, 0.025, 0.03]
Rerrs = [0.39, 0.004, 0.51,  0.04,  0.08,  0.14]

print('For hard-wired values of M and R, the densities in solar units MULTIPLIED BY 1000 are:')

for KIC, M, R, Merr, Rerr in zip(KICs, Ms, Rs, Merrs, Rerrs):
    density, err = densitycalc(M=M*u.Msun, R=R*u.Rsun, Merr=Merr*u.Msun, Rerr=Rerr*u.Rsun)
    print(KIC, (density/solardensity).value*1e3, (err/solardensity).value*1e3)

