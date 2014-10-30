ELCtools
========
Various python programs for interacting with Jerry Orosz's ELC code.

The mid-2014 version of ELC itself is included here, but the documentation is incomplete and this repo is not primarily concerned with compiling / running / understanding all the tweakable settings in ELC. All I will say is, you can try to compile geneticELC and/or markovELC with a command like this:

> ifort -O2 -f77rtl -132 -mcmodel=medium -shared-intel markovELC.for -o markovELC

But if you care deeply about ELC itself, you should be in touch with Jerry (jorosz@mail.sdsu.edu).

Tools in this repo include:
  - <b>ELCgapmaker.py</b>: creates the ELCgap.inp file by identifying "gap" regions outside of eclipses. This can be used to specify which time chunks of a light curve to ignore when running ELC.
  - <b>ELCplotter_new.py</b>: plots light curve, RV curve, and best fits thereof after running geneticELC or markovELC. Includes "zoom panels" for both the primary and secondary eclipses.
  - <b>2Dhistgrid.py</b>: plots a whole bunch of 2D histograms of fit variables vs. each other so you can see the parameter space explored by either geneticELC or markovELC.
  - ... and more to come

You might also be interested in my kepler-makelc repo, which contains:
  - <b>ELClcprep.py</b>: turns output from kepler-makelc/makelc.py into an ELC-friendly light curve text file
  - (among other great things)
