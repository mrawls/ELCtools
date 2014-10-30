ELCtools
========
Various python programs for interacting with Jerry Orosz's ELC code.

The mid-2014 version of ELC itself is included here, but the documentation is incomplete and this repo is not primarily concerned with compiling / running / understanding all the tweakable settings in ELC. It is in FORTRAN. All I will say is, you can try to compile geneticELC and/or markovELC with a command like this:

> ifort -O2 -f77rtl -132 -mcmodel=medium -shared-intel markovELC.for -o markovELC

But if you care deeply about ELC itself, you should be in touch with Jerry (jorosz@mail.sdsu.edu).

Tools in this repo include:
  - ELCgapmaker.py
  - ELCplotter_new.py
  - 2Dhistgrid.py
  - ... and more to come
