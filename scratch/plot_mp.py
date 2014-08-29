#! /usr/bin/env python 

from globs import results as r

import sys


### SETUP

det = 'H1'
run = 'S5'

kind = sys.argv[1]

if len(sys.argv)==2:
    extra_name = ''
else:
    extra_name = sys.argv[2]

results_dir = '/home/misi/results/' + det + run + extra_name + '/' 

# Define result objects and load from disk
mpGR = r.ResultsMP('GR', det=det, run=run)
mpGR.load(results_dir, extra_name=extra_name, verbose=True)

mpG4v = r.ResultsMP('G4v', det=det, run=run)
mpG4v.load(results_dir, extra_name=extra_name)

### PLOT
mpGR.plot(kind, log=True, title=False)
mpG4v.plot(kind, log=True, title=False)
