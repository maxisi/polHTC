#! /usr/bin/env python

"""Generates all plots for paper
"""
import os
import sys
sys.path.append(os.getcwd())
from lib import results as r

###############################################################################
# LOAD RESULTS
p = '/home/misi/results/'
detectors = ['H1']#, 'L1', 'H2']
runs = ['S5', 'S6']
injections = ['GR', 'G4v']

for d in detectors:
    for run in runs:
        for injkind in injections:
            # load results for all pulsars
            results = r.ResultsMP(injkind, det=d, run=run, path=p)
            # pick results for Crab pulsar
            crab = results._results['J0534+2200']
            #if d == 'H1' and injkind == 'GR':
                # plot Crab:
                #crab.plot_hs(title=False)
                #crab.plot_p('s', title=False, star=10, starsize=12,
                #            xlim=(0,12), legend=False, methods=['G4v'])
            crab.plot('s', title=False, suffix='ex')