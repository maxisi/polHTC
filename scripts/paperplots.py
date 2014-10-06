#! /usr/bin/env python

"""Generates all plots for paper
"""
import os
import sys
import numpy as np
from scipy.interpolate import interp1d
import matplotlib

sys.path.append(os.getcwd())
from lib import results as r

###############################################################################
# CONFIG
p = '/home/misi/results/'

###############################################################################
# METHODS PLOTS (Crab)
# crab_gr = r.Results('H1', 'S5', 'J0534+2200', 'GR')
# crab_gr.load(path=p)
# # s vs hrec
# crab_gr.plot_hs(title=False)
# # p--value
# crab_gr.plot_p('s', title=False, star=10, starsize=12,
#             xlim=(0, 12), legend=False, methods=['G4v'])
# # s/hrec vs hinj example
# crab_gr.plot('s', methods=['GR'], title=False, legend=False, suffix='_ex')
# crab_gr.plot('h', methods=['GR'], title=False, legend=False, suffix='_ex',
#              ylim=(0, .9e-24))

###############################################################################
# CRAB RESULTS PLOTS

# GR
# crab_gr.plot('s', methods=['GR', 'Sid'], band_conf=0, title=False,
#              suffix='_sid')
# crab_gr.plot('s', methods=['GR', 'G4v'], band_conf=0, title=False,
#              suffix='_g4v')
# # G4v
# crab_g4v = r.Results('H1', 'S5', 'J0534+2200', 'G4v')
# crab_g4v.load(path=p)
#
# crab_g4v.plot('s', methods=['GR', 'Sid'], band_conf=0, title=False,
#              suffix='_sid')
# crab_g4v.plot('s', methods=['GR', 'G4v'], band_conf=0, title=False,
#              suffix='_gr')
#
###############################################################################
# ALL PSRs RESULTS PLOTS
detectors = ['H1']#, 'L1', 'H2']
runs = ['S5', 'S6']
injections = ['GR', 'G4v']

for d in detectors:
    for run in runs:
        for injkind in injections:
            results = r.ResultsMP(injkind, det=d, run=run, path=p)
            results.plot_ref(filetype='png')