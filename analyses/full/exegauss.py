#! /usr/bin/env python

import random
import os
import sys
import h5py
import cPickle as pickle
import logging
import numpy as np

'''
Performs a FULL analysis for a given PSR, detector and injection kind and phase
USES GAUSSIAN NOISE. Loops over instantiations and saves results.
'''

# unpack
process_name, psr, det, run, injkind, ninstSTR, ninjSTR = sys.argv

#note: arguments are taken as strings: need to convert to numerical values.
ninst = int(ninstSTR)
ninj = int(ninjSTR)

sys.path.append(os.path.expanduser('~') + '/polHTC/')
sys.path.append(os.getcwd())

from lib import general as g
g.setuplog('fullGauss_%(det)s%(run)s_%(psr)s_%(injkind)s' % locals())
from lib import results as res

# setup log
log = logging.getLogger('InjSrch Prep')

log.info('Preparing inj-search for PSR %(psr)s on %(det)s %(run)s data with'
         ' %(injkind)s injections.' % locals())
log.info('Performing ' + str(ninj) + ' injections on ' + str(ninst) +
         ' instantiations.')

## CHECK FINEHET DATA EXISTS AND EXTRACT TIME SERIES
pair = g.Pair(psr, det)
pair.load_finehet(run)
pair.det.load_vectors(pair.time, filename=psr)

# obtain overall standard deviation of data
# (will be used to create Gaussian noise)
datastd = np.std(pair.data)

## ANALYSIS PARAMETERS
# frequencies for re-heterodynes
frange = [1.0e-7, 0.85-3]
# injection strengths proportional to overall noise magnitude
hinj_magnitude = int(np.round(np.log10(abs(datastd))) - 1)
hinjrange = [1.0E-27, 10**hinj_magnitude]

## GET SEARCH AND INJECTION RANDOM PARAMETERS

## RE-HETERODYNES
log.info('Producing frequency list for re-heterodyne')

freq_lst = np.linspace(frange[0], frange[1], ninst)

log.info('Preparing search parameters.')

pol_range = [
    pair.psr.param['POL'] - pair.psr.param['POL error'],
    pair.psr.param['POL'] + pair.psr.param['POL error']
]

inc_range = [
    pair.psr.param['INC'] - pair.psr.param['INC error'],
    pair.psr.param['INC'] + pair.psr.param['INC error']
    ]

# phi0_range = [0., np.pi/2]

## INJECTION LOCATIONS AND PARAMETERS
log.info('Preparing injection strengths in range: ' + str(hinjrange))

# list of injection strengths:
injections = np.linspace(hinjrange[0], hinjrange[1], ninj)
# indices where injections will be placed:
injLocations = [int(x) for x in np.linspace(0, ninst, ninj, endpoint=False)]
# empty vector (most instantiations won't have injections):
hinj_lst = np.zeros(ninst)
# injection strength index (indicates hinj_lst for each inst):
hinj_lst[injLocations] = injections

# create ninj random values for pol and inc:
random.seed(2)

pols = [random.uniform(pol_range[0], pol_range[1]) for ip in range(ninj)]
# empty vector (most instantiations won't have injections):
pol_lst = np.zeros(ninst)
# injection strength index (indicates hinj_lst for each inst):
pol_lst[injLocations] = pols

incs = [random.uniform(inc_range[0], inc_range[1]) for ii in range(ninj)]
# empty vector (most instantiations won't have injections):
inc_lst = np.zeros(ninst)
# injection strength index (indicates hinj_lst for each inst):
inc_lst[injLocations] = incs

ph0s = [random.uniform(0, 2*np.pi) for ip0 in range(ninj)]
# empty vector (most instantiations won't have injections):
ph0_lst = np.zeros(ninst)
# injection strength index (indicates hinj_lst for each inst):
ph0_lst[injLocations] = ph0s

###############################################################################
# setup results
results = res.Results(det, run, psr, injkind, prefix='gauss_')
results.hinj = []


for n in range(ninst):

    freq = freq_lst[n]

    hinj = hinj_lst[n]
    pol = pol_lst[n]
    inc = inc_lst[n]
    ph0 = ph0_lst[n]

    ## Create Gaussian noise
    log.info('Create white noise.')
    inst = np.array([random.gauss(0., datastd) for n in range(len(pair.data))])

    ## SEARCH
    # inject if needed
    if hinj != 0:
        log.info('Injecting.')
        inst += hinj * pair.signal(injkind, pol, inc, ph0)
        # assuming already rescaled by 1/2
        h = [hinj * ap(inc, ph0) for _, ap in
             pair.Signal.templates[injkind].iteritems()]
        hinj = np.linalg.norm(h)
    # search
    results_n = pair.search(data=inst, pol=pair.psr.param['POL'], std=datastd)
    # unpack results
    results.hinj.append(hinj)
    for m in g.SEARCHMETHODS:
        results.hrec[m].append(results_n[m]['h'])
        results.srec[m].append(results_n[m]['s'])
        results.arec[m].append(results_n[m]['a'])

log.info('Saving results')
results.export()
