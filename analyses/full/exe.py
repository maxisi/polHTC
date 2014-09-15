#! /usr/bin/env python

import random
import sys
import logging
import numpy as np
import os

"""
Performs a FULL analysis for a given PSR, detector and injection kind & phase.
Loops over instantiations and saves results.
"""

# unpack
process_name, psr, det, run, injkind, pdif, ninstSTR, ninjSTR = sys.argv

#note: arguments are taken as strings: need to convert to numerical values.
ninst = int(ninstSTR)
ninj = int(ninjSTR)

sys.path.append(os.path.expanduser('~') + '/polHTC/')
sys.path.append(os.getcwd())

from lib import general as g
g.setuplog('full_%(det)s%(run)s_%(psr)s_%(injkind)s%(pdif)s' % locals())
from lib import results as res

# setup log
log = logging.getLogger('InjSrch Prep')

log.info('Preparing inj-search for PSR %(psr)s on %(det)s %(run)s data with'
         '%(injkind)s %(pdif)s injections.' % locals())
log.info('Performing ' + str(ninj) + ' injections on ' + str(ninst) +
         ' instantiations.')


## ANALYSIS PARAMETERS
frange = [1.0e-7, 1.0e-5]  # frequencies for re-heterodynes
hinjrange=[1.0E-27, 1.0E-23]  # injection strengths IMP: MIGHT NEED TO BE TUNED

## CHECK FINEHET DATA EXISTS AND EXTRACT TIME SERIES
pair = g.Pair(psr, det)
pair.load_finehet(run)
pair.get_sigma()
pair.det.load_vectors(pair.time, filename=psr)


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

pols = [random.uniform(pol_range[0], pol_range[1])
        for i in np.arange(0, ninj)]
# empty vector (most instantiations won't have injections):
polinj_lst = np.zeros(ninst)
# injection strength index (indicates hinj_lst for each inst):
polinj_lst[injLocations] = pols

incs = [random.uniform(inc_range[0], inc_range[1])
        for i in np.arange(0, ninj)]
# empty vector (most instantiations won't have injections):
incinj_lst = np.zeros(ninst)
# injection strength index (indicates hinj_lst for each inst):
incinj_lst[injLocations] = incs


###############################################################################

## PRELUDE

# setup results
results = res.Results(det, run, psr, injkind, pdif)
results.hinj = hinj_lst


for n in np.arange(0, ninst):

    freq = freq_lst[n]

    hinj = hinj_lst[n]
    polinj = polinj_lst[n]
    incinj = incinj_lst[n]

    ## RE-HETERODYNE

    log.info('Reheterodyne.')

    inst = g.het(freq, pair.data, pair.time)

    ## SEARCH

    # inject if needed
    if hinj != 0:
        log.info('Injecting.')
        inst += (hinj/2.) * pair.signal(injkind, pdif, polinj, incinj)
        # note factor of 1/2, see MP (2.12)

    results_n = pair.search(data=inst, pol=pair.psr.param['POL'])

    # unpack results
    for m in g.SEARCHMETHODS:
        results.hrec[m] += [results_n[m]['h']]
        results.srec[m] += [results_n[m]['s']]
        results.arec[m] += [results_n[m]['a']]


log.info('Saving results')
results.export()
