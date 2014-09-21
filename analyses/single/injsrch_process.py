#! /usr/bin/env python

import sys
import os
import cPickle as pickle
import logging

import h5py
sys.path.append(os.getcwd())
from lib import general as g


"""
Analyzes ONE instantiation of data: starting with finehet seed, reheterodynes,
if necessary injects and searches.

Required input arguments:
psr         pulsar name code
det         detector name
run         detector run name (S5, S6)
injkind     kind of injection (GR, G4v)
pdif        phase difference between injection components (p, m, 0)
n           instantiation index: int in (0, ninst-1)
"""

## PRELUDE

# unpack
process_name, psr, det, run, injkind, pdif, nSTR = sys.argv

# note: arguments are taken as strings: need to convert to numerical values.
n = int(nSTR)

pathname = g.analysis_path(det, run, psr, injkind, pdif)

# setup log
g.setuplog('p' + nSTR, logpath=pathname+'/logs/')
log = logging.getLogger('InjSrch Process')


## SETUP

log.info('Import search-injection data.')

try:
    with h5py.File(pathname + '/info.hdf5', 'r') as f:
        freq = f['srch/freq'][n]

        hinj = f['inj/h'][n]
        polinj = f['inj/pol'][n]
        incinj = f['inj/inc'][n]
except:
    raise IOError('FATAL: did not find injsrch info in: ' + pathname)

log.info('Create PSR-det pair and load data and std.')

pair = g.Pair(psr, det)
pair.load_finehet(run)
pair.get_sigma()
pair.det.load_vectors(pair.time, filename=psr)

## RE-HETERODYNE

log.info('Reheterodyne.')

inst = g.het(freq, pair.data, pair.time)

## SEARCH

# inject if needed
if hinj != 0:
    log.info('Injecting.')
    inst += (hinj/2.) * pair.signal(injkind, pdif, polinj, incinj)
    # note factor of 1/2, see MP (2.12)

# search
results_n = pair.search(data=inst, pol=pair.psr.param['POL'])

# unpack results
results = {}
for m in g.SEARCHMETHODS:
    results[m] = {
        'h': 2 * results_n[m]['h'],
        's': results_n[m]['s'],
        'a': results_n[m]['a']
    }

log.info('Saving results')
pickle.dump(results, open(pathname + '/results/r' + str(n) + '.p', 'wb'))
