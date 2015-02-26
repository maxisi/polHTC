#! /usr/bin/env python

import sys
import os
import cPickle as pickle
import logging

import h5py
sys.path.append(os.getcwd())
from polHTC import general as g


"""
Analyzes ONE instantiation of data: starting with finehet seed, reheterodynes,
if necessary injects and searches.

Required input arguments:
psr         pulsar name code
det         detector name
run         detector run name (S5, S6)
injkind     kind of injection (GR, G4v)
n           instantiation index: int in (0, ninst-1)
"""

## PRELUDE

# unpack
process_name, psr, det, run, injkind, nSTR = sys.argv

# note: arguments are taken as strings: need to convert to numerical values.
n = int(nSTR)

pathname = g.analysis_path(det, run, psr, injkind)

# setup log
g.setuplog('p' + nSTR, logpath=pathname+'/logs/')
log = logging.getLogger('InjSrch Process')


## SETUP

log.info('Import search-injection data.')

try:
    with h5py.File(pathname + '/info.hdf5', 'r') as f:
        freq = f['srch/freq'][n]

        hinj = f['inj/h'][n]
        pol = f['inj/pol'][n]
        inc = f['inj/inc'][n]
        ph0 = f['inj/ph0'][n]
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
    inst += hinj * pair.signal(injkind, pol, inc, phi0=ph0)
    # assuming already rescaled by 1/2

# search
results_n = pair.search(data=inst, inc=pair.psr.param['INC'] or None)
# inc=... was added so that if INC is not the default (0), the search is locked

# unpack results
results = {}
for m in g.SEARCHMETHODS:
    results[m] = {
        'h': results_n[m]['h'],
        's': results_n[m]['s'],
        'a': results_n[m]['a']
    }

log.info('Saving results')
pickle.dump(results, open(pathname + '/results/r' + str(n) + '.p', 'wb'))
