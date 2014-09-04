#! /usr/bin/env python 

import sys
import h5py
import cPickle as pickle
import logging
import numpy as np

from lib import general as g

'''
Analyzes ONE instantiation of data: starting with finehet seed, reheterodynes, if
necessary injects and searches.

Required input arguments:
psr         pulsar name code
det         detector name
run         detector run name (S5, S6)
injkind     kind of injection (GR, G4v)
pdif        phase difference between injection components (p, m, 0)
n           instantiation index: int in (0, ninst-1)
'''

## PRELUDE

# unpack
process_name, psr, det, run, injkind, pdif, nSTR = sys.argv

# note: arguments are taken as strings: need to convert to numerical values.
n = int(nSTR)

pathname = g.analysis_path(det, run, psr, injkind, pdif)

# setup log
g.setuplog('p'+ nSTR, logpath=pathname+'/logs/')
log = logging.getLogger('InjSrch Process')


## SETUP

log.info('Import search-injection data.')

try:
    with h5py.File(pathname + '/info.hdf5', 'r') as f:
        freq = f['srch/freq'][n]
        polsrch = f['srch/pol'][n]
        incsrch = f['srch/inc'][n]
    
        hinj = f['inj/h'][n]
        polinj = f['inj/pol'][n]
        incinj = f['inj/inc'][n]
except:
    log.error('FATAL: did not find injsrch info in: ' + pathname, exc_info=True)
    
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
    inst += (hinj/2.) * pair.signal(injkind, pdif, polinj, incinj) # note factor of 1/2, see MP (2.12)

# setup results
results = {}

# search
for m in g.search_methods:
        log.info('Searching: ' + m)
        
        # obtain design matrix and divide by standard deviation
        A = pair.design_matrix(m, polsrch, incsrch) / pair.sigma
        # note that dm will be complex-valued, but the imaginary part is 0.
        # this is useful later when dotting with b
        
        # define data vector
        b = inst / pair.sigma
        
        # perform SVD decomposition (http://web.mit.edu/be.400/www/SVD/Singular_Value_Decomposition.htm)
        U, s, V = np.linalg.svd(A.T, full_matrices=False)
        W = np.diag(1./s)
        # Note that np.linalg.svd returns Vt, not V. in NR notation
        # (http://docs.scipy.org/doc/numpy/reference/generated/numpy.linalg.svd.html)

        
        # define covariance matrix
        cov = np.dot(np.dot(V.T, W**2), V) # See 'Covariance' page in Polarizations tab of LIGO 2013 Notebook
        
        VtW = np.dot(V.T, W)
        Utb = np.dot(U.T, b)
        
        # results:
        a = np.dot(VtW, Utb.T)
        
        # strength:
        if injkind in ['GR', 'G4v']:
            h = 2 * (abs(a).sum()) / len(a)
        else:
            h = 2 * np.linalg.norm(a)
        
        # significance:
        s = np.sqrt(abs(np.dot(a.conj(), np.linalg.solve(cov, a))))
        
        results[m] = {
                        'a' : a,
                        'h' : h,
                        's' : s
                        }

log.info('Saving results')
pickle.dump(results, open(pathname + '/results/r' + str(n) + '.p', 'wb') )
