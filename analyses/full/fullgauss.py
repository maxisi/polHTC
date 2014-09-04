#! /usr/bin/env python 

import random
import os
import sys
import h5py
import cPickle as pickle
import logging
import numpy as np

'''
Performs a FULL analysis for a given PSR, detector and injection kind and phase. USES GAUSSIAN NOISE.
Loops over instantiations and saves results.
'''

# unpack
process_name, psr, det, run, injkind, pdif, ninstSTR, ninjSTR = sys.argv

#note: arguments are taken as strings: need to convert to numerical values.
ninst = int(ninstSTR)
ninj = int(ninjSTR)

from lib import general as g
g.setuplog('fullGauss_%(det)s%(run)s_%(psr)s_%(injkind)s%(pdif)s' % locals()) # argument added to log filename
from lib import results as res

# setup log
log = logging.getLogger('InjSrch Prep')

log.info('Preparing inj-search for PSR %(psr)s on %(det)s %(run)s data with %(injkind)s %(pdif)s injections.' % locals())
log.info('Performing ' + str(ninj) + ' injections on ' + str(ninst) + ' instantiations.')


## ANALYSIS PARAMETERS
hinjrange=[1.0E-27, 1.0E-23] # injection strengths IMP! MIGHT NEED TO BE TUNED!

if 'J0534+2200' in psr:
    hinjrange=[1.0E-27, 1.0E-24]

## CHECK FINEHET DATA EXISTS AND EXTRACT TIME SERIES
pair = g.Pair(psr, det)
pair.load_finehet(run)
pair.det.load_vectors(pair.time, filename=psr)

# obtain overall standard deviation of data (will be used to create Gaussian noise)
datastd   = np.std(pair.data)

## GET SEARCH AND INJECTION RANDOM PARAMETERS
def params(src, hinjrange, ninj, ninst, log):

    ## SEARCH PARAMETERS
    log.info('Preparing search parameters.')

    pol_range = [
                    src.param['POL'] - src.param['POL error'],
                    src.param['POL'] + src.param['POL error']
                    ]

    inc_range = [
                    src.param['INC'] - src.param['INC error'],
                    src.param['INC'] + src.param['INC error']
                    ]

    # phi0_range = [0., np.pi/2]

    # create ninst random values for pol and inc
    random.seed(1)
    polsrch = [random.uniform(pol_range[0], pol_range[1]) for n in np.arange(0, ninst)]
    incsrch = [random.uniform(inc_range[0], inc_range[1]) for n in np.arange(0, ninst)]


    ## INJECTION LOCATIONS AND PARAMETERS
    log.info('Preparing injection strengths in range: ' + str(hinjrange))

    injections = np.linspace(hinjrange[0], hinjrange[1], ninj) #list of injection strengths
    injLocations = [int(x) for x in np.linspace(0, ninst, ninj, endpoint=False)] # indices where injections will be placed
    hinj = np.zeros(ninst) # empty vector (most instantiations won't have injections)
    hinj[injLocations] = injections # injection strength index (indicates hinj for each inst)

    log.info('Preparing injection parameters.')

    # create ninj random values for pol and inc
    random.seed(2)

    pols = [random.uniform(pol_range[0], pol_range[1]) for n in np.arange(0, ninj)]
    polinj = np.zeros(ninst) # empty vector (most instantiations won't have injections)
    polinj[injLocations] = pols # injection strength index (indicates hinj for each inst)

    incs = [random.uniform(inc_range[0], inc_range[1]) for n in np.arange(0, ninj)]
    incinj = np.zeros(ninst) # empty vector (most instantiations won't have injections)
    incinj[injLocations] = incs # injection strength index (indicates hinj for each inst)
    
    return polsrch, incsrch, hinj, polinj, incinj

polsrch_lst, incsrch_lst, hinj_lst, polinj_lst, incinj_lst = params(pair.psr, hinjrange, ninj, ninst, log)


##########################################################################################

## PRELUDE

# setup results
results = res.Results(det, run, psr, injkind, pdif, extra_name='gauss')
results.hinj = hinj_lst

# setup random seed
random.seed(3)

for n in np.arange(0, ninst):
    
    # get search parameters
    polsrch = polsrch_lst[n]
    incsrch = incsrch_lst[n]
    
    # get injection parameters
    hinj = hinj_lst[n]
    polinj = polinj_lst[n]
    incinj = incinj_lst[n]
    
    ## Create Gaussian noise

    log.info('Create white noise')

    inst = np.array([random.gauss(0., datastd) for n in range(len(pair.data))])

    ## SEARCH

    # inject if needed
    if hinj != 0:
        log.info('Injecting.')
        inst += (hinj/2.) * pair.signal(injkind, pdif, polinj, incinj) # note factor of 1/2, see MP (2.12)

    # search
    for m in g.search_methods:
            log.info('Searching: ' + m)
        
            # obtain design matrix and divide by standard deviation
            A = pair.design_matrix(m, polsrch, incsrch) / datastd
            # note that dm will be complex-valued, but the imaginary part is 0.
            # this is useful later when dotting with b
        
            # define data vector
            b = inst / datastd
        
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
            if m in ['GR', 'G4v']:
                results.hrec[m] += [2 * (abs(a).sum()) / len(a)]
            else:
                results.hrec[m] += [2 * np.linalg.norm(a)]
        
            # significance:
            results.srec[m] += [np.sqrt(abs(np.dot(a.conj(), np.linalg.solve(cov, a))))]
        

log.info('Saving results')
results.export()
