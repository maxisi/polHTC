#! /usr/bin/env python 

import os
import sys
import h5py
import cPickle as pickle
import logging
import numpy as np

from globals import general as g

'''
Performs a FULL analysis for a given PSR, detector and injection kind and phase.
Loops over instantiations and saves results.
'''

# unpack
process_name, psr, det, run, injkind, pdif, ninstSTR, ninjSTR = sys.argv

#note: arguments are taken as strings: need to convert to numerical values.
ninst = int(ninstSTR)
ninj = int(ninjSTR)

# setup log
g.setuplog('single_%(det)s%(run)s_%(psr)s_%(injkind)s%(pdif)s' % locals()) # argument added to log filename

log = logging.getLogger('InjSrch Prep')

log.info('Preparing inj-search for PSR %(psr)s on %(det)s %(run)s data with %(injkind)s %(pdif)s injections.' % locals())
log.info('Performing ' + str(ninj) + ' injections on ' + str(ninst) + ' instantiations.')


## ANALYSIS PARAMETERS
frange = [1.0e-7, 1.0e-5] # frequencies for re-heterodynes
hinjrange=[1.0E-27, 1.0E-23] # injection strengths IMP! MIGHT NEED TO BE TUNED!


## STRUCTURE
log.info('Creating file structure.')

pathname = g.analysis_path(det, run, psr, injkind, pdif)

for dir in g.localpaths:
    try:
        os.makedirs(pathname + '/' + dir)
        log.debug(dir)
    except OSError:
        log.debug('"globals/%(dir)s" already exists' % locals())
    except:
        log.debug('Unable to create analysis file structure.' % locals(), exc_info=True) 

log.info('Writing ID record.')

idtext = [
            'INJECTION-SEARCH ANALYSIS\n',
            str(datetime.now()) + '\n',
            'PSR:\t'     + psr,
            'Det:\t'     + det,
            'Run:\t'     + run,
            'Inj:\t'     + injkind,
            '\nninst:\t' + str(ninst),
            'frange:\t'  + str(frange),
            '\nninj:\t'  + str(ninj),
            'hinj:\t'    + str(hinjrange)   
            ]

with open(pathname + '/id.txt', 'w') as f:
    for line in idtext: f.write(line + '\n')
    f.write('\n###\nREFERENCE: PSR, pulsar name; det, detector name; run, data run (S5, S6); injkind, kind of injection (GR, G4v); pdif, phase difference between injection components (p, m, 0); ninst, number of background instantiations; ninj, number of injections.')
    

## CHECK FINEHET DATA EXISTS AND EXTRACT TIME SERIES
pair = g.Pair(psr, det)
pair.load_finehet(run)
pair.get_sigma()
pair.det.load_vectors(pair.time, filename=psr)


## GET SEARCH AND INJECTION RANDOM PARAMETERS
def params(src, frange, hinjrange, ninj, ninst, log):
    ## RE-HETERODYNES
    log.info('Producing frequency list for re-heterodyne')

    freq = np.linspace(frange[0], frange[1], ninst)


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
    
    return freq, polsrch, incsrch, hinj, polinj, incinj

freq, polsrch, incsrch, hinj, polinj, incinj = params(pair.psr, frange, hinjrange, ninj, ninst, log)


##########################################################################################

## PRELUDE

for n in np.arange(0, ninst):

    ## RE-HETERODYNE

    log.info('Reheterodyne.')

    inst = g.het(freq, pair.data, pair.time)

    ## SEARCH

    # inject if needed
    if hinj != 0:
        log.info('Injecting.')
        inst += (hinj/2.) * pair.signal(injkind, pdif, polinj, incinj) # note factor of 1/2, see MP (2.12)

    # setup results
    results = g.Results(det, run, psr, injkind, pdif)
    results.hinj = hinj

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
                results.hrec[m] = 2 * (abs(a).sum()) / len(a)
            else:
                results.hrec[m] = 2 * np.linalg.norm(a)
        
            # significance:
            results.srec[m] = np.sqrt(abs(np.dot(a.conj(), np.linalg.solve(cov, a))))
        

    log.info('Saving results')
    results.export()
