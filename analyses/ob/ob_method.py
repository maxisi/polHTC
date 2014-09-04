#! /usr/bin/env python 

import random
import cPickle as pickle
import logging
import sys


'''
OPENS BOX for a given PSR, detector and search method.
'''

##########################################################################################

## PRELUDE

# unpack
process_name, psr, det, run, method, gridsizeSTR = sys.argv

from lib import general as g
g.setuplog('ob_%(det)s%(run)s_%(psr)s_%(method)s' % locals()) # argument added to log filename
from lib import results as res

#note: arguments are taken as strings: need to convert to numerical values.
gridsize = int(gridsizeSTR)
# n is the size of the grid over which psi and inc parameters will be drawn

# setup log
log = logging.getLogger('OB Prep')

log.info('Preparing to open box for PSR %(psr)s %(det)s %(run)s data with %(method)s template.' % locals())
log.info('Using ' + gridsizeSTR + 'values of INC and POL (total ' + str(gridsize**2) + 'searches).')


## LOAD FINEHET DATA
pair = g.Pair(psr, det)
pair.load_finehet(run)
pair.get_sigma()
pair.det.load_vectors(pair.time, filename=psr)


##########################################################################################

## PROCESS

if method in ['GR', 'G4v']:

    src = pair.psr

    ## SETUP
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

    # create n random values for pol and inc
    random.seed(1)
    polsrch = [random.uniform(pol_range[0], pol_range[1]) for n in range(gridsize)]
    incsrch = [random.uniform(inc_range[0], inc_range[1]) for n in range(gridsize)]

    # setup partial results
    h = []
    s = []
    a = []
    pols = []
    incs = []
    
    ## SEARCH
    log.info('Searching.')
    
    for pol in polsrch:
        for inc in incsrch:
            r_ob = pair.search_finehet(methods=[method], pol=pol, inc=inc)
            h += [r_ob['h']]
            s += [r_ob['s']]
            a += [r_ob['a']]
            pols += [pol]
            incs += [inc]
    
    log.info('Taking best values.')
    
    # take set of values that returned highest significance
    best = sorted(zip(s, h, a, pols, incs))[-1]
    
    # pack into final results
    results = {
                's'  : best[0],
                'h'  : best[1],
                'a'  : best[2],
                'pol': best[3],
                'inc': best[4]
                }
else:
    log.info('Searching.')
    results = pair.search_finehet(methods=[method])
            

log.info('Saving results')
try:
    filename = 'ob_' + det + run + '_' + psr + '_' + method
    with open(g.paths['ob'] + filename + '.p', 'wb') as f:
        pickle.dump(results, f)
except:
    log.error('Unable to save OB results', exc_info=True)
