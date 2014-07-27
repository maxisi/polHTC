#! /usr/bin/env python 

import sys
from datetime import datetime
import os
import logging
import numpy as np
import random
import h5py
import socket

import general as g

'''
Sets up analysis given PSR, detector, injection kind, search method, number of
instantiations and number of injections. Builds directory structure and writes DAG file.

Required input arguments:
psr         pulsar name code
det         detector name
run         detector run name (S5, S6)
kind     kind of injection (GR, G4v)
pdif        phase difference between injection components (p, m, 0)
ninst       number of background instantiations
ninj        number of injections
'''

# unpack
process_name, psr, det, run, kind, pdif, ninstSTR, ninjSTR = sys.argv

#note: arguments are taken as strings: need to convert to numerical values.
ninst = int(ninstSTR)
ninj = int(ninjSTR)

# setup log
g.setuplog('master_%(det)s%(run)s_%(psr)s_%(kind)s%(pdif)s_' % locals()) # argument added to log filename
log = logging.getLogger('InjSrch Prep')

log.info('Preparing inj-search for PSR %(psr)s on %(det)s %(run)s data with %(kind)s %(pdif)s injections.' % locals())
log.info('Performing ' + str(ninj) + ' injections on ' + str(ninst) + ' instantiations.')


## ANALYSIS PARAMETERS
frange = [1.0e-7, 1.0e-5] # frequencies for re-heterodynes
hinjrange=[1.0E-27, 1.0E-24] # injection strengths IMP! MIGHT NEED TO BE TUNED!


## STRUCTURE
log.info('Creating file structure.')

pathname = g.analysis_path(det, run, psr, kind, pdif)

for dir in g.localpaths:
    try:
        os.makedirs(pathname + '/' + dir)
        log.debug(dir)
    except OSError:
        log.debug('"globals/%(dir)s" already exists' % locals())
        
log.info('Writing ID record.')

idtext = [
            'INJECTION-SEARCH ANALYSIS\n',
            str(datetime.now()) + '\n',
            'PSR:\t'     + psr,
            'Det:\t'   + det,
            'Run:\t'        + run,
            'Inj:\t'+ kind,
            '\nninst:\t'   + str(ninst),
            'frange:\t' + str(frange),
            '\nninj:\t'     + str(ninj),
            'hinj:\t' + str(hinjrange)   
            ]

with open(pathname + '/id.txt', 'w') as f:
    for line in idtext: f.write(line + '\n')
    f.write('\n###\nREFERENCE: PSR, pulsar name; det, detector name; run, data run (S5, S6); kind, kind of injection (GR, G4v); pdif, phase difference between injection components (p, m, 0); ninst, number of background instantiations; ninj, number of injections.')
    

## CHECK FINEHET DATA EXISTS AND EXTRACT TIME SERIES
psrdet = g.Pair(psr, det)
psrdet.load_finehet(run, check_vectors=True) # also checks vectors for det are existent and healthy.


## RE-HETERODYNES
log.info('Producing frequency list for re-heterodyne')

freq = np.linspace(frange[0], frange[1], ninst)


## SEARCH PARAMETERS
log.info('Preparing search parameters.')

src = psrdet.psr

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
polinj = [random.uniform(pol_range[0], pol_range[1]) for n in np.arange(0, ninj)]
incinj = [random.uniform(inc_range[0], inc_range[1]) for n in np.arange(0, ninj)]


## SAVE ARRAYS
log.info('Saving rehetereodyne, injection and search parameters')

try:
    f = h5py.File(pathname + '/info.hdf5', 'w')
    
    f.create_dataset('srch/freq', data = freq)
    f.create_dataset('srch/pol', data = polsrch)
    f.create_dataset('srch/inc', data = incsrch)
    
    f.create_dataset('inj/h', data = hinj)
    f.create_dataset('inj/pol', data = polinj)
    f.create_dataset('inj/inc', data = incinj)
    
    f.close()
except:
    log.error('FATAL: could not save injsrch info in: '+pathname + '/info.hdf5', exc_info=True)
    sys.exit()


# WRITE SUBMIT FILE

home = os.path.expanduser('~')
project_dir = home + '/polHTC'

# get hostname to determine what server we are on
hostname = socket.gethostname()

if 'ldas' in hostname:
    scratch_dir = '/usr1/max.isi/'
elif 'atlas' in hostname:
    # get result of hostname -s command in bash
    hs = hostname.split('.')[0]
    scratch_dir = '/atlas/user/' + hs + '/max.isi/'
else:
    scratch_dir = project_dir + '/logs/'
    
subfile_lines = [
                'Universe = Vanilla',
                'Executable = ' + project_dir + '/injsrch_process.py',
                'initialdir = ' + project_dir + '/' + pathname + '/',
                'arguments = "%(psr)s, %(det)s, %(run)s, %(kind)s, %(pdif)s, $(Process)"' % locals(),
                'Output = ' + scratch_dir + 'injsrch-$(Process).out',
                'Error = ' + scratch_dir + 'injsrch-$(Process).err',
                'Log = ' + scratch_dir + 'injsrch-$(Process).log',
                'getenv = true',
                'Queue ' + ninstSTR
                ]

with open(g.submit_path(det, run, psr, kind, pdif), 'w') as f:
    for l in subfile_lines: f.write(l + '\n')