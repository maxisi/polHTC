#! /usr/bin/env python

import sys
from datetime import datetime
import os
import logging
import numpy as np
import random
import h5py

sys.path.append(os.getcwd())
from polHTC import general as g

'''
Sets up analysis given PSR, detector, injection kind, search method, number of
instantiations and number of injections. Builds directory structure and writes
DAG file.

Required input arguments:
psr         pulsar name code
det         detector name
run         detector run name (S5, S6)
kind     kind of injection (GR, G4v)
ninst       number of background instantiations
ninj        number of injections
'''

# unpack
process_name, psr, det, run, kind, ninstSTR, ninjSTR = sys.argv

# note: arguments are taken as strings: need to convert to numerical values.
ninst = int(ninstSTR)
ninj = int(ninjSTR)

# setup log
g.setuplog(
    'master_%(det)s%(run)s_%(psr)s_%(kind)s' % locals())

log = logging.getLogger('InjSrch Prep')

log.info(
    'Preparing inj-search for PSR %(psr)s on %(det)s %(run)s data with'
    '%(kind)s injections.' % locals())
log.info('Performing ' + str(ninj) + ' injections on ' + str(
    ninst) + ' instantiations.')

## CHECK FINEHET DATA EXISTS AND EXTRACT TIME SERIES
pair = g.Pair(psr, det)
pair.load_finehet(run, check_vectors=True)

## ANALYSIS PARAMETERS
# frequencies for re-heterodynes
frange = [1.0e-7, 0.85e-3]
# injection strengths proportional to overall noise magnitude
hinj_magnitude = np.ceil(np.log10(abs(max(pair.data)))) - 3
hinjrange = [1.0E-27, 10**hinj_magnitude]

## STRUCTURE
log.info('Creating file structure.')

pathname = g.analysis_path(det, run, psr, kind)

for d in g.localpaths:
    try:
        os.makedirs(pathname + '/' + d)
        log.debug(d)
    except OSError:
        log.debug('"%(pathname)s/%(d)s" already exists' % locals())

log.info('Writing ID record.')

idtext = [
    'INJECTION-SEARCH ANALYSIS\n',
    str(datetime.now()) + '\n',
    'PSR:\t' + psr,
    'Det:\t' + det,
    'Run:\t' + run,
    'Inj:\t' + kind,
    '\nninst:\t' + str(ninst),
    'frange:\t' + str(frange),
    '\nninj:\t' + str(ninj),
    'hinj:\t' + str(hinjrange)
]

with open(pathname + '/id.txt', 'w') as f:
    for line in idtext:
        f.write(line + '\n')
    f.write(
        '\n###\nREFERENCE: PSR, pulsar name; det, detector name; run, data run'
        ' (S5, S6); kind, kind of injection (GR, G4v);'
        ' between injection components (p, m, 0); ninst, number of background'
        ' instantiations; ninj, number of injections.')


## RE-HETERODYNES
log.info('Producing frequency list for re-heterodyne')

freq = np.linspace(frange[0], frange[1], ninst)


## SEARCH PARAMETERS
log.info('Preparing search parameters.')

src = pair.psr

pol_range = [
    src.param['POL'] - src.param['POL error'],
    src.param['POL'] + src.param['POL error']
]

inc_range = [
    src.param['INC'] - src.param['INC error'],
    src.param['INC'] + src.param['INC error']
]


## INJECTION LOCATIONS AND PARAMETERS
log.info('Preparing injection strengths in range: ' + str(hinjrange))

#list of injection strengths:
injections = np.linspace(hinjrange[0], hinjrange[1], ninj)
# indices where injections will be placed:
injLocations = [int(x) for x in np.linspace(0, ninst, ninj, endpoint=False)]
# empty vector (most instantiations won't have injections):
hinj = np.zeros(ninst)
# injection strength index (indicates hinj for each inst):
hinj[injLocations] = injections

log.info('Preparing injection parameters.')

# create ninj random values for pol and inc
random.seed(2)

pols = [random.uniform(pol_range[0], pol_range[1]) for n in range(ninj)]
# empty vector (most instantiations won't have injections):
pol = np.zeros(ninst)
# injection strength index (indicates hinj for each inst):
pol[injLocations] = pols

incs = [random.uniform(inc_range[0], inc_range[1]) for ni in range(ninj)]
# empty vector (most instantiations won't have injections):
inc = np.zeros(ninst)
# injection strength index (indicates hinj for each inst):
inc[injLocations] = incs

ph0s = [random.uniform(0, 2*np.pi) for np0 in range(ninj)]
# empty vector (most instantiations won't have injections):
ph0 = np.zeros(ninst)
# injection strength index (indicates hinj_lst for each inst):
ph0[injLocations] = ph0s


## SAVE ARRAYS
log.info('Saving rehetereodyne, injection and search parameters')

f = h5py.File(pathname + '/info.hdf5', 'w')

srch_grp = f.create_group('srch')

srch_grp.create_dataset('freq', data=freq)

inj_grp = f.create_group('inj')

inj_grp.create_dataset('h', data=hinj)
inj_grp.create_dataset('pol', data=pol)
inj_grp.create_dataset('inc', data=inc)
inj_grp.create_dataset('ph0', data=ph0)

f.close()

# WRITE SUBMIT FILE

log.info('Write submit file.')

try:
    home = os.path.expanduser('~')
    project_dir = home + '/polHTC'

    analysis_folder = '/analyses/single'

    # get hostname to determine what server we are on
    cluster = g.Cluster()

    # define scratch space
    scratch_dir = cluster.scratch_dir

    subfile_lines = [
        'Universe = Vanilla',
        'Executable = ' + project_dir + analysis_folder +
        '/injsrch_process.py',

        'initialdir = ' + project_dir + '/',
        'arguments = "%(psr)s %(det)s %(run)s %(kind)s $(Process)"'
        % locals(),
        'Output = ' + scratch_dir + 'injsrch-$(Process).out',
        'Error = ' + scratch_dir + 'injsrch-$(Process).err',
        'Log = ' + scratch_dir + 'injsrch-$(Process).log',
        'getenv = true',
        'Queue ' + ninstSTR
    ]

    with open(g.submit_path(det, run, psr, kind), 'w') as f:
        for l in subfile_lines:
            f.write(l + '\n')
except IOError:
    log.error('Unable to write submit file!', exc_info=True)
