#! /usr/bin/env python

import random
import sys
import numpy as np
import os

sys.path.append(os.path.expanduser('~') + '/polHTC/')
sys.path.append(os.getcwd())

from polHTC import general as g

"""
Injects signal on Crab data, then reheterodynes and checks the signal is gone.
"""

# PSR and detector data
psr = 'J0534+2200'
det = 'H1'
run = 'S5'

# Reheterodyne info
ninst = 1e4
frange = [1.0e-7, 1.0e-5]
frequencies = np.linspace(frange[0], frange[1], ninst)
freq = frequencies[random.randint(0, ninst)]

# Injection parameters
hinjrange = [1.0E-27, 1.0E-23]
injections = np.linspace(hinjrange[0], hinjrange[1], ninst)
hinj = injections[random.randint(0, ninst)]

# Analyze
pair = g.Pair(psr, det)
pair.load_finehet(run)
pair.get_sigma()
pair.det.load_vectors(pair.time, filename=psr)

# Inject
injected = pair.data + hinj * pair.signal('GR', pair.psr.param['POL'],
                                          pair.psr.param['INC'])
h = [hinj * ap(pair.psr.param['POL'], 0) for _, ap in
     pair.Signal.templates['GR'].iteritems()]
hinj = np.linalg.norm(h)

# Pre-Search
results_1 = pair.search(data=injected, pol=pair.psr.param['POL'],
                        methods=['GR'])

# Reheterodyne
inst = g.het(freq, injected, pair.time)

# Search
results_2 = pair.search(data=inst, pol=pair.psr.param['POL'], methods=['GR'])

# Print
print '=================='
print 'RE-HETERODYNE TEST'
print '----------'
print 'Injected strength: %r' % hinj
print 'Frequency: %r' % freq
print '----------'
print 'BEFORE HETERODYNE'
print 'Recovered strength: %r' % results_1['GR']['h']
print 'Recovered significance: %2f' % results_1['GR']['s']
print '----------'
print 'AFTER HETERODYNE'
print 'Recovered strength: %r' % results_2['GR']['h']
print 'Recovered significance: %2f' % results_2['GR']['s']
print '----------'
print 'SIGNIFICANCE RATIO (2/1): %2f' % (results_2['GR']['s']/results_1['GR']['s'])
print '=================='