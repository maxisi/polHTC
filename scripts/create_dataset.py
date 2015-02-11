#! /usr/bin/env python

import sys
import os
import numpy as np
import argparse
import random

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.append(os.getcwd())
from lib import general as g

# PARSE ARGUMENTS
parser = argparse.ArgumentParser()
parser.add_argument("output_path")
parser.add_argument("-d", "--det", default='H1')
parser.add_argument("-r", "--run", default='S5')
parser.add_argument("-i", "--inject", default=None)
parser.add_argument("--h0", default=1e-24, type=float)
parser.add_argument("-p", "--psr", default="J0534+2200")
parser.add_argument("--real_noise", action="store_true")
parser.add_argument("--start_time", default=630720013.0, type=float)
parser.add_argument("--days", default=1, type=int)
parser.add_argument("--end_time", default=0, type=float)
parser.add_argument("-v", "--verbose", action="store_true")

args = parser.parse_args()
verbose = args.verbose
run = args.run
detname = args.det
psrname = args.psr
injtmp = args.inject

REHETFREQ = 0.85-3
NOISESTD = 1e-23
PHI0 = 0

if args.real_noise or injtmp:
    pair = g.Pair(psrname, detname)
    pair.load_finehet(run)

# produce data
if not args.real_noise:
    # create time vector
    t0 = args.start_time
    tf = args.end_time
    if tf:
        # use provided final time
        if not tf > t0:
            print 'Error: end time must be after start time.'
            sys.exit(1)
        elif verbose:
            print "Time vector: (%.f, %.f)" % (t0, tf)
        time = np.arange(t0, tf, 60)
    else:
        # use number of days requested
        time = np.arange(t0, t0 + args.days * g.SS, 60)
        if verbose:
            print "Time vector: %.1f days starting at GPS %.f." % (args.days,
                                                                   t0)

    if verbose:
        print "Generating Gaussian noise of std %.1e." % NOISESTD
    data = np.array([random.gauss(0., NOISESTD) for n in range(len(time))]) +\
           1j*np.array([random.gauss(0., NOISESTD) for n in range(len(time))])
    description_noise = "Gaussian"
else:
    if verbose:
        print "Using actual %s %s noise for PSR %s." % (run, detname, psrname)
    data = g.het(REHETFREQ, pair.data, pair.time)
    time = pair.time
    description_noise = "Actual"

# inject
if injtmp:
    if verbose:
        print "Injecting %s with h0 = %.1e" % (injtmp, args.h0)
    pair.det.load_vectors(time, filename=psrname)
    inc = pair.psr.param['INC']
    pol = pair.psr.param['POL']
    data += args.h0 * pair.signal(injtmp, pol, inc, PHI0)
    # find effective injection strength
    h = [args.h0 * ap(inc, PHI0) for _, ap in
         pair.Signal.templates[injtmp].iteritems()]
    hinj = np.linalg.norm(h)
    description_inj = "%s (h = %.1e)" % (injtmp, hinj)
else:
    description_inj = "no"

# save
with open(args.output_path, 'w') as f:
    f.write("# %s noise with %s injection.\n" % (description_noise,
                                                 description_inj))
    for (t, re, im) in zip(time, data.real, data.imag):
        f.write("%.f\t%r\t%r\n" % (t, re, im))

print "ASCII data saved: %s" % args.output_path