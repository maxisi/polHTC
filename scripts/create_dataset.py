#! /usr/bin/env python

import sys
import os
import numpy as np
import argparse
from scipy.interpolate import interp1d

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
        print "Time vector: %.1f days starting at GPS %.f" % (args.days, t0)
tf = time[-1]

# produce data
if not args.real_noise:
