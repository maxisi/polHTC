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
from lib import results as r
from lib import general as g

# PARSE ARGUMENTS
parser = argparse.ArgumentParser()
parser.add_argument("injkind")
parser.add_argument("-d", "--det", default='H1')
parser.add_argument("-r", "--run", default='S5')
parser.add_argument("-t", "--tmp", nargs='+')
parser.add_argument("-fl", "--flim", nargs='+', type=int, default=[30, 1500])
parser.add_argument("-dp", "--datapath", default='/home/misi/results/')
parser.add_argument("-lp", "--linepath", default='/home/misi/Documents/'
                                                 'P1200104-v4/')
parser.add_argument("-pp", "--plotpath", default='tmp/plots/')
parser.add_argument("--prefix", default='')
parser.add_argument("--suffix", default='')
parser.add_argument("--linecolor", default='c')
parser.add_argument("--linewidth", default=1)
parser.add_argument("--linealpha", default=1, type=float)
parser.add_argument("--marker", default='*')
parser.add_argument("--markercolor")
parser.add_argument("--markersize", default=10)
parser.add_argument("--filetype", default='pdf')
parser.add_argument("--legendloc", default='lower right')

args = parser.parse_args()
run = args.run
det = args.det
kind = args.injkind
linepath = args.linepath
flim = args.flim
methods = args.tmp or [kind]

###############################################################################
# LOAD SENSITIVITY CURVES
if run == 'S5':
    if det == 'H1':
        ts = 527.*86400.
        asd = np.loadtxt(linepath + 'lho4k_070318_strain.txt', comments='%')
    elif det == 'H2':
        asd = np.loadtxt(linepath + 'lho2k_070514_strain.txt', comments='%')
        ts = 535.*86400.
    elif det == 'L1':
        ts = 405.*86400.
        asd = np.loadtxt(linepath + 'llo_060604_strain.txt', comments='%')
    else:
        print 'ERROR: invalid detector run combination: %r and %r' % (det, run)
        sys.exit(1)
elif run == 'S6':
    if det == 'H1':
        ts = 238.*86400.
        asd = np.loadtxt(linepath + 'lho4k_15May2010_05hrs17min45secUTC_strain'
                                    '.txt', comments='%')
    elif det == 'L1':
        ts = 225.*86400.
        asd = np.loadtxt(linepath + 'llo4k_31May2010_09hrs05min45secUTC_strain'
                                    '.txt', comments='%')
    else:
        print 'ERROR: invalid detector run combination: %r and %r' % (det, run)
        sys.exit(1)
else:
    print "ERROR: %r injections are not supported." % run
    sys.exit(1)

# get a range of frequencies
fs = np.arange(flim[0], flim[1], 0.1)
# interpolate amplitude spectral densities to same freq range
intf = interp1d(asd[:, 0], asd[:, 1], kind='linear')
# convert to power spectra
sint = np.square(intf(fs))

# scale factor for sensitivity estimate (Dupuis 2005)
sf = 10.8
# get the harmonic mean of the S5 time weighted sensitivities
sens = sf * np.sqrt(sint/ts)  # ! MIGHT BE WRONG

###############################################################################
# LOAD MY DATA
mp = r.ResultsMP(kind, det=det, run=run)
mp.load(path=args.datapath, prefix=args.prefix, suffix=args.suffix)
freq = 2 * np.array([psr.param['FR0'] for psr in mp._psrs]).astype(float)
hmin = mp.getstat('hmin', det_thrsh=.99, det_conf=0)

###############################################################################
# PRINT

for m in methods:
    scale_vector = hmin[m]/ np.sqrt(intf(freq)/ts)
    scale_factor = np.mean(scale_vector)
    print scale_factor