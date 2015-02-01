#! /usr/bin/env python

"""
Compares results for a given detector, run and injection kind with their
Gaussian counterpart and the corresponding rescaled detector ASD curve.
"""

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
parser.add_argument("--linecolor", default='k')
parser.add_argument("--linewidth", default=1)
parser.add_argument("--linealpha", default=.3, type=float)
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
hmin = mp.getstat('hmin', det_thrsh=.999, det_conf=.95)

mp_gauss = r.ResultsMP(kind, det=det, run=run)
mp_gauss.load(path=args.datapath, prefix='gauss_', suffix=args.suffix)
freq_gauss = 2 * np.array([psr.param['FR0']
                           for psr in mp_gauss._psrs]).astype(float)
hmin_gauss = mp_gauss.getstat('hmin', det_thrsh=.999, det_conf=.95)

###############################################################################
# PLOT
mplparams = {
    'text.usetex': True,  # use LaTeX for all text
    'axes.linewidth': 1,  # set axes linewidths to 0.5
    'axes.grid': True,  # add a grid
    'axes.axisbelow': True,
    'grid.linestyle': '-',  # grid as solid lines not dashed
    'grid.color': 'gray',  # grid lines grey
    'grid.linewidth': 0.5,
    'font.family': 'serif',
    'font.size': 18
}
matplotlib.rcParams.update(mplparams)

fig = plt.figure(figsize=(14, 10), dpi=300)

plt.loglog(fs, sens, color=args.linecolor, alpha=args.linealpha,
           linewidth=args.linewidth, label='Expected: %s %s' % (det, run))

plt.hold(True)

for m in methods:
    plt.loglog(freq_gauss, hmin_gauss[m], '.k',
               label="%s template (Gaussian)" % m, markersize=args.markersize)
               # markerfacecolor='none',
               # markeredgecolor=args.markercolor or g.plotcolor[m])
    plt.loglog(freq, hmin[m], args.marker, label="%s template" % m,
               markersize=args.markersize,
               color=args.markercolor or g.plotcolor[m])

plt.xlim(flim[0], flim[1])
plt.ylim(1e-27, 1e-23)
plt.xlabel(r'Gravitational-wave Frequency (Hz)', fontsize=20, fontweight=100)
plt.ylabel(r'Strain Sensitivity', fontsize=20, fontweight=100)

plt.legend(numpoints=1, loc=args.legendloc)
plt.setp(plt.gca().get_legend().get_texts(), fontsize='12')

#fig.subplots_adjust(left=0.18, bottom=0.15)

fig.savefig(args.plotpath + args.prefix + 'gausscomp_senscurve_%s%s_%s.%s'
            % (det, run, kind, args.filetype), bbox_inches='tight')

print "Plot saved: " + args.plotpath + args.prefix +\
      'gausscomp_senscurve_%s%s_%s.%s' % (det, run, kind, args.filetype)
