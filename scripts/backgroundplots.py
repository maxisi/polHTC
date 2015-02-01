#! /usr/bin/env python

import sys
import os
import numpy as np
import scipy.stats

import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = "stix"

sys.path.append(os.getcwd())
from lib import general as g


###############################################################################
# SETUP

mplparams = {
    'text.usetex': True,  # use LaTeX for all text
    'axes.linewidth': 1,  # set axes linewidths to 0.5
    'axes.grid': False,  # add a grid
    'axes.labelweight': 'normal',
    #'axes.axisbelow': True,
    #'grid.linestyle': '-',  # grid as solid lines not dashed
    #'grid.color': 'gray',  # grid lines grey
    #'grid.linewidth': 0.5,
    'font.family': 'serif',
    'font.size': 24
}
matplotlib.rcParams.update(mplparams)
plots_dir = g.paths['plots']

psr = 'J0534+2200'
ifo = 'H1'
run = 'S5'

# ###############################################################################
# # PLOT APs
#
# # create time vector
# days = 1
# t = np.array([x+ 630720013 for x in range(0, int(days*g.SS), 60)])
# # GPS 630720013 --> UTC 	Jan 01, 2000	00:00:00
#
# # create source and signal objects
# src = g.Pulsar(psr)
# signal = g.Signal.from_names(ifo, psr, t)
#
# ### PLOTS
# tick_locs = np.linspace(t[0], t[-1], 5)
# tick_name = [r'00:00', r'06:00', r'12:00', r'18:00', r'24:00']
#
# for p in g.pols.keys():
#     fig, ax = plt.subplots(1)
#     ax.plot(t, signal(p, pol=src.param['POL']), label='$%s$' % g.pol_sym[p])
#
#     ax.set_xticks(tick_locs)
#     ax.set_xticklabels(tick_name)
#     ax.set_ylim(-1,1)
#     ax.set_xlim(t[0], t[-1])
#
#     ax.grid(b=True, axis='y')
#
#     ax.set_ylabel('$A_{%s}$' % g.pol_sym[p])
#     ax.set_xlabel('UTC Time (h)')
#
#     fontProperties = {'family':'serief','weight': 'light'}
#     a = plt.gca()
#     a.set_xticklabels(tick_name, fontProperties)
#     a.set_yticklabels(a.get_yticks(), fontProperties)
#     matplotlib.rc('font', size=16)
# #    matplotlib.rc('axes', labelsize=20)
#
# #    matplotlib.rc('font', weight='normal')
#
#     figname = '%(plots_dir)spol_%(p)s_%(ifo)s_%(psr)s.pdf' % locals()
#
#     fig.savefig(figname, bbox_inches='tight')
#     plt.close()
#     print '%s AP plot: %s' % (p, figname)
# #    sys.exit(0)
#
#
# ###############################################################################
# # PLOT EXPECTED SIGNALS
#
# for tmpkind in ['GR', 'G4v']:
#     # Build signal
#     sig = signal(tmpkind,pol=src.param['POL'], inc=src.param['INC'])
#
#     plt.plot(t, 2.*sig.real, 'b', label='Re')
#     plt.plot(t, 2.*sig.imag, 'r', label='Im')
#     plt.plot(t, 2.*abs(sig), linestyle='dashed', color='0.5', label='Norm')
#
#     plt.ylim(-1,1)
#     plt.xlim(t[0], t[-1])
#
#     plt.grid(b=True, axis='y')
#
#     plt.legend(numpoints=1, loc='lower right')
#
#     plt.ylabel('$\Lambda$ (' + tmpkind + ')')
#     plt.xlabel('UTC Time (h)')
#
#     plt.xticks(tick_locs, tick_name)
#
#     a = plt.gca()
#     a.set_xticklabels(tick_name, fontProperties)
#     a.set_yticklabels(a.get_yticks(), fontProperties)
#
#     figname = '%(plots_dir)stem_%(tmpkind)s_%(psr)s.pdf' % locals()
#     plt.savefig(figname, bbox_inches='tight')
#     plt.close()
#
#     print '%s signal plot: %s' % (tmpkind, figname)

###############################################################################
# PLOT FINEHET
p = g.Pair(psr, ifo)
p.load_finehet(run)
t = p.time
d = p.data
#
# plt.figure(figsize=(16, 3), dpi=127)
#
# plt.plot(t, d.real)
#
# plt.xlabel('GPS time')
# plt.ylabel('Het. strain (Re)')
# plt.xlim(t[0], t[-1])
#
# figname = '%(plots_dir)sreFinehet_%(ifo)s_%(run)s_%(psr)s' % locals()
#
# plt.savefig(figname, bbox_inches='tight')
# plt.close()
#
# print 'Finehet plot: ' + figname


# ###############################################################################
# # PLOT DAILY STANDARD DEVIATION
# s = p.get_sigma()
#
# fig, ax = plt.subplots(1)
# ax.plot(t, s, '+')
# ax.grid(b=True, axis='both')
# ax.set_xlabel('GPS time', fontsize=20)
# ax.set_ylabel(r'$\sigma$', fontsize=22)
#
# figname = '%(plots_dir)sstd_%(ifo)s_%(run)s_%(psr)s' % locals()
# fig.savefig(figname, bbox_inches='tight')
# plt.close()
#
# print 'Daily STD plot: ' + figname

###############################################################################
# FINEHET HISTOGRAMS

m = np.mean(d.real)
s = np.std(d.real)
gaussian = np.random.normal(m, s, len(d))

for l in [True, False]:
    # produce histograms
    _, _, _ = plt.hist(d.real, 100, color='b', log=l, histtype='step', normed=1,
                       label='Data')
    _, _, _ = plt.hist(gaussian, 100, color='r', log=l, histtype='step', normed=1,
                       label='Gaussian')

    plt.xlabel('Het. strain (Re)')
    plt.ylabel('Normed count')
    plt.legend(numpoints=1)

    figname = '%(plots_dir)shist_%(det)s_%(run)s_%(psr)s' % locals()
    if l:
        figname += '_log'
    plt.savefig(figname + '.pdf', bbox_inches='tight')
    plt.close()

    print 'Finehet histogram: %s' % figname

###############################################################################
# GAUSSIANITY TESTS

## (0) Setup tests:

def normedks(segment):
    # Performs KS test after norming data (see below)
    m = np.mean(segment)
    s = np.std(segment)
    normed_data = (d-m)/s # KS assumes m=0, s=1; we correct for that.
    return scipy.stats.normaltest(segment)

# Anderson-Darling test
# http://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.anderson.html
ad = {
    'stat': lambda x : scipy.stats.anderson(x)[0],
    'yname': 'AD Stat',
    'tname': 'Anderson-Darling test'
}

# Kolmogorov-Smirnov test
ks = {
    'stat': lambda x : normedks(x)[1],
    'yname': 'p',
    'tname': 'Kolmogorov-Smirnov test'
}

## (1) Data

d = p.data.real
# generate gaussian noise
d_gauss = np.random.randn(len(t))
# find number of days which the data spans
ndays = int(np.ceil((t[-1]-t[0])/g.ss))
# this will also be the number of bins over which the data will be split
count = np.histogram(t, bins=ndays)[0]
# Note from histogram manual: All but the last (righthand-most) bins are
# half-open.

for test in [ad, ks]:
    stat = test['stat']
    yname = test['yname']
    tname = test['tname']

    # Split data series based on the count. Define initial and final indices for
    # each day, i_0 and i_f.
    i_0 = 0
    a = []
    a_gauss = []
    for c in count:
        i_f = i_0 + c  # create final index
        segment = d[i_0:i_f] # pick day-worth of data
        segment_gauss = d_gauss[i_0:i_f]
        if c != 0:
            a += [stat(segment)]
            a_gauss += [stat(segment_gauss)]
        else:
            a += [0]
        # update initial index
        i_0 = i_f

    plt.plot(a, '+', label='Data')

    plt.plot(a_gauss, 'r+', label='Gaussian')

    plt.xlabel('Days')
    figname = '%(plots_dir)s%(kind)s_%(det)s_%(run)s_%(psr)s.pdf' % locals()

    plt.ylabel(yname)
    plt.legend(numpoints=1)
    plt.savefig(figname, bbox_inches='tight')
    plt.close()

    print 'Gauss test plot: ' + figname