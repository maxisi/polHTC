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
from polHTC import general as g


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
plots_dir = "/Users/maxisi/Desktop/"

psr = 'J0534+2200'
ifo = 'H1'
run = 'S5'
fontProperties = {'family':'serief','weight': 'light'}

# create time vector
days = 1
t = np.array([x+ 630720013 for x in range(0, int(days*g.SS), 60)])
# GPS 630720013 --> UTC 	Jan 01, 2000	00:00:00
#
# create source and signal objects
src = g.Pulsar(psr)
signal = g.Signal.from_names(ifo, psr, t)
#
### PLOTS
tick_locs = np.linspace(t[0], t[-1], 5)
tick_name = [r'00:00', r'06:00', r'12:00', r'18:00', r'24:00']

for tmpkind in ['GR', 'G4v']:
    # Build signal
    sig = signal(tmpkind,pol=src.param['POL'], inc=src.param['INC'])

    lw=10
    alpha=.6
#
    plt.plot(t, 2.*sig.real, 'b', label='Re', lw=lw, color="#005851", alpha=alpha)
    plt.plot(t, 2.*sig.imag, 'r', label='Im', lw=lw, color="#E53E30", alpha=alpha)
    plt.plot(t, 2.*abs(sig), linestyle='dashed', color='0.5', label='Norm', lw=lw, alpha=alpha)
#
    plt.ylim(-1,1)
    plt.xlim(t[0], t[-1])
#
    plt.grid(b=True, axis='y')
#
    plt.legend(numpoints=1, loc='lower right')
#
    plt.ylabel('$\Lambda$ (' + tmpkind + ')')
    plt.xlabel('UTC Time (h)')
#
    plt.xticks(tick_locs, tick_name)
#
    a = plt.gca()
    a.set_xticklabels(tick_name, fontProperties)
    a.set_yticklabels(a.get_yticks(), fontProperties)
#
    figname = '%(plots_dir)stem_%(tmpkind)s_%(psr)s.pdf' % locals()
    plt.savefig(figname, bbox_inches='tight')
    plt.close()
#
    print '%s signal plot: %s' % (tmpkind, figname)    
