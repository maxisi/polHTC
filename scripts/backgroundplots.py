#! /usr/bin/env python

import sys
import os
import numpy as np

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

###############################################################################
# PLOT APs

# create time vector
days = 1
t = np.array([x+ 630720013 for x in range(0, int(days*g.SS), 60)])
# GPS 630720013 --> UTC 	Jan 01, 2000	00:00:00

# create source and signal objects
src = g.Pulsar(psr)
signal = g.Signal.from_names(ifo, psr, t)

### PLOTS
tick_locs = np.linspace(t[0], t[-1], 5)
tick_name = ['00:00', '06:00', '12:00', '18:00', '24:00']

for p in g.pols.keys():
    plt.plot(t, signal(p, pol=src.param['POL']), label='$%s$' % g.pol_sym[p])

    plt.ylim(-1,1)
    plt.xlim(t[0], t[-1])

    plt.grid(b=True, axis='y')

    plt.ylabel('$A_{%s}$' % g.pol_sym[p])
    plt.xlabel('UTC Time (h)')

    plt.xticks(tick_locs, tick_name)

#    matplotlib.rc('font', size=14)
#    matplotlib.rc('axes', labelsize=20)

    figname = '%(plots_dir)spol_%(p)s_%(ifo)s_%(psr)s.pdf' % locals()

    plt.savefig(figname, bbox_inches='tight')
    plt.close()
    print '%s AP plot: %s' % (p, figname)


###############################################################################
# PLOT EXPECTED SIGNALS

for tmpkind in ['GR', 'G4v']:
    # Build signal
    sig = signal(tmpkind,pol=src.param['POL'], inc=src.param['INC'])

    plt.plot(t, sig.real, 'b', label='Re')
    plt.plot(t, sig.imag, 'r', label='Im')
    plt.plot(t, abs(sig), linestyle='dashed', color='0.5', label='Norm')

    plt.ylim(-1,1)
    plt.xlim(t[0], t[-1])

    plt.grid(b=True, axis='y')

    plt.legend(numpoints=1, loc='lower right')

    plt.ylabel('$\Lambda$ (' + tmpkind + ')')
    plt.xlabel('UTC Time (h)')

    plt.xticks(tick_locs, tick_name)

    figname = '%(plots_dir)stem_%(tmpkind)s_%(kind)s_%(psr)s.pdf' % locals()
    plt.savefig(figname, bbox_inches='tight')
    plt.close()

    print '%s signal plot: %s' % (tmpkind, figname)

###############################################################################
# PLOT FINEHET
p = g.Pair(psr, ifo)
p.load_finehet(run)
t = p.time
d = p.data

plt.figure(figsize=(16, 3), dpi=127)

plt.plot(t, d.real)

plt.xlabel('GPS time')
plt.ylabel('Het. strain (Re)')
plt.xlim(t[0], t[-1])

figname = '%(plots_dir)sreFinehet_%(det)s_%(run)s_%(psr)s' % locals()

plt.savefig(figname, bbox_inches='tight')
plt.close()

print 'Finehet plot: ' + figname


###############################################################################
# PLOT DAILY STANDARD DEVIATION
s = p.get_sigma()

plt.plot(t, s, '+')

plt.xlabel('GPS time')
plt.ylabel(r'$\sigma$')

figname = '%(plots_dir)sstd_%(det)s_%(run)s_%(psr)s' % locals()
plt.savefig(figname, bbox_inches='tight')
plt.close()

print 'Daily STD plot: ' + figname