#! /usr/bin/env python

import sys
import os
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.append(os.getcwd())
from polHTC import general as g


### SETUP

processname, kind = sys.argv

psr = 'J0534+2200'
detname = 'H1'

plots_dir = '/Users/maxisi/Desktop/'#g.paths['plots']

days=1

# create time vector

t = np.array([x+ 630720013 for x in range(0, int(days*g.SS), 60)])
# GPS 630720013 --> UTC 	Jan 01, 2000	00:00:00

# create detector and source objects

src = g.Pulsar(psr)
#
# det = g.Detector(ifo)
# det.load_vectors(t)

signal = g.Signal.from_names(detname, psr, t)

### PLOTS

tick_locs = np.linspace(t[0], t[-1], 5)
tick_name = ['00:00', '06:00', '12:00', '18:00', '24:00']

if 'all' in kind:
    #plots_dir += 'pols/'

    for p in g.pols.keys():
        plt.plot(t, signal(p, pol=src.param['POL']), label='$' + g.pol_sym[p] + '$')

        if kind == 'all':
            #plt.title(g.detnames(ifo) + ' response to ' + g.pols[p] + ' signals from ' + psr)

            plt.ylim(-1,1)
            plt.xlim(t[0], t[-1])

            plt.grid(b=True, axis='y')

            plt.ylabel('$A_{' + g.pol_sym[p] + '}$')
            plt.xlabel('UTC Time (h)')

            plt.xticks(tick_locs, tick_name)

            matplotlib.rc('font', size=14)
            matplotlib.rc('axes', labelsize=20)

            plt.savefig('%(plots_dir)spol_%(p)s_%(ifo)s_%(psr)s.pdf' % locals(), bbox_inches='tight')
            plt.close()

    if kind != 'all':
        plt.title(g.detnames(detname) + ' response to signals from ' + psr)

        plt.ylim(-1,1)
        plt.xlim(t[0], t[-1])

        plt.grid(b=True, axis='y')

        plt.ylabel('$A_{p}$')
        plt.xlabel('UTC Time (h)')

        plt.legend(numpoints=1, loc='lower right', ncol=3)

        plt.xticks(tick_locs, tick_name)

        plt.savefig('%(plots_dir)spol_all_%(ifo)s_%(psr)s.pdf' % locals())#, bbox_inches='tight')
        plt.close()

elif kind in g.pols.keys():

    plots_dir += 'pols/'

    p = kind

    plt.plot(t, signal(kind, pol=src.param['POL']))

    plt.title(g.detnames(detname) + ' response to ' + g.pols[p] + ' signals from ' + psr)

    plt.ylim(-1,1)
    plt.xlim(t[0], t[-1])

    plt.grid(b=True, axis='y')

    plt.ylabel('$A_{' + g.pol_sym[p] + '}$')
    plt.xlabel('UTC Time (h)')

    plt.xticks(tick_locs, tick_name)

    plt.savefig('%(plots_dir)spol_%(p)s_%(ifo)s_%(psr)s.pdf' % locals(), bbox_inches='tight')
    plt.close()

elif kind in ['scalar', 'vector', 'tensor']:
    plots_dir += 'pols/'

    for p in g.pol_kinds[kind]:
        plt.plot(t, signal(p, pol=src.param['POL']), label='$'+g.pol_sym[p]+'$')

    plt.title(g.detnames(detname) + ' response to ' + kind + ' signals from ' + psr)

    plt.ylim(-1,1)
    plt.xlim(t[0], t[-1])

    plt.grid(b=True, axis='y')

    plt.ylabel('$A_{p}$')
    plt.xlabel('UTC Time (h)')

    plt.xticks(tick_locs, tick_name)

    plt.legend(numpoints=1, loc='lower left')

    plt.savefig('%(plots_dir)spol_%(kind)s_%(ifo)s_%(psr)s.pdf' % locals(), bbox_inches='tight')
    plt.close()

elif kind in ['GRp', 'GRm', 'GR0', 'G4vp', 'G4vm', 'G4v0']:
    # Parse arguments
    tmpkind = kind[0:-1]
    pdif = kind[-1]

    # Build signal
    sig = signal(tmpkind,pol=src.param['POL'],
                 inc=src.param['INC'])

    # Format plot

    # plt.title(g.detnames(ifo) + ' response to ' + kind + ' signals from ' + psr)

    plt.plot(t, sig.real, 'b', label='Re', lw=2, color="#005851")
    plt.plot(t, sig.imag, 'r', label='Im', lw=2, color="7A303F")
    plt.plot(t, abs(sig), linestyle='dashed', color='0.5', label='Norm', lw=2)

    plt.ylim(-1,1)
    plt.xlim(t[0], t[-1])

    plt.grid(b=True, axis='y')

    plt.legend(numpoints=1, loc='lower right')

    plt.ylabel('$\Lambda$ (' + tmpkind + ')')
    plt.xlabel('UTC Time (h)')

    plt.xticks(tick_locs, tick_name)

    plt.savefig('%(plots_dir)stem_%(tmpkind)s_%(kind)s_%(psr)s.pdf' % locals(), bbox_inches='tight')
    plt.close()

print 'Plot saved to: ' + plots_dir# + figname
