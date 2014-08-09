#! /usr/bin/env python 

from globals import general as g

import sys
import numpy as np
import scipy.stats

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


### SETUP

processname, kind = sys.argv

psr = 'J0534+2200'
detname = 'H1'

plots_dir = 'scratch/plots/'

days=1

# create time vector

t = np.array([x+ 630720013 for x in range(0, int(days*g.ss), 60)])
# GPS 630720013 --> UTC 	Jan 01, 2000	00:00:00

# obtain detector and source vectors

src = g.Pulsar(psr)
wx, wy, wz = src.vectors()

det = g.Detector(detname)
det.load_vectors(t)
dx = det.dx; dy = det.dy

### PLOTS

tick_locs = np.linspace(t[0], t[-1], 5)
tick_name = ['00:00', '06:00', '12:00', '18:00', '24:00']

if 'all' in kind:
    plots_dir += 'pols/'
        
    for p in g.pols:
        plt.plot(t, getattr(g, p)(dx,dy,wx,wy,wz), label='$' + g.pol_sym[p] + '$')
        
        if kind=='all':
            plt.title(g.detnames(detname) + ' response to ' + g.pol_names[p] + ' signals from ' + psr)

            plt.ylim(-1,1)
            plt.xlim(t[0], t[-1])
        
            plt.grid(b=True, axis='y')
        
            plt.ylabel('$A_{' + g.pol_sym[p] + '}$')
            plt.xlabel('UTC Time (h)')
        
            plt.xticks(tick_locs, tick_name)
        
            plt.savefig('%(plots_dir)spol_%(p)s_%(detname)s_%(psr)s.pdf' % locals(), bbox_inches='tight')
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
    
        plt.savefig('%(plots_dir)spol_all_%(detname)s_%(psr)s.pdf' % locals(), bbox_inches='tight')
        plt.close()

elif kind in g.pols:
    
    plots_dir += 'pols/'
    
    p = kind
    
    plt.plot(t, getattr(g, p)(dx,dy,wx,wy,wz))
    
    plt.title(g.detnames(detname) + ' response to ' + g.pol_names[p] + ' signals from ' + psr)

    plt.ylim(-1,1)
    plt.xlim(t[0], t[-1])
    
    plt.grid(b=True, axis='y')
    
    plt.ylabel('$A_{' + g.pol_sym[p] + '}$')
    plt.xlabel('UTC Time (h)')
    
    plt.xticks(tick_locs, tick_name)
    
    plt.savefig('%(plots_dir)spol_%(p)s_%(detname)s_%(psr)s.pdf' % locals(), bbox_inches='tight')
    plt.close()

elif kind in ['scalar', 'vector', 'tensor']:
    plots_dir += 'pols/'
        
    for p in g.pol_kinds[kind]:
        plt.plot(t, getattr(g, p)(dx,dy,wx,wy,wz), label='$'+g.pol_sym[p]+'$')
        
    plt.title(g.detnames(detname) + ' response to ' + kind + ' signals from ' + psr)

    plt.ylim(-1,1)
    plt.xlim(t[0], t[-1])
    
    plt.grid(b=True, axis='y')
    
    plt.ylabel('$A_{p}$')
    plt.xlabel('UTC Time (h)')
    
    plt.xticks(tick_locs, tick_name)
    
    plt.legend(numpoints=1, loc='lower left')
    
    plt.savefig('%(plots_dir)spol_%(kind)s_%(detname)s_%(psr)s.pdf' % locals(), bbox_inches='tight')
    plt.close()

elif kind in ['GRp', 'GRm', 'GR0', 'G4vp', 'G4vm', 'G4v0']:
    # Build signal. (copied straight from general)
    
    tmpkind = kind[0:-1]
    pdif = kind[-1]
    
    signal = np.zeros(len(t))+1j*np.zeros(len(t))
    inc = src.param['INC']
    
    for A, a in g.templateinfo[tmpkind].iteritems():
        signal += a(inc, pdif) * (A(dx,dy,wx,wy,wz) + 0j)
    
    # Format plot
    
    plt.title(g.detnames(detname) + ' response to ' + kind + ' signals from ' + psr)


    plt.plot(t, signal.real, 'b', label='Re')
    plt.plot(t, signal.imag, 'r', label='Im')
    plt.plot(t, abs(signal), linestyle='dashed', color='0.5', label='Norm')
    
    plt.ylim(-1,1)
    plt.xlim(t[0], t[-1])
    
    plt.grid(b=True, axis='y')
    
    plt.legend(numpoints=1, loc='lower left')
    
    plt.ylabel('$A_{' + tmpkind + '}$')
    plt.xlabel('UTC Time (h)')
    
    plt.xticks(tick_locs, tick_name)
    
    plt.savefig('%(plots_dir)stem_%(tmpkind)s_%(kind)s_%(psr)s.pdf' % locals(), bbox_inches='tight')
    plt.close()

print 'Plot saved to: ' + plots_dir# + figname