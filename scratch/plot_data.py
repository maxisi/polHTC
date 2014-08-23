#! /usr/bin/env python 

from globs import general as g

import sys
import numpy as np
import scipy.stats

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


### SETUP

processname, type = sys.argv

psr = 'J0534+2200'
det = 'H1'
run = 'S5'

plots_dir = 'scratch/plots/'

# Load data
p = g.Pair(psr, det)
p.load_finehet(run)

t = p.time
d = p.data


### PLOTS

if type == 'data':
    # PLOT FINEHET DATA (RE)
    
    plt.figure(figsize=(16,3), dpi=127)
    # Plot
    plt.plot(t, d.real)

    #plt.title('%(det)s %(run)s data heterodyned for PSR %(psr)s' % locals() )

    plt.xlabel('GPS time')
    plt.ylabel('Het. strain (Re)')
    plt.xlim(t[0], t[-1])

    figname = 'reFinehet_%(det)s_%(run)s_%(psr)s' % locals()

    plt.savefig(plots_dir + figname, bbox_inches='tight')
    plt.close()

    print 'Finehet plot saved to: ' + plots_dir + figname

elif type in ['hist', 'histlog']:
    # HISTOGRAM FINEHET
    
    figname = 'hist_%(det)s_%(run)s_%(psr)s' % locals()
    
    if type == 'histlog':
        l = True
        figname += '_log'
    else:
        l = False
    
    _, _, _ = plt.hist(d.real, 100, color='b', log=l, histtype='step', normed=1, label='Data')
    
    # Fit gaussian
    m = np.mean(d.real)
    s = np.std(d.real)

    gaussian = np.random.normal(m, s, len(d))
    
    _, _, _ = plt.hist(gaussian, 100, color='r', log=l, histtype='step', normed=1, label='Gaussian')
    
    plt.title('Histogram of %(det)s %(run)s data heterodyned for PSR %(psr)s' % locals() )

    plt.xlabel('Real part of het. data')
    plt.ylabel('Normed count')

    plt.legend(numpoints=1)

    plt.savefig(plots_dir + figname + '.pdf', bbox_inches='tight')
    plt.close()

elif type == 'dailystd':
    # PLOT DAILY STD

    s = p.get_sigma()

    plt.plot(t, s, '+')
    
    plt.title('Daily std. dev. of %(det)s %(run)s data heterodyned for PSR %(psr)s' %locals() )
    
    plt.xlabel('GPS time')
    plt.ylabel('$\sigma$')

    figname = 'std_%(det)s_%(run)s_%(psr)s' % locals()
    plt.savefig(plots_dir + figname, bbox_inches='tight')
    plt.close()


print 'Plot saved to: ' + plots_dir + figname
