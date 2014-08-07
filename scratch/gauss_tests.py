#! /usr/bin/env python 

import general as g

import sys
import numpy as np
import scipy.stats

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


### SETUP

processname, x, kind = sys.argv

psr = 'J0534+2200'
det = 'H1'
run = 'S5'

plots_dir = 'scratch/plots/gaussianity/'

# Load data
p = g.Pair(psr, det)
p.load_finehet(run)

t = p.time
d = p.data.real #abs(p.data)



### STAT FUNCTIONS

def normedks(segment):
    # Performs KS test after norming data (see below)
    
    m = np.mean(segment)
    s = np.std(segment)
    
    normed_data = (d-m)/s # KS assumes m=0, s=1; we correct for that.
    
    return scipy.stats.normaltest(segment)


if kind == 'ad':
    # Anderson-Darling test
    # http://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.anderson.html
    
    stat = lambda x : scipy.stats.anderson(x)[0] # return stat value
    yname = 'AD Stat'
    tname = 'Anderson-Darling test'
    
elif kind == 'dak2':
    # D'Agostino's K-squared test
    # http://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.normaltest.html
    
    stat = lambda x : scipy.stats.normaltest(x)[1]
    yname = 'p'
    tname = "D'Agostino's K-squared test"
    
elif kind == 'ks':
    # Kolmogorov-Smirnov test
    
    stat = lambda x : normedks(x)[1]
    yname = 'p'
    tname = 'Kolmogorov-Smirnov test'

else:
    print 'Error: Unknown stat kind "%(kind)s".' % locals()
    sys.exit()


## PLOT

d_gauss = np.random.randn(len(t))

# find number of days which the data spans
ndays = int( np.ceil( (t[-1]-t[0])/g.ss) )

# this will also be the number of bins over which the data will be split
count = np.histogram(t, bins=ndays)[0]

# Note from histogram manual: All but the last (righthand-most) bin are half-open.

# Split data series based on the count. Define initial and final indices for each
# day, i_0 and i_f.
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
        
    i_0 = i_f # update initial index

if x == 'std':
    # x-axis will be the standard deviation

    # Split data series based on the count. Define initial and final indices for each
    # day, i_0 and i_f.
    i_0 = 0
    s = []
    for c in count:

        i_f = i_0 + c  # create final index

        segment = d[i_0:i_f] # pick day-worth of data

        if c != 0:
            s += [np.std(segment)]
        else:
            s += [0]
        
        i_0 = i_f # update initial index    plt.plot(s,a, '+', label='Data')

    plt.plot(s, a, '+', label='Data')
    
    plt.xlabel('$\sigma$')
    figname = '%(kind)s_%(det)s_%(run)s_%(psr)s_std' % locals()
    plt.title('%(tname)s vs $\sigma$ for daily %(det)s %(run)s %(psr)s data' % locals())   

else:
    plt.plot(a, '+', label='Data')

    plt.plot(a_gauss, 'r+', label='Gaussian')

    plt.xlabel('Days')
    figname = '%(kind)s_%(det)s_%(run)s_%(psr)s' % locals()
    plt.title('%(tname)s for daily %(det)s %(run)s %(psr)s data' % locals())
    
plt.ylabel(yname)



plt.legend(numpoints=1)


plt.savefig(plots_dir + figname + '.pdf', bbox_inches='tight')
plt.close()

print 'Plot saved to: ' + plots_dir + figname