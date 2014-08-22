#! /usr/bin/env python 

from globals import general as g

import sys
import numpy as np
import scipy.stats

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


### SETUP

processname, x, y = sys.argv

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

def givestat(kind):
    if kind == 'ad':
        # Anderson-Darling test
        # http://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.anderson.html
    
        stat = lambda x : scipy.stats.anderson(x)[0] # return stat value
        tname = 'Anderson-Darling test'
    
    elif kind == 'dak2':
        # D'Agostino's K-squared test
        # http://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.normaltest.html
    
        stat = lambda x : scipy.stats.normaltest(x)[1]
        tname = "D'Agostino's K-squared test"
    
    elif kind == 'ks':
        # Kolmogorov-Smirnov test
    
        stat = lambda x : normedks(x)[1]
        tname = 'Kolmogorov-Smirnov test'
    
    elif kind == 'std':
        stat = lambda x : np.std(x)
        tname = 'Standard deviation'

    else:
        print 'Error: Unknown stat kind "%(kind)s".' % locals()
        sys.exit()
    
    return stat, tname
    
# get functions to correlate
xstat, xname = givestat(x)

ystat, yname = givestat(y)

# find number of days which the data spans
ndays = int( np.ceil( (t[-1]-t[0])/g.ss) )

# this will also be the number of bins over which the data will be split
count = np.histogram(t, bins=ndays)[0]
# Note from histogram manual: All but the last (righthand-most) bins are half-open.

# Split data series based on the count. Define initial and final indices for each
# day, i_0 and i_f.
i_0 = 0
a = []
b = []
for c in count:

    i_f = i_0 + c  # create final index

    segment = d[i_0:i_f] # pick day-worth of data

    if c != 0:
        a += [xstat(segment)]
        b += [ystat(segment)]
    else:
        a += [0]
        b += [0]
        
    i_0 = i_f # update initial index
    
# There is an error that causes the AD test to return an inf value in some cases
# Just replace these entries with 0
a = np.array(a)
a[np.isinf(a)] = 0

b = np.array(b)
b[np.isinf(b)] = 0

# Correlate and print
matrix = np.array([a, b])

print 'Correlation coefficient (normalized covariance matrix) between:' 
print xname + ' and '+ yname
print '\t' + str(np.corrcoef(matrix)[0][1])