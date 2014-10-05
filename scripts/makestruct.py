#! /usr/bin/env python

'''
Creates basic file structure for polarization analysis (if not already existent).
'''

import os

# Create global directory structure
globalpaths = (
            'lib/',
            'data/',
            'tmp/vectors/',
	    'tmp/plots/',
	    'tmp/htc/',
	    'tmp/injsrch/',
            'logs/',
            'config/',
            )

print 'Creating global file structure:'

for d in globalpaths:
    try:
        os.makedirs(d)
        print d
    except OSError:
        print 'Error: %r already exists' % d

# (maybe next section can be moved to script that imports the data)

# for det, runs in general.detruns.iteritems():
#     for run in runs:
#         p = 'globals/data/' + det + '/' + run
#         try:
#             os.makedirs(p)
#         except OSError:
#             print 'Error: "%(p)s" already exists' % locals()
