#! /usr/bin/env python 

'''
Creates basic file structure for polarization analysis (if not already existent).
'''

import os

# Create global directory structure
globalpaths = (
            'logs/',
            'globs/data/',
            'globs/vectors/',
            'config/',
            'subs/',
            'dags/'
            )
            
print 'Creating global file structure:'
    
for dir in globalpaths:
    try:
        os.makedirs(dir)
        print dir
    except OSError:
        print 'Error: "globals/%(dir)s" already exists' % locals()
        
# (maybe next section can be moved to script that imports the data)

# for det, runs in general.detruns.iteritems():
#     for run in runs:
#         p = 'globals/data/' + det + '/' + run
#         try:
#             os.makedirs(p)
#         except OSError:
#             print 'Error: "%(p)s" already exists' % locals()