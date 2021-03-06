#! /usr/bin/env python

import sys
from polHTC import results as res

'''
Collect results of single PSR analysis and export to public directory.
'''

processname, det, run, psr, kind = sys.argv

# setup log
#g.setuplog('results_%(det)s%(run)s_%(psr)s_%(kind)s' % locals() )
#log = logging.getLogger('Results')


## SETUP

# Create results object
r = res.Results(det, run, psr, kind)

# Collect results
r.collect()

# Export results
r.export()

print 'Results for %s %s with %s %s injections exported to:' \
      % (det, run, psr, kind)
print r.paths['export']
