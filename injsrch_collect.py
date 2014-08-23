#! /usr/bin/env python 

import sys
from globs import results as res

'''
Collect results of single PSR analysis and export to public directory.
'''

processname, det, run, psr, kind, pdif = sys.argv

# setup log
#g.setuplog('results_%(det)s%(run)s_%(psr)s_%(kind)s%(pdif)s' % locals() )
#log = logging.getLogger('Results')


## SETUP

# Create results object
r = res.Results(det, run, psr, kind, pdif)

# Collect results
r.collect()

# Export results
r.export()

print 'Results for %(det)s %(run)s %(psr)s with %(kind)s %(pdif)s injections expported to:' % locals()
print r.paths['export']
