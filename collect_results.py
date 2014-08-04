#! /usr/bin/env python 

import sys
import logging
import general as g
from datetime import datetime

'''
Collect results of single PSR analysis and export to public directory.
'''

processname, det, run, psr, kind, pdif = sys.argv

# setup log
date = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
g.setuplog('%(date)s_results_%(det)s_%(run)s_%(psr)s_%(kind)s%(pdif)s' % locals() )
log = logging.getLogger('Results')


## SETUP

# Create results object
r = g.Results(det, run, psr, kind, pdif)

# Collect results
r.collect()

# Export results
r.export()

print 'Results for %(det)s %(run)s %(psr)s with %(kind)s %(pdf)s injections expported to:' % locals()
print r.paths['export']
