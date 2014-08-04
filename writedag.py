#! /usr/bin/env python 

import sys
import general as g

pname, det, run, psrIN, ninstSTR, ninjSTR = sys.argv

'''
Writes DAG file for a complete one-PSR analysis; namely, injections of all kinds.
'''

dagname = g.dag_path(det, run, psrIN)

# determine what PSRs to analyze from argument
if psrIN == 'all':
    # assuming all pulsars in psrlist should be included
    psrlist = g.read_psrlist()
else:
    # assuming one single psr was requested
    psrlist = [psrIN]

# lines() helps write DAG
def lines(det, run, psr, injkind, pdif, ninstSTR, ninjSTR):
    # return DAG lines corresponding to one injection kind
    
    jobname = injkind + pdif
    
    l = [
        '# ' + injkind + pdif + '\n',
        'JOB ' + jobname + ' ' + g.submit_path(det, run, psr, injkind, pdif),
        'SCRIPT PRE %(jobname)s injsrch_master.py %(psr)s %(det)s %(run)s %(injkind)s %(pdif)s %(ninstSTR)s %(ninjSTR)s' % locals(),
        'SCRIPT POST %(jobname)s injsrch_collect.py %(det)s %(run)s %(psr)s %(injkind)s %(pdif)s' % locals(),
        '\n'
        ]
    
    return l

# write dag
for psr in psrlist:
    with open(dagname, 'w') as f:
        for injkind in ['GR', 'G4v']:
            for pdif in ['p', 'm']:
                txt_lines = lines(det, run, psr, injkind, pdif, ninstSTR, ninjSTR)
                for l in txt_lines:
                    f.write(l+'\n')

print 'DAG written to: ' + dagname
print 'Submit using: condor_submit_dag -maxidle 5000 %(dagname)s' % locals()
