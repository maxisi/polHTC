#! /usr/bin/env python 

import sys
import general as g

pname, det, run, psr, ninstSTR, ninjSTR = sys.argv

'''
Writes DAG file for a complete one-PSR analysis; namely, injections of all kinds.
'''

dagname = g.dag_path(det, run, psr)

def lines(det, run, psr, injkind, pdif, ninstSTR, ninjSTR):
    # return DAG lines corresponding to one injection kind
    
    jobname = injkind + pdif
    
    l = [
        '# ' + injkind + pdif + '\n',
        'JOB ' + jobname + ' ' + g.submit_path(det, run, psr, injkind, pdif),
        'SCRIPT PRE %(jobname)s injsrch_master.py %(psr)s %(det)s %(run)s %(injkind)s %(pdif)s %(ninstSTR)s %(ninjSTR)s' % locals(),
        'SCRIPT POST %(jobname)s collect_results.py %(det)s %(run)s %(psr)s %(injkind)s %(pdif)s' % locals(),
        '\n'
        ]
    
    return l

with open(dagname, 'w') as f:
    for injkind in ['GR', 'G4v']:
        for pdif in ['p', 'm']:
            txt_lines = lines(det, run, psr, injkind, pdif, ninstSTR, ninjSTR)
            for l in txt_lines:
                f.write(l+'\n')

print 'DAG written to: ' + dagname