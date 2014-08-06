#! /usr/bin/env python 

import os
import sys
import general as g

pname, det, run, psrIN, ninstSTR, ninjSTR = sys.argv

'''
Writes DAG file for a complete one-PSR analysis; namely, injections of all kinds.
'''

dagname = g.dag_path(det, run, psrIN)

# determine what PSRs to analyze from argument
try:
    # Determine whether psrIN is a chunk index (e.g. '2'). 
    int(psrIN)
    
    # If it is, assume the requested PSRs are those in the list named psrlist_run_psrIN.txt
    psrlist = g.read_psrlist(name = det + '_' + run + '_' + psrIN)

except ValueError:
    if psrIN == 'all':
        # assuming all pulsars in psrlist should be included
        psrlist = g.read_psrlist()
    else:
        # assuming one single psr was requested
        psrlist = [psrIN]

# load PSR exclusion list (if it exists):
try:
    with open(g.paths['badpsrs'], 'r') as f:
        badpsrs=[]
        for line in f.readlines():
            badpsrs += [line.strip()] # (.strip() removes \n character)
except:
    print 'WARNING: no PSR exclusion list found'

goodpsrs = list( set(psrlist) - set(badpsrs) )

# lines() helps write DAG
def lines(det, run, psr, injkind, pdif, ninstSTR, ninjSTR):
    # return DAG lines corresponding to one injection kind

    home = os.path.expanduser('~')
    project_dir = home + '/polHTC/'

    jobname = injkind + pdif + '_' + psr
    
    l = [
        '# ' + psr + ' ' + injkind + pdif + '\n',
        'JOB ' + jobname + ' ' + project_dir + g.submit_path(det, run, psr, injkind, pdif),
        #'SCRIPT PRE %(jobname)s %(project_dir)sinjsrch_master.py %(psr)s %(det)s %(run)s %(injkind)s %(pdif)s %(ninstSTR)s %(ninjSTR)s' % locals(),
        'SCRIPT POST %(jobname)s %(project_dir)sinjsrch_collect.py %(det)s %(run)s %(psr)s %(injkind)s %(pdif)s' % locals(),
        'CATEGORY %(jobname)s analysis' % locals(),
        '\n'
        ]
    
    return l

# write dag
with open(dagname, 'w') as f:
    for psr in goodpsrs:
            for injkind in ['GR', 'G4v']:
                for pdif in ['p']: #,'m']:
                    txt_lines = lines(det, run, psr, injkind, pdif, ninstSTR, ninjSTR)
                    for l in txt_lines:
                        f.write(l+'\n')
    f.write('MAXJOBS analysis 1\n')
        # this prevents more than 1 jobs to be submitted at once, limiting the max numb of
        # queued processes to 1 * ninst

print 'DAG written to: ' + dagname
print 'Submit using: condor_submit_dag %(dagname)s' % locals()
