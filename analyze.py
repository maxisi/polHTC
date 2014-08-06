#! /usr/bin/env python 

'''
Pre-processes Condor DAG submission:
    - Calls injsrch_master for each analysis
    - Calls makedag
    - Submits
'''

import os
import sys

import general as g

pname, det, run, psrIN, ninstSTR, ninjSTR = sys.argv

##########################################################################################
# Determine what PSRs to analyze from argument

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

print 'Analizying ' + str(len(goodpsrs)) + ' PSRs.'


##########################################################################################
# Execute injsrch_master

print 'Pre-processing.'

for psr in goodpsrs:
        for injkind in ['GR', 'G4v']:
            for pdif in ['p']: #,'m']:
                master_call = './injsrch_master.py %(psr)s %(det)s %(run)s %(injkind)s %(pdif)s %(ninstSTR)s %(ninjSTR)s' % locals()
                os.system(master_call)


##########################################################################################
# Execute writedag

print 'Writing DAG.'

dag_call = './writedag.py %(det)s %(run)s %(psrIN)s %(ninstSTR)s %(ninjSTR)s' % locals() 
os.system(dag_call)


##########################################################################################
# Submit
# 
# print 'Submitting.'
# 
# dagname = g.dag_path(det, run, psrIN)
# os.system('condor_submit_dag %(dagname)s' % locals())