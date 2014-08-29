#! /usr/bin/env python 

import os
import sys
from globs import general as g

processname, det, run, psrIN, gridsize = sys.argv

'''
Writes DAG file for a complete OB analysis; namely, searching for GR, G4v & Sid on PSRs from
the pulsar list.
'''

##########################################################################################

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

##########################################################################################
# PATHS

subname = 'subs/ob.sub'

dagname = 'dags/ob_%(det)s%(run)s_%(psrIN)s.dag' % locals()


##########################################################################################
# WRITE SUBMIT

home = os.path.expanduser('~')
project_dir = home + '/polHTC'

analysis_folder = '/analysis_ob' # this serves to separate unrelated executables

# get hostname to determine what server we are on
cluster = g.Cluster()

# define scratch space
scratch_dir = cluster.scratch_dir + 'ob_$(psr)_$(det)$(run)s'

subfile_lines = [
                'Universe = Vanilla',
                'Executable = ' + project_dir + analysis_folder + '/ob_method.py',
                'initialdir = ' + project_dir + '/',
                'arguments = "$(psr) $(det) $(run) $(method) $(gridsize)"' % locals(),
                'Output = ' + scratch_dir + '.out',
                'Error = ' + scratch_dir + '.err',
                'Log = ' + scratch_dir + '.log',
                'getenv = true',
                'Queue'
                ]

with open(subname, 'w') as f:
    for l in subfile_lines: f.write(l + '\n')


##########################################################################################
# WRITE DAG

# lines() helps write DAG
def lines(det, run, psr, method, gridsize):
    # return DAG lines corresponding to one injection kind

    home = os.path.expanduser('~')
    project_dir = home + '/polHTC/'

    jobname = injkind + pdif + '_' + psr

    l = [
        '# ' + psr + ' ' + injkind + pdif + '\n',
        'JOB ' + jobname + ' ' + project_dir + subname,
        'VARS %(jobname)s psr="%(psr)s"' % locals(),
        'VARS %(jobname)s det="%(det)s"' % locals(),
        'VARS %(jobname)s run="%(run)s"' % locals(),
        'VARS %(jobname)s method="%(method)s"' % locals(),
        'VARS %(jobname)s gridsize="%(gridsize)s"' % locals(),
        'RETRY %(jobname)s 3' % locals(), # retry job 3 times
        '\n'
        ]
    
    return l


with open(dagname, 'w') as f:
    f.write('# %(det)s %(run)s %(psrIN)s %(method)s %(gridsize)s\n\n' % locals() )
    f.write('CONFIG dags/dagman.config\n\n') # points to configuration file with DAGman variables
    
    for psr in goodpsrs:
            for method in ['GR', 'G4v', 'Sid']:
                for l in lines(det, run, psr, method, gridsize):
                    f.write(l+'\n')


# Configure Dagman to not limit the number of times a node is put on hold
with open('dags/dagman.config', 'w') as f:
    f.write('DAGMAN_MAX_JOB_HOLDS = 0')

print 'DAG written to: ' + dagname
print 'Submit using: condor_submit_dag %(dagname)s' % locals()
