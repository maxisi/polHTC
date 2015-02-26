#! /usr/bin/env python

import os
import sys

sys.path.append(os.path.expanduser('~') +'/polHTC/')
from polHTC import general as g

processname, det, run, psrIN, gridsize = sys.argv

'''
Writes DAG file for a complete OB analysis; namely, searching for GR, G4v & Sid on PSRs from
the pulsar list.
'''

##########################################################################################

# determine what PSRs to analyze from argument
psrlist= g.read_psrlist(det=det, run=run, name=psrIN)

# load PSR exclusion list (if it exists):
badpsrs = g.read_psrlist(name='bad')

goodpsrs = list( set(psrlist) - set(badpsrs) )

##########################################################################################
# PATHS

analysis_name = 'ob'

analysis_folder = '/analyses/ob/' # this serves to separate unrelated executables

home = os.path.expanduser('~')
project_dir = home + '/polHTC'

# define scratch space
scratch_dir = g.Cluster().scratch_dir + analysis_name + '_$(psr)_$(det)$(run)s'

subname =  g.paths['subs'] + analysis_name + '.sub'

dagname = g.dag_path(det, run, psrIN, name=analysis_name)

config  = g.paths['dags'] + 'dagman.config'


##########################################################################################
# WRITE SUBMIT

subfile_lines = [
                'Universe = Vanilla',
                'Executable = ' + project_dir + analysis_folder + 'exe.py',
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

    jobname = method + '_' + psr

    l = [
        '# ' + psr + ' ' + method + '\n',
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
    f.write('# %(det)s %(run)s %(psrIN)s %(gridsize)s\n\n' % locals() )
    f.write('CONFIG '+ config + '\n\n') # points to configuration file with DAGman variables

    for psr in goodpsrs:
            for method in ['GR', 'G4v', 'Sid']:
                for l in lines(det, run, psr, method, gridsize):
                    f.write(l+'\n')


# Configure Dagman to not limit the number of times a node is put on hold
with open(config, 'w') as f:
    f.write('DAGMAN_MAX_JOB_HOLDS = 0')

print 'DAG written to: ' + dagname
print 'Submit using: condor_submit_dag %(dagname)s' % locals()
