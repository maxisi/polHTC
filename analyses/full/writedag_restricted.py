#! /usr/bin/env python

import os
import sys

sys.path.append(os.path.expanduser('~') + '/polHTC/')
sys.path.append(os.getcwd())
from polHTC import general as g

"""
Writes DAG file for a complete MANY-PSR analysis; namely, injecting GR, G4v on
PSRs from the pulsar list.
"""

###############################################################################
# PARSE ARGUMENTS

args = sys.argv
if len(args) == 6:
    processname, det, run, psrIN, ninstSTR, ninjSTR, injtype = args

    name = ''
elif len(args) == 7:
    processname, det, run, psrIN, ninstSTR, ninjSTR, injtype, name = args

else:
    print 'Too many arguments to unpack.'
    sys.exit(1)

analysis_name = 'full' + name

if injtype == "both":
    injkinds = ["GR", "G4v"]
else:
    injkinds = [injtype]

# determine what PSRs to analyze from argument
psrlist = g.read_psrlist(det=det, run=run, name=psrIN)

# load PSR exclusion list (if it exists):
badpsrs = g.read_psrlist(name='bad')
print type(psrlist)
print type(badpsrs)
goodpsrs = list(set(psrlist) - set(badpsrs))

###############################################################################
# PATHS

subname = g.paths['subs'] + analysis_name + '.sub'

dagname = g.dag_path(det, run, psrIN, name=analysis_name)

config  = g.paths['dags'] + 'dagman.config'


###############################################################################
# WRITE SUBMIT

home = os.path.expanduser('~')
project_dir = os.path.join(home, 'polHTC/')

analysis_folder = '/analyses/full/'
# (this serves to separate unrelated executables)

# get hostname to determine what server we are on:
cluster = g.Cluster()

# define scratch space
scratch_dir = cluster.scratch_dir + analysis_name +\
              '_$(psr)_$(det)$(run)_$(injkind)'

subfile_lines = [
    'Universe = Vanilla',
    'Executable = ' + project_dir + analysis_folder + 'exe' + name + '.py',
    'initialdir = ' + project_dir + '/',
    'arguments = "$(psr) $(det) $(run) $(injkind) $(ninst) $(ninj)"',
    'Output = ' + scratch_dir + '.out',
    'Error = ' + scratch_dir + '.err',
    'Log = ' + scratch_dir + '.log',
    'getenv = true',
    'request_memory = 1 GB',
    'Queue'
]

with open(subname, 'w') as f:
    for l in subfile_lines:
        f.write(l + '\n')


###############################################################################
# WRITE DAG

with open(dagname, 'w') as f:
    # Write initial comment:
    f.write('# %(det)s %(run)s %(psrIN)s %(ninstSTR)s %(ninjSTR)s\n\n'
            % locals())
    # Point to configuration file with DAGman variables:
    f.write('CONFIG ' + config + '\n\n')

    for psr in goodpsrs:
            for injkind in injkinds:
                jobname = injkind + '_' + psr

                txt_lines = [
                    '# ' + psr + ' ' + injkind + '\n',
                    'JOB %s %s' % (jobname, subname),
                    'VARS %s psr="%s"' % (jobname, psr),
                    'VARS %(jobname)s det="%(det)s"' % locals(),
                    'VARS %(jobname)s run="%(run)s"' % locals(),
                    'VARS %(jobname)s injkind="%(injkind)s"' % locals(),
                    'VARS %(jobname)s ninst="%(ninstSTR)s"' % locals(),
                    'VARS %(jobname)s ninj="%(ninjSTR)s"' % locals(),
                    'RETRY %(jobname)s 3' % locals(),  # retry job 3 times
                    '\n'
                ]

                for l in txt_lines:
                    f.write(l + '\n')


# Configure Dagman to not limit the number of times a node is put on hold
with open(config, 'w') as f:
    f.write('DAGMAN_MAX_JOB_HOLDS = 0')

print 'DAG written to: ' + dagname
print 'Submit using: condor_submit_dag %s' % dagname
