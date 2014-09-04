#! /usr/bin/env python 

import os, sys
#sys.path.append(os.path.expanduser('~') +'/polHTC/')
from lib import general as g

'''
Writes DAG file for a complete MANY-PSR analysis; namely, injecting GR, G4v on PSRs from
the pulsar list.
'''

##########################################################################################
# PARSE ARGUMENTS

args = sys.argv
if len(args) == 6:
    processname, det, run, psrIN, ninstSTR, ninjSTR = args
    analysis_name = 'full'
elif len(args) == 7:
    processname, det, run, psrIN, ninstSTR, ninjSTR, name = args
    analysis_name = 'full' + name
else:
    print 'Too many arguments to unpack.'
    sys.exit(1)


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

subname = g.paths['subs'] + 'generic_' + analysis_name + '.sub'

dagname = g.dag_path(det, run, psrIN, name=analysis_name)

config  = g.paths['dags'] + 'dagman.config'


##########################################################################################
# WRITE SUBMIT

home = os.path.expanduser('~')
project_dir = home + '/polHTC'

analysis_folder = '/analyses/full/' # this serves to separate unrelated executables

# get hostname to determine what server we are on
cluster = g.Cluster()

# define scratch space
scratch_dir = cluster.scratch_dir + analysis_name + '_$(psr)_$(det)$(run)_$(injkind)$(pdif)s'

subfile_lines = [
                'Universe = Vanilla',
                'Executable = ' + project_dir + analysis_folder + analysis_name + '.py',
                'initialdir = ' + project_dir + '/',
                'arguments = "$(psr) $(det) $(run) $(injkind) $(pdif) $(ninst) $(ninj)"' % locals(),
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
def lines(det, run, psr, injkind, pdif, ninstSTR, ninjSTR):
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
        'VARS %(jobname)s injkind="%(injkind)s"' % locals(),
        'VARS %(jobname)s pdif="%(pdif)s"' % locals(),
        'VARS %(jobname)s ninst="%(ninstSTR)s"' % locals(),
        'VARS %(jobname)s ninj="%(ninjSTR)s"' % locals(),
        'RETRY %(jobname)s 3' % locals(), # retry job 3 times
        '\n'
        ]
    
    return l


with open(dagname, 'w') as f:
    f.write('# %(det)s %(run)s %(psrIN)s %(ninstSTR)s %(ninjSTR)s\n\n' % locals() )
    f.write('CONFIG '+ config + '\n\n') # points to configuration file with DAGman variables
    
    for psr in goodpsrs:
            for injkind in ['GR', 'G4v']:
                for pdif in ['p']: #,'m']:
                    txt_lines = lines(det, run, psr, injkind, pdif, ninstSTR, ninjSTR)
                    for l in txt_lines:
                        f.write(l+'\n')


# Configure Dagman to not limit the number of times a node is put on hold
with open(config, 'w') as f:
    f.write('DAGMAN_MAX_JOB_HOLDS = 0')

print 'DAG written to: ' + dagname
print 'Submit using: condor_submit_dag %(dagname)s' % locals()
