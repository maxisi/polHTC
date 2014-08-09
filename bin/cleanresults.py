#! /usr/bin/env python

# Removes exported results. If a an integer is provided, removes results corresponding to pulsars IN THAT LIST ONLY.

import sys
import os
import numpy as np
from globals import general as g

if len(sys.argv)==1:
    # deleting results for ALL pulsars in main list
    list_name = ''

else:
    # assuming first argument is list name
    list_name = sys.argv[1]
    
# Import psr list 
psr_list = g.read_psrlist(list_name)

# Get path were files are located
c = g.Cluster()

# Get names of files in path
fnames = [x for x in os.listdir(c.public_path) if 'results' in x]

# Get PSR names from filenames
psr_file = []
delete_list = []
for f in fnames:
    for det in ['H1', 'H2', 'L1', 'V1']:
        for run in (g.ligoruns + g.virgoruns):
            for inj in ['GR', 'G4V']:
                for p in ['p', 'm', '0']:
                    psr = f.strip('results_%(det)s%(run)s_' % locals())
                    psr = psr.strip('_%(inj)s%(p)s.hdf5' % locals() )
                    psr_file += [psr]
                    delete_list += [f]

# Get confirmation from user.
message = 'Removing ' + str(len(delete_list)) + ' result files from ' + c.public_path + ' corresponding to list ' + list_name + '. Confirm? (y/n)'

proceed = input(message)

if proceed in ['y', 'yes']:
    # Delete files in delete_list
    for i in np.arange(0, len(psr_file)):
        if psr_file[i] in psr_list:
            os.remove(c.public_path + delete_list[i])
