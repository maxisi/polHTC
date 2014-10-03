#! /usr/bin/env python

import os
import sys
import numpy as np

""" Import PSR names from local data path.
"""
if len(sys.argv) == 3:
    process_name, det, run = sys.argv
    split = 1
else:
    process_name, det, run, split = sys.argv

destination_path = 'config/'
origin_path = 'data/%s/%s/' % (det, run)

try:
    os.makedirs(destination_path)
except OSError:
    print 'Warning: "%(destination_path)s" already exists' % locals()


print 'Importing PSR list from: ' + origin_path

pre = 'finehet_'
pst = '_' + det

# fnames is the list of files in the directory for which:
#   1. filename contains 'finehet'
#   2. filename ends with detector name
# (thus discarding aux files)
fnames = [x for x in os.listdir(origin_path) if 'finehet' in x]

# obtain names of PSRs in directory
allpsrs_set = set([f.strip(pre).strip('.hdf5').strip(det).strip('_')
                   for f in fnames])
#set removes duplicates
allpsrs = list(allpsrs_set)

# split list into groups if necessary
avg = len(allpsrs) / float(split)
psrgroups = []
last = 0.0
while last < len(allpsrs):
    psrgroups.append(allpsrs[int(last):int(last + avg)])
    last += avg

# save full list to file (this is good for PSRcat)
path = destination_path + 'psrlist_' + det + run + '.txt'

with open(path, 'a') as f:
    for psr in allpsrs:
        f.write(psr + '\n')

print str(len(allpsrs)) + ' PSRs for ' + det + ' ' + run + ' added to: ' + path

# save list to file
if split != 1:
    for i in np.arange(0,len(psrgroups)):
        path = destination_path + 'psrlist_' + det + '_' + run + '_' + str(i)\
               + '.txt'

        if i == 0:
            o_mode = 'w'
        else:
            o_mode = 'a'

        with open(path, o_mode) as f:
            for psr in psrgroups[i]: f.write(psr + '\n')

        print str(len(psrgroups[i])) + ' PSRs for ' + det + ' ' + run +\
              ' added to: ' + path