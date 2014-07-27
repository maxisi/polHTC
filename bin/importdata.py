#! /usr/bin/env python 

import os
import sys
import h5py
import numpy as np

'''
Import finehet data from M. Pitkin. Can take absolute path of CSV data as argument.
Otherwise assumes running on Atlas.
'''

# DETERMINE WHAT ORIGIN DIRECTORY TO USE

if len(sys.argv)>1:
    # assuming arguments are: detector run path
    detectors = [sys.argv[1]]
    run = sys.argv[2]
    directory = sys.argv[3]
    
    print 'Importing finehet %(detectors)s %(run)s data from $(directory)s\n' % locals()
    
else:
    detectors = ['H1','L1','V1']
    directory = g.paths['originalData'] + '/data'
    print 'Assuming ATLAS. Importing %(directory)s for H1, L1, V1' % locals()


# CREATE DIRECTORIES
# assuming main file structure exists, create subfolders for Data

ligoruns = ('S5', 'S6')
virgoruns = ('S1', 'S2')

detruns = {
            'H1' : ligoruns,
            'H2' : ligoruns,
            'L1' : ligoruns,
            'V1' : virgoruns,
            }

for det, runs in detruns.iteritems():
    for run in runs:
        p = 'globals/data/' + det + '/' + run
        try:
            os.makedirs(p)
        except OSError:
            print 'Warning: "%(p)s" already exists' % locals()


# IMPORT DATA

pre = 'finehet_'
pst = '_'

# loop over detectors
for det in detectors:
    
    path = directory
    
    # in ATLAS case append detname to path and set run to default
    if len(detectors)>1:
        path += det
        
        if det == 'V1':
            run = 'S1'
        else:
            run = 'S6'
    
    # fnames is the list of files in the directory for which:
    #   1. filename contains 'finehet'
    #   2. filename ends with detector name
    fnames = [x for x in os.listdir(path) if ('finehet' in x and x.find(det)==(len(x)-2))]
    
    # obtain names of PSRs in directory
    allpsrs = [f.strip(pre).strip(pst + det) for f in fnames]
    
    # open each finehet file and import to hdf5
    for psr in allpsrs:
        filename = pre + psr + pst + det
        with open(path + '/' + filename, 'r') as f_txt:
            
            # extract data from file
            t = []
            re = []
            im = []
            
            for row in f_txt.readlines():
                rowlist = row.split('\t')
                t += [float(rowlist[0])]
                re += [float(rowlist[1])]
                im += [float(rowlist[2])]
            
            d = np.array(re) + 1j*np.array(im, dtype=np.float64)
                        
            # create hdf5 file
            try:
                f_hdf5 = h5py.File('globals/data/'+det+'/'+run+'/'+filename+'.hdf5', 'w')

                f_hdf5.create_dataset('time', data=np.array(t, dtype=np.int))
                f_hdf5.create_dataset('data', data=d, dtype=np.complex_)

            except:
                print 'Could not write finehet to HDF5:'
                print sys.exc_info()[0]
                sys.exit()