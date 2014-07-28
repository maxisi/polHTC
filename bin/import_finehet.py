#! /usr/bin/env python 

import os
import sys
import h5py
import numpy as np

'''
Import finehet data from M. Pitkin. Can take absolute path of CSV data as argument.
Otherwise assumes running on Atlas.
'''

process_name, det, run = sys.argv


# CREATE DESTINATION DIRECTORIES
# assuming main file structure exists, create subfolders for Data

destination_path = 'globals/data/' + det + '/' + run + '/'

try:
    os.makedirs(destination_path)
except OSError:
    print 'Warning: "%(p)s" already exists' % locals()


# IMPORT DATA

def imp(det, run, origin_path, destination_path):
    '''
    Imports det, run data for all pulsars in origin_path to destination_path (hdf5)
    '''
    
    pre = 'finehet_'
    pst = '_' + det
    
    # fnames is the list of files in the directory for which:
    #   1. filename contains 'finehet'
    #   2. filename ends with detector name
    # (thus discarding aux files)
    fnames = [x for x in os.listdir(destination_path) if ('finehet' in x and x.find(det)==(len(x)-2))]
    
    # obtain names of PSRs in directory
    allpsrs = [f.strip(pre).strip(pst) for f in fnames]
    
    # open each finehet file and import to hdf5
    for psr in allpsrs:
        print psr #
        filename = pre + psr + pst
        
        with open(origin_path + filename, 'r') as f_txt:
            
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
                f_hdf5 = h5py.File(destination_path + filename + '.hdf5', 'w')

                f_hdf5.create_dataset('time', data=np.array(t, dtype=np.int))
                f_hdf5.create_dataset('data', data=d, dtype=np.complex_)
                
                print det + run + 'data for' + str(len(allpsrs)) + ' PSRs imported to: ' + destination_path
                
            except:
                print 'Could not write finehet to HDF5:'
                print sys.exc_info()[0]
                sys.exit()


if run == 'S5':

    origin_path = '/home/matthew/analyses/S5/V4/fine/' + det + '/total/'

    imp(det, run, origin_path, destination_path)

elif run in ['S6', 'VSR']:

    root_path = '/home/matthew/analyses/S6_all/results/'
    root_end  = 'data' + det + '/',
    o_paths = [
                root_path + root_end,
                root_path + 'J0534+2200/' + root_end,
                root_path + 'J0537-6910/full/'+ root_end,
                root_path + 'J2022+3842/'+ root_end,
                ]
    if run == 'VSR': o_paths += [root_path + 'J1833-1034/' + root_end]

    # Special cases (inco/indep):
    # J1952+3252
    # J0835-4510
    # J1813-1246
    # J0537-6910 (took full)
    
    for p in o_paths: imp(det, run, p, destination_path)