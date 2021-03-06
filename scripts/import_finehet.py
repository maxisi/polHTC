#! /usr/bin/env python

import os
import sys
import h5py
import numpy as np

'''
Import finehet data from M. Pitkin. Can take absolute path of CSV data as
argument. Otherwise assumes running on Atlas.
'''

process_name, det, run = sys.argv


# CREATE DESTINATION DIRECTORIES
# assuming main file structure exists, create subfolders for Data

destination_path = 'data/' + det + '/' + run + '/'

try:
    os.makedirs(destination_path)
except OSError:
    print 'Warning: %r already exists' % destination_path


# IMPORT DATA

def imp(det, run, origin_path, destination_path):
    """
    Imports det, run data for all pulsars in origin_path to destination_path
    (hdf5)
    """

    print('Importing finehet from: ' + origin_path + '\nPSRs: ')

    pre = 'finehet_'
    pst = '_' + det

    # fnames is the list of files in the directory for which:
    #   1. filename contains 'finehet'
    #   2. filename ends with detector name
    # (thus discarding aux files)
    fnames = [x for x in os.listdir(origin_path) if
              ('finehet' in x and x.find(det) == (len(x)-2))]

    # obtain names of PSRs in directory
    allpsrs = [f.strip(pre).strip(det).strip('_') for f in fnames]

    imported = []

    # open each finehet file and import to hdf5
    for psr in allpsrs:

        filename = pre + psr + pst

        try:
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
                    f_hdf5 = h5py.File(destination_path + filename + '.hdf5',
                                       'w')

                    f_hdf5.create_dataset('time', data=np.array(t,
                                                                dtype=np.int))
                    f_hdf5.create_dataset('data', data=d, dtype=np.complex_)

                    imported += [psr]  # mark as successfully imported
                    print(psr + ', ')

                except IOError:
                    print 'Could not write finehet to HDF5:'
                    print sys.exc_info()[0]
                    sys.exit(1)
        except:
            print '\nCould not import PSR %s from %r' % (psr, origin_path + filename)
	    print sys.exc_info()


    print '\n' + det + ' ' + run + ' data for ' + str(len(imported)) +\
          ' PSRs imported to: ' + destination_path


if run == 'S5':

    origin_path = '/home/matthew/analyses/S5/V4/fine/%s/total/' % det
    imp(det, run, origin_path, destination_path)

    # Data for some PSRs needs to be replaced with data from another directory
    origin_path2 = '/home/matthew/analyses/S5/extra/results/data%s/' % det
    imp(det, run, origin_path2, destination_path)

elif run in ['S6', 'VSR']:

    root_path = '/home/matthew/analyses/S6_all/results/'
    root_end = 'data%s/' % det
    o_paths = [
        root_path + root_end,
        root_path + 'J0534+2200/' + root_end,
        root_path + 'J0537-6910/full/' + root_end,
        root_path + 'J2022+3842/' + root_end,
    ]
    if run == 'VSR':
        o_paths += [root_path + 'J1833-1034/' + root_end]

    # Special cases (inco/indep):
    # J1952+3252
    # J0835-4510
    # J1813-1246
    # J0537-6910 (took full)

    for p in o_paths:
        imp(det, run, p, destination_path)

else:
    print 'Did not recognize "' + run + '" as a valid run'
