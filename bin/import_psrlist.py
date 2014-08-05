#! /usr/bin/env python 

import os
import sys
import h5py
import numpy as np

'''
Import PSR names from M. Pitkin. Can take absolute path of CSV data as argument.
Otherwise assumes running on Atlas.
'''
if len(sys.argv) == 3:
    process_name, det, run = sys.argv
    split = 1
else:
    process_name, det, run, split = sys.argv

# CREATE DESTINATION DIRECTORIES
# assuming main file structure exists, create subfolders for Data

destination_path = 'config/'

try:
    os.makedirs(destination_path)
except OSError:
    print 'Warning: "%(destination_path)s" already exists' % locals()


# IMPORT NAMES

def imp(det, run, split, origin_path, destination_path):
    '''
    Imports det, run PSR names for pulsars in origin_path to destination_path.
    '''
    
    print 'Importing PSR list from: ' + origin_path
    
    pre = 'finehet_'
    pst = '_' + det
    
    # fnames is the list of files in the directory for which:
    #   1. filename contains 'finehet'
    #   2. filename ends with detector name
    # (thus discarding aux files)
    fnames = [x for x in os.listdir(origin_path) if ('finehet' in x and x.find(det)==(len(x)-2))]

    # obtain names of PSRs in directory
    allpsrs = [f.strip(pre).strip(pst) for f in fnames]
    
    # split list into groups if necessary
    avg = len(allpsrs) / float(split)
    psrgroups = []
    last = 0.0
    while last < len(allpsrs):
        psrgroups.append(allpsrs[int(last):int(last + avg)])
        last += avg
        
    # save full list to file (this is good for PSRcat)
    path = destination_path + 'psrlist.txt'
    
    with open(path, 'a') as f:
        for psr in allpsrs: f.write(psr + '\n')

    print str(len(psrgroups[i])) + ' PSRs for ' + det + ' ' + run + ' added to: ' + path
    
    # save list to file
    if split!=1:
        for i in np.arange(0,len(psrgroups)):
            path = destination_path + 'psrlist_' + run + '_' + str(i) + '.txt'
        
            if i == 0:
                o_mode = 'w'
            else:
                o_mode = 'a'

            with open(path, o_mode) as f:
                for psr in psrgroups[i]: f.write(psr + '\n')
    
            print str(len(psrgroups[i])) + ' PSRs for ' + det + ' ' + run + ' added to: ' + path


if run == 'S5':

    origin_path = '/home/matthew/analyses/S5/V4/fine/' + det + '/total/'

    imp(det, run, split, origin_path, destination_path)

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
    #   J1952+3252
    #   J0835-4510
    #   J1813-1246
    #   J0537-6910 (took full)
    
    for p in o_paths: imp(det, run, split, p, destination_path)

else:
    print 'Did not recognize "' + run + '" as a valid run'
