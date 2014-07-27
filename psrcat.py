#! /usr/bin/env python

import logging
from datetime import datetime
from bs4 import BeautifulSoup
from urllib2 import urlopen
import cPickle as pickle
import numpy as np

import general as g

'''
Gets location parameters for pulsars in input list from ATNF online catalogue.
Creates and pickles corresponding pandas DataFrame 'pulsar_parameters'.
'''

# SET UP LOG
date = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
g.setuplog(date + '_psrcat')

log = logging.getLogger('psrcat')


log.info('Loading PSR list.')

# open files with PSR list and read each line into the list 'psrs'.
# (assuming 1 psr per line)

psrs = []
try:
    f = open('config/psrlist.txt', 'r')
    for line in f.readlines(): psrs += [line.strip()] # (.strip() removes \n character)
    f.close()
except IOError:
    log.error('Could not open psrlist text in: config/psrlist.txt', exec_info=True)
    

log.info('Building PSR catalogue.')

def atnfurl(psr_list):
    pre = 'http://www.atnf.csiro.au/research/pulsar/psrcat/proc_form.php?version=1.47&JName=JName&RaJ=RaJ&DecJ=DecJ&startUserDefined=true&c1_val=&c2_val=&c3_val=&c4_val=&sort_attr=jname&sort_order=asc&condition=&'
    post = '&ephemeris=short&coords_unit=raj%2Fdecj&radius=&coords_1=&coords_2=&style=Long+with+errors&no_value=*&nohead=nohead&fsize=3&x_axis=&x_scale=linear&y_axis=&y_scale=linear&state=query&table_bottom.x=36&table_bottom.y=16'
    names= 'pulsar_names='

    for psr in psr_list:
        names+=psr.replace('+', '%2B')

        if psr != psr_list[-1]:
            names+='%0D%0A'
    url = "%(pre)s%(names)s%(post)s" % locals()

    if len(url)>2000:
        self.log.warning('WARNING! URL %d characters!' % len(url))

    return url

# Get data
url = atnfurl(psrs)

log.debug('Parsing webpage.')
soup = BeautifulSoup(urlopen(url)) # get webpage and parse
text = str(soup.pre.string)

log.debug('Creating catalogue dictionary.')
psrcat={}
for u in text.split('\n'): # take each line
    if u!='': # if the line is not empty
        # separate using spaces and remove empty elements
        a = [x for x in u.split(' ') if x !='']
        # add PSR entry to dictionary, containing parameter dictionary
        psrcat[a[1]] = {g.paramNames[n] : a[n] for n in np.arange(2,len(a))}
        # the array starts from 2 to exclude the list number of the pulsar and its name

log.debug('Reformatting RAS and DEC into decimals.')

for psr, paramdict in psrcat.iteritems(): # for each PSR
    for paramName, paramData in paramdict.iteritems(): # for each parameter
        # replace value with result of formatting function for parameter
        psrcat[psr][paramName] = g.paramFormat[paramName](paramData)
        
log.debug('Adding extra parameters.')

# try to open file with extra parameters
try:
    psrsExtra = [] # list of PSRs with extra parameters
    with open(g.paths['extrapsrparam'],'r') as f:
        # each line contains info for one PSR
        for row in f.readlines():
            # separate arguments
            extraParamList = row.split(',')
            psrsExtra += [extraParamList[0]]
            # if PSR (first element) in the list of pulsars
            if extraParamList[0] in psrs:
                # for each parameter in list, add element to PSR entry in psrcat dictionary
                for n in np.arange(1,len(g.extraParamNames)):
                    psrcat[extraParamList[0]][g.extraParamNames[n]] = float(extraParamList[n])
except IOError:
    log.warning('No extra PSR info found in: ' + g.paths['extrapsrparam'], exc_info=True)
    
# if there was no extra information, fill missing parameters with standard values
for psr in [src for src in psrs if src not in psrsExtra]: # for PSR in list for which there's no extra info
    for param in g.extraParamNames: # for each parameter
        if param!= None: psrcat[psr][param] = g.extraParamStandard[param]

try:
    pickle.dump(psrcat, open(g.paths['psrcat'], 'wb'))
except IOError:
    log.error('Could not save psrcat in: ' + g.paths['psrcat'], exec_info=True)