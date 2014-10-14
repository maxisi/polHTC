#! /usr/bin/env python

"""Generates all plots for paper
"""
import os
import sys
from tabulate import tabulate

sys.path.append(os.getcwd())
from lib import results as r

###############################################################################
# CONFIG
p = '/home/misi/results/'
tablepath = 'tmp/plots/tableOB.txt'

# # ALL PSRs RESULTS PLOTS
detectors = ['H1', 'L1']#, 'H2']
runs = ['S5']#,'S6']
injections = ['GR', 'G4v']

obcontainer = {}
for d in detectors:
    for run in runs:
        obcontainer[d+run] = {}
        # create MP results object (enough to do for 1 injkind when getting p)
        results = r.ResultsMP(injections[0], det=d, run=run, path=p)
        # operations need be done only once:
        # plot
        results.plot_ob('fgw', 'p_ob', title=False, ylim=(1e-3, 1),
                        log=True)
        # obtain PSR list
        for n in range(len(results.psrlist)):
            obcontainer[d+run][results.psrlist[n]] = {}
        # create hmin table
        pob = results.openboxes(det_thrsh=.999, p_fitorder=2)[2]
        for injkind in injections:
            for n in range(len(results.psrlist)):
                obcontainer[d+run][results.psrlist[n]][injkind] = pob[injkind][n]
                if injkind == 'GR':
                    obcontainer[d+run][results.psrlist[n]]['Sid'] = pob['Sid'][n]

tableob = {}
for d in detectors:
    for run in runs:
        t = [('PSR', 'p GR', 'p G4v', 'p Sid')]
        for psr, pobs in obcontainer[d+run].iteritems():
            t.append((psr, r'\num{%.3f}' % pobs['GR'],
                      r'\num{%.3f}' % pobs['G4v'],
                      r'\num{%.3f}' % pobs['Sid']))
        tableob[d+run] = tabulate(t, headers='firstrow', tablefmt='latex')
        tableob[d+run] = tableob[d+run].replace(r'\textbackslash{}', '\\').\
            replace(r'\{', r'{').replace(r'\}', r'}')

with open(tablepath, 'w') as f:
    for d in detectors:
        for run in runs:
            f.write('\n\n')
            f.write('%%' + d + run + ' sensitivity\n')
            f.write(tableob[d + run])
    print 'Table saved: %r' % tablepath