#! /usr/bin/env python

"""Generates all OB plots for paper
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

latextable = {}
asciitable = {}
for d in detectors:
    for run in runs:
        detbest = (1, 'Z', 'Z')
        tex = [('PSR', 'p GR', 'p G4v', 'p Sid')]
        asciitable[d + run] = 'PSR ROT_FREQ H_DEP H_INDEP\n'

        for psr, pobs in obcontainer[d+run].iteritems():
            tex.append((psr, r'\num{%.3f}' % pobs['GR'],
                       r'\num{%.3f}' % pobs['G4v'],
                       r'\num{%.3f}' % pobs['Sid']))
            asciitable[d + run] += '%s %.3f %.3f %.3f\n' % (psr, pobs['GR'],
                                                          pobs['G4v'],
                                                          pobs['Sid'])
	    best = min([(pobs[m], m, psr) for m in ['GR', 'G4v', 'Sid']])
	    detbest = min(detbest, best)

        print 'Best OB %s %s: %r' % (run, d, detbest)
	latextable[d+run] = tabulate(tex, headers='firstrow', tablefmt='latex')
        latextable[d+run] = latextable[d+run].replace(r'\textbackslash{}', '\\').\
            replace(r'\{', r'{').replace(r'\}', r'}')

with open(tablepath, 'w') as f:
    for d in detectors:
        for run in runs:
            f.write('\n\n')
            f.write('%%' + d + run + ' sensitivity\n')
            f.write(latextable[d + run])
    print 'LaTex table saved: %r' % tablepath

for d in detectors:
    for run in runs:
        with open('tmp/plots/results%s%s.txt' % (d, run), 'w') as f:
            f.write(asciitable[d + run])
        print 'ASCII tables saved: tmp/plots/results%s%s.txt' % (d, run)
