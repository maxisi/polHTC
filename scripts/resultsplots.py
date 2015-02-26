#! /usr/bin/env python

"""Generates all plots for paper
"""
import os
import sys
from tabulate import tabulate

sys.path.append(os.getcwd())
from polHTC import results as r

###############################################################################
# CONFIG
p = '/home/misi/results/'

###############################################################################
# METHODS PLOTS (Crab)
crab_gr = r.Results('H1', 'S5', 'J0534+2200', 'GR')
crab_gr.load(path=p)
# s vs hrec
crab_gr.plot_hs(title=False)
# p--value
crab_gr.plot_p('s', title=False, star=10, starsize=12,
            xlim=(0, 12), legend=False, methods=['G4v'])
# s/hrec vs hinj example
crab_gr.plot('s', methods=['GR'], title=False, legend=False, suffix='_ex')
crab_gr.plot('h', methods=['GR'], title=False, legend=False, suffix='_ex',
             ylim=(0,1.2e-24),)

###############################################################################
# CRAB RESULTS PLOTS

# GR
crab_gr.plot('s', methods=['GR', 'Sid'], band_conf=0, title=False,
             suffix='_sid')
crab_gr.plot('s', methods=['GR', 'G4v'], band_conf=0, title=False,
             suffix='_g4v')
# G4v
crab_g4v = r.Results('H1', 'S5', 'J0534+2200', 'G4v')
crab_g4v.load(path=p)

crab_g4v.plot('s', methods=['G4v', 'Sid'], band_conf=0, title=False,
             suffix='_sid')
crab_g4v.plot('s', methods=['GR', 'G4v'], band_conf=0, title=False,
             suffix='_gr')

cgrob = crab_gr.openbox(methods=['GR', 'Sid'], det_thrsh=.999, p_fitorder=2)
cg4vob = crab_g4v.openbox(methods=['G4v', 'Sid'], det_thrsh=.999, p_fitorder=2)

crab_gr.plot_p('s', methods=['GR','Sid'], star=cgrob[1], fit=0, manyfiles=True)
crab_g4v.plot_p('s', methods=['G4v','Sid'], star=cg4vob[1], fit=0, manyfiles=1)

###############################################################################
# ALL PSRs RESULTS PLOTS
detectors = ['H1', 'L1']#, 'H2']
runs = ['S5']#,'S6']
injections = ['GR', 'G4v']
mpbest = {}
mp_rho = {}
cr_rho = {}

for d in detectors:
    for run in runs:
        mpbest[d + run] = [(' ', 'PSR', r'$\nu$', r'$h_{\rm dep}$',
                            r'$h_{\rm indep}$')]
	mp_rho[d + run] = [(' ', 'GR', 'G4v', 'Sid')]
	cr_rho[d + run] = [(' ', 'GR', 'G4v', 'Sid')]

        for injkind in injections:
            results = r.ResultsMP(injkind, det=d, run=run, path=p)
            # produce plot
            results.plot_ref()

            if d == 'H1':
                results.plot('fgw', 's_slope', methods=['GR', 'G4v', 'Sid'],
                             logx=1, logy=1, title=0, xlim=(30, 1500))
                results.plot('fgw', 's_noise', methods=['GR', 'G4v', 'Sid'],
                             logx=1, title=0, xlim=(30, 1500))

            # create best-hmin tables
            names, hmin = results.sortby('psrlist', 'hmin')
            freq, _ = results.sortby('FR0', 'hmin')

            # LaTex
            mpbest[d + run].append((injkind, names[injkind][0],
                                     r'\num{%.2f}' % freq[injkind][0],
                                     r'\num{%.3g}' % hmin[injkind][0],
                                     r'\num{%.3g}' % hmin['Sid'][0]))

    	    # scale factors
	    sf_avg = results.scalefactor(methods=['GR', 'G4v', 'Sid'])
	    sf_crab = results.scalefactor(psr='J0534+2200',
	                methods=['GR', 'G4v', 'Sid'])

            mp_rho[d+run].append((injkind,
                            r'\num{%.2f}' % sf_avg['GR'],
		    	    r'\num{%.2f}' % sf_avg['G4v'],
			    r'\num{%.2f}' % sf_avg['Sid']))
            cr_rho[d+run].append((injkind,
                            r'\num{%.2f}' % sf_crab['GR'],
		   	    r'\num{%.2f}' % sf_crab['G4v'],
			    r'\num{%.2f}' % sf_crab['Sid']))


###############################################################################
# GAUSSIAN NOISE RESULTS PLOTS

#for d in detectors:
#    for run in runs:
#        for injkind in injections:
#            results = r.ResultsMP(injkind, det=d, run=run)
#            results.load(path=p, prefix='gauss_')
#            # produce plot
#            results.plot_ref()

###############################################################################
# TABLES

# Crab sensitivity
crab_gr_hmin = crab_gr.min_h_det(.999)
crab_g4v_hmin = crab_g4v.min_h_det(.999)

crabhmin = {' ': ('GR', 'G4v')}
for m in ['GR', 'G4v', 'Sid']:
    crabhmin[m] = (r'\num{%.3g}' % crab_gr_hmin[m],
                   r'\num{%.3g}' % crab_g4v_hmin[m])

t1 = tabulate(crabhmin, headers=crabhmin.keys(), tablefmt='latex')
t1 = t1.replace(r'\textbackslash{}', '\\').replace(r'\{', r'{').replace(r'\}',
                                                                        r'}')

# Crab sensitivity relative to matching template
crabhmin_rel = {' ': ('GR', 'G4v')}
for m in ['GR', 'G4v', 'Sid']:
    crabhmin_rel[m] = (r'num{%.2f}' % (crab_gr_hmin[m]/crab_gr_hmin['GR']),
                       r'num{%.2f}' % (crab_g4v_hmin[m]/crab_g4v_hmin['G4v']))
t2 = tabulate(crabhmin_rel, headers=crabhmin_rel.keys(), tablefmt='latex')
t2 = t2.replace(r'\textbackslash{}', '\\').replace(r'\{', r'{').replace(r'\}',
                                                                        r'}')
# MP sensitivity
mpbest_table = {}
mprho_table = {}
crrho_table = {}
for d in detectors:
    for run in runs:
        t = tabulate(mpbest[d+run], headers='firstrow', tablefmt='latex')
        mpbest_table[d + run] = t.replace(r'\textbackslash{}', '\\').\
            replace(r'\{', r'{').replace(r'\}', r'}').replace(r'\$', '$').\
            replace(r'\_', '_')
	t3 = tabulate(mp_rho[d+run], headers='firstrow', tablefmt='latex')
	mprho_table[d + run] = t3.replace(r'\textbackslash{}', '\\').\
            replace(r'\{', r'{').replace(r'\}', r'}').replace(r'\$', '$').\
            replace(r'\_', '_')
	# crab
	t4 =tabulate(cr_rho[d+run], headers='firstrow', tablefmt='latex')
	crrho_table[d + run] = t4.replace(r'\textbackslash{}', '\\').\
            replace(r'\{', r'{').replace(r'\}', r'}').replace(r'\$', '$').\
            replace(r'\_', '_')


with open('tmp/plots/tables.txt', 'w') as f:
    f.write('%% Crab sensitivity\n')
    f.write(t1)
    f.write('\n\n')
    f.write('%% Crab relative sensitivity\n')
    f.write(t2)
    for d in detectors:
        for run in runs:
            f.write('\n\n')
            f.write('%%' + d + run + ' sensitivity\n')
            f.write(mpbest_table[d + run])
            f.write('\n\n')
            f.write('%%' + d + run + ' avg rho\n')
            f.write(mprho_table[d + run])
            f.write('\n\n')
            f.write('%%' + d + run + ' Crab rho\n')
            f.write(crrho_table[d + run])

