"""Generates MP reference plots, all in one for a given detector and run.
"""

import sys
import os
import numpy as np
from scipy.interpolate import interp1d
sys.path.append(os.getcwd())
from polHTC import results as r
# set up plotting backend
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = "stix"

mplparams = {
    'text.usetex': True,  # use LaTeX for all text
    'axes.linewidth': 1,  # set axes linewidths to 0.5
    'axes.grid': False,  # add a grid
    #'axes.axisbelow': True,
    #'grid.linestyle': '-',  # grid as solid lines not dashed
    #'grid.color': 'gray',  # grid lines grey
    #'grid.linewidth': 0.5,
    'font.family': 'serif',
    'font.size': 18
}
matplotlib.rcParams.update(mplparams)

det = sys.argv[1]
run = sys.argv[2]

linecolor='gray'
linewidth=1
linealpha=.5
prefix=''
suffix=''
marker='*'
markersize=10
markercolor=None
legendloc = 1
filetype = 'pdf'
path = '/home/misi/polHTC/tmp/plots/'

###############################################################################
p = '/home/misi/results/'
resultsGR = r.ResultsMP('GR', det=det, run=run, path=p)
resultsG4v = r.ResultsMP('G4v', det=det, run=run, path=p)

freqGR = 2 * np.array([psr.param['FR0']
                       for psr in resultsGR._psrs]).astype(float)
hminGR = resultsGR.getstat('hmin')

freqG4v = 2 * np.array([psr.param['FR0']
                       for psr in resultsG4v._psrs]).astype(float)
hminG4v = resultsG4v.getstat('hmin')


###############################################################################

lp = '/home/misi/Documents/P1200104-v4/'


if run == 'S5':
    if det == 'H1':
        ts = 527.*86400.
        asd = np.loadtxt(lp + 'lho4k_070318_strain.txt', comments='%')
    elif det == 'H2':
        asd = np.loadtxt(lp + 'lho2k_070514_strain.txt', comments='%')
        ts = 535.*86400.
    elif det == 'L1':
        ts = 405.*86400.
        asd = np.loadtxt(lp + 'llo_060604_strain.txt', comments='%')
    else:
        print 'ERROR: invalid detector run combination: %r and %r'\
              % (det, run)
        sys.exit(1)
elif run == 'S6':
    if det == 'H1':
        ts = 238.*86400.
        asd = np.loadtxt(lp + 'lho4k_15May2010_05hrs17min45secUTC'
                              '_strain.txt', comments='%')
    elif det == 'L1':
        ts = 225.*86400.
        asd = np.loadtxt(lp + 'llo4k_31May2010_09hrs05min45secUTC'
                              '_strain.txt', comments='%')
    else:
        print 'ERROR: invalid detector run combination: %r and %r'\
              % (det, run)
        sys.exit(1)
else:
    print "ERROR: %r injections are not supported." % run
    sys.exit(1)


# get a range of frequencies
fs = np.arange(30, 1500, 0.1)
# interpolate amplitude spectral densities to same freq range
intf = interp1d(asd[:, 0], asd[:, 1], kind='linear')
# convert to power spectra
sint = np.square(intf(fs))

# scale factor for sensitivity estimate (Dupuis 2005)
sf = 10.8
# get the harmonic mean of the S5 time weighted sensitivities
sens = sf * np.sqrt(sint/ts)  # ! MIGHT BE WRONG

###############################################################################

fig = plt.figure(figsize=(14, 10), dpi=300)

plt.loglog(freqGR, hminGR['Sid'], '^', label="GR injections (indep)" ,
                       markersize=markersize,
                       markeredgecolor='g', markerfacecolor='none')

plt.loglog(freqG4v, hminG4v['Sid'], '^', label="G4v injections (indep)" ,
                       markersize=markersize,
                       markeredgecolor='r', markerfacecolor='none')

plt.loglog(freqGR, hminGR['GR'], '*', label="GR injections (dep)" ,
                       markersize=markersize,
                       color='g')

plt.loglog(freqG4v, hminG4v['G4v'], '*', label="G4v injections (dep)" ,
                       markersize=markersize,
                       color='r')

plt.loglog(fs, sens, color=linecolor, alpha=linealpha,
           linewidth=linewidth, label='%s %s expected sens.'
                                      % (det, run))

plt.xlim(30, 1500)
plt.ylim(1e-27, 1e-23)
plt.xlabel(r'Gravitational-wave Frequency (Hz)', fontsize=20,
           fontweight=100)
plt.ylabel(r'Strain Sensitivity ($1/\sqrt{{\rm Hz}}$)', fontsize=20, fontweight=100)

plt.legend(numpoints=1, loc=legendloc)

plt.setp(plt.gca().get_legend().get_texts(), fontsize='12')
#fig.subplots_adjust(left=0.18, bottom=0.15)
p = path + '%ssenscurveINDP_%s%s%s.%s' \
           % (prefix, det, run, suffix, filetype)
fig.savefig(p)
plt.close()
print "Plot saved: %r" % p