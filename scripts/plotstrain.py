import numpy as np
from scipy.interpolate import interp1d
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

det = 'H1'
run = 'S5'

linecolor='gray'
linewidth=1
linealpha=1
prefix=''
suffix=''
marker='*'
markersize=10
markercolor=None
legendloc = 1
filetype = 'png'
path = '/home/misi/polHTC/tmp/plots/'

lp = '/home/misi/Documents/P1200104-v4/'
# H1 S5
ts = 527.*86400.
asd = np.loadtxt(lp + 'lho4k_070318_strain.txt', comments='%')

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

fig = plt.figure(figsize=(14, 10), dpi=300)

plt.loglog(fs, intf(fs), color='b', alpha=linealpha,
           linewidth=linewidth, label='%s %s ASD'
                                      % (det, run))

plt.loglog(fs, sens, color=linecolor, alpha=linealpha,
           linewidth=linewidth, label='%s %s expected sens.'
                                      % (det, run))

plt.xlim(30, 1500)
#plt.ylim(1e-27, 1e-23)
plt.xlabel(r'Gravitational-wave Frequency (Hz)', fontsize=20,
           fontweight=100)
plt.ylabel(r'Strain Sensitivity ($1/\sqrt{{\rm Hz}}$)', fontsize=20, fontweight=100)

plt.legend(numpoints=1, loc=legendloc)

plt.setp(plt.gca().get_legend().get_texts(), fontsize='12')
#fig.subplots_adjust(left=0.18, bottom=0.15)
p = path + '%ssenscurve_%s%s%s.%s' \
           % (prefix, det, run, suffix, filetype)
fig.savefig(p)
plt.close()
print "Plot saved: %r" % p
