import logging
import cPickle as pickle
import sys
import os

import numpy as np
from scipy.interpolate import interp1d
from collections import OrderedDict

import h5py


# set up plotting backend
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = "stix"

from . import general as g

g.setuplog('results')

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


class Results(object):
    def __init__(self, det, run, psr, kind, methods=g.SEARCHMETHODS,
                 prefix='', suffix='', verbose=False):

        self.log = logging.getLogger('results.' + kind + '.' + psr)

        self.verbose = verbose

        self.det = det
        self.run = run
        self.psr = psr
        self.injkind = kind
        self.pair = None

        # Defining local srchmethods to avoid relying on lib.general.
        # This also allows one to load only results of a specific method.
        if isinstance(methods, basestring):
            # assuming only one method was requested, e.g. 'GR'
            self.search_methods = [methods]
        elif isinstance(methods, list):
            # assuming a list of methods was requested, e.g. ['GR', 'G4v']
            self.search_methods = methods
        else:
            self.log.warning('Search method: %r is not a string or a list'
                             % methods)

        self.paths = {
            'analysis': g.analysis_path(det, run, psr, kind),
            'export': prefix + 'results_' + det + run + '_' + psr + '_' +
                      kind + suffix + '.hdf5'
        }

        self.hinj = None
        self.ninst = None
        self.ninj = None

        # Initialize result containers
        # (DELIBERATELY NOT CONCISE TO SUPPORT NUMPY 2.6.6)
        self.hrec = {}  # recovered strength
        self.srec = {}  # recovery significance
        self.arec = {}  # recovery significance
        for m in self.search_methods:
            self.hrec[m] = []
            self.srec[m] = []
            self.arec[m] = []

        self.hob = {}  # ob strength
        self.sob = {}  # ob significance
        for m in self.search_methods:
            self.hob[m] = None
            self.sob[m] = None

    #--------------------------------------------------------------------------
    # IO operations
    def collect(self):
        """Collects results from multiple files in analysis_path/results as
        produced by submitting several injsrch_process jobs to Condor.
        NOT USED BY NEWEST VERSION.

        Note that "single" analyses do not rescale h when saving, so that is
        done here
        """

        self.log.info('Collecting results.')

        path = self.paths['analysis']

        try:
            with h5py.File(path + '/info.hdf5', 'r') as f:
                hinj = f['inj/h'][:]
                incinj = f['inj/inc'][:]
        except:
            raise IOError('Did not find injsrch info in: ' + path)

        # rescale hinj (note h_effective is independent of psi)
        hinj_new = []
        for h, i in zip(hinj, incinj):
            h = [h * ap(i, 0) for _, ap in
        g.Signal().templates[self.injkind].iteritems()]
            hinj_new.append(np.linalg.norm(h))

        self.hinj = np.array(hinj_new)

        self.ninst = len(self.hinj)

        # count nonzero elements in hinj:
        self.ninj = len(np.flatnonzero(self.hinj))

        self.log.debug('Looping over files.')

        search_methods = set()
        for n in range(self.ninst):
            self.log.debug('File ' + str(n))
            filename = path + '/results/r' + str(n) + '.p'

            try:
                with open(filename, 'rb') as f:
                    # load results dictionary
                    results = pickle.load(f)
                    # for each method retrieve h and s
                    for m in results.keys():
                        self.hrec[m].append(results[m]['h'])
                        self.srec[m].append(results[m]['s'])
                    search_methods = search_methods.union(set(results.keys()))
            except IOError:
                message = 'Unable to load result info from: ' + filename
                self.log.error(message, exc_info=True)
                print message
                print sys.exc_info()

        for m in search_methods:
            self.hrec[m] = np.array(self.hrec[m])
            self.srec[m] = np.array(self.srec[m])

        self.search_methods = list(search_methods)

    def export(self, path=None, verbose=None):
        """Exports results to HDf5 file.

        Arguments
        ---------
        path: `str` [optional]
            If no destination path is provided, takes cluster  public access
            folder as destination. This defaults to pwd when cluster is not
            identified.
        verbose: `boolean`
            Determines whether to print error message details.
        """

        self.log.info('Exporting results.')

        # Determine export destination.
        path = path or g.Cluster().public_dir
        export_path = path + self.paths['export']

        with h5py.File(export_path, 'w') as f:
            # save injected h
            f.create_dataset('hinj', data=self.hinj)
            # save recovered h and s
            for m in self.search_methods:
                grp = f.create_group(m)
                grp.create_dataset('hrec', data=self.hrec[m])
                grp.create_dataset('srec', data=self.srec[m])
                # grp.create_dataset('arec', data=self.arec[m])
        self.log.info('Results exported: ' + export_path)

    def load(self, path=None):
        """Loads results from HDf5 file.

        Arguments
        ---------
        path: `str` [optional]
            If no destination path is provided, takes cluster  public access
            folder as destination. This defaults to pwd when cluster is not
            identified.
        """
        self.log.info('Loading results.')
        # Determine origin destination.
        path = path or g.Cluster().public_dir
        export_path = path + self.paths['export']
        with h5py.File(export_path, 'r') as f:
            # load injected h
            self.hinj = f['hinj'][:]
            self.ninst = len(self.hinj)
            self.ninj = len(self.hinj[self.hinj > 0])
            # load recovered h and s
            for m in self.search_methods:
                self.hrec[m] = f[m + '/hrec'][:]
                self.srec[m] = f[m + '/srec'][:]

    def pickseries(self, kind, verbose=False):
        """
        Returns already-loaded recovered data of 'kind', making sure it exists.
        (Small snippet of code use multiple times below.)

        Arguments
        ---------
        kind: `basestring`
            Recovered parameter h or significance s.
        """
        if kind in ['s', 'srec', 'sig']:
            y = self.srec
            name = '$s$'
        elif kind in ['h', 'hrec', 'h0']:
            y = self.hrec
            name = r'$h_{\rm rec}$'
        else:
            self.log.error('Did not recognize value "' + str(kind) + '".',
                           exc_info=verbose)
            return
        return y, name

    #--------------------------------------------------------------------------
    # Statistics

    def get_detection_threshold(self, kind, threshold, methods=None):
        """
        Returns the value of hrec/srec which is above a threshold of the noise
        only instantiations.

        Arguments
        ---------
        kind: `str`
            Recovered parameter h or significance s.
        threshold: `float`
            Detection threshold.

        This percentage can be adapted by providing a different
        'threshold' value.
        Takes one argument that indicates whether the hrec or srec stat is
        computed.
        """
        methods = methods or self.search_methods
        d, _ = self.pickseries(kind)

        detection_threshold = {}
        for m in methods:
            # array of false positives
            noise_only = d[m][self.hinj == 0]
            n_noise_only = len(noise_only)
            # sort by loudness
            noise_only_sorted = np.sort(noise_only)
            # get threshold location
            t_loc = int(threshold * n_noise_only)
            # issue warning if threshold is above all false-positives
            if t_loc >= (n_noise_only - 1):
                self.log.warning('Threshold placed at loudest false positive.')
                t_loc = n_noise_only - 1
            # get threshold value
            detection_threshold[m] = noise_only_sorted[t_loc]
        # return dictionary of thresholds
        return detection_threshold

    def quantify(self, kind, detection_threshold, methods=None, band_conf=None,
                 det_conf=None, rollwindow=None):
        """
        Performs a linear fit ignoring values of hinj=0 and those under
        noise threshold.

        Returns
        -------
        m: `float`
            such that y = m*x from lstsq fit
        rmse: `float`
            sum of residuals = Sum[(x-y)^2], so RMSE = sqrt(s.o.r/N)
        ymax: `float`
            point defining line parallel to fit that encloses b_c% of
            pts OVER fit (inj, rec)
        ymin: `float`
            point defining line parallel to fit that encloses b_c% of
            pts UNDER fit (inj,rec)
        noise: `float`
            result of get_detection_threshold()
        """

        self.log.info('Performing linear fit of ' + str(kind) + ' data.')
        methods = methods or self.search_methods
        # obtain data
        d, _ = self.pickseries(kind)
        # obtain noise levels
        noise = self.get_detection_threshold(kind, detection_threshold,
                                             methods=methods)
        # get linear fit
        slope = {}
        rmse = {}
        ymax = {}
        ymin = {}
        y1side = {}
        x_short = {}
        rollmean = {}
        rollstd = {}
        for m in methods:
            self.log.debug('Selecting fit data.')
            # pick false postives above noise threshold
            x = self.hinj[(self.hinj != 0) & (d[m] > noise[m])]
            y = d[m][(self.hinj != 0) & (d[m] > noise[m])]
            # if not enough datapoints above noise to fit, raise error
            if len(y) < 2:
                raise RuntimeError('Not enough datapoints above noise.')
            # put vectors in proper shape for lstsq function
            x_vertical = np.reshape(x, (len(x), 1))
            y_vertical = np.reshape(y, (len(y), 1))

            self.log.debug('Performing fit.')
            # Fit datapoints to a straight line with y--intersect = 0
            # http://docs.scipy.org/doc/numpy/reference/generated/
            # numpy.linalg.lstsq.html
            p, residualsum, _, _ = np.linalg.lstsq(x_vertical, y_vertical)
            slope[m] = p[0][0]
            # the sum of residuals is the squared Euclidean 2-norm for each
            # column in b - a*x; i.e., sum of residuals = Sum[(x-y)^2],
            # so RMSE = sqrt(s.o.r./N)
            rmse[m] = np.sqrt(residualsum[0]/len(y))

            deviations = y - slope[m]*x

            if det_conf:
                # Get one sided confidence limit
                # (point to the left of b_c% of all points)
                oneside_id = int(np.floor((1-det_conf) * len(deviations)))
                y1side[m] = [(yi, yr) for (dev, yi, yr) in
                             sorted(zip(deviations, x, y))][oneside_id]

            if band_conf:
                self.log.debug('Computing bands at %.2f confidence.'
                               % band_conf)
                # Compute lines parallel to fit line that enclose `band_conf`
                # fraction of points.
                # 1. subtract fit from data and sort (`deviations` above)
                # 2a. find deviation above (band_conf %) of all points
                dev_top = deviations[deviations > 0]
                top_id = int(np.ceil(band_conf * len(dev_top)))
                if top_id > len(dev_top)-1:
                    self.log.warning('Top band above all points.')
                    top_id = -1
                # 2b. find corresponding (yinj, yrec)
                rec_top = y[deviations > 0]
                inj_top = x[deviations > 0]
                ymax[m] = [(yi, yr) for (dev, yi, yr) in
                           sorted(zip(dev_top, inj_top, rec_top))][top_id]
                # 3a. find value that is below (band_conf %) of all datapoints
                dev_bot = deviations[deviations <= 0]
                min_id = int(np.floor((1-band_conf) * len(dev_top)))
                # (note (1-band_conf) due to sorting of negative values.)
                if min_id > len(dev_top)-1:
                    self.log.warning('Bottom band below all points.')
                    min_id = -1
                # 3b. find corresponding (yinj, yrec)
                rec_bot = y[deviations <= 0]
                inj_bot = x[deviations <= 0]
                ymin[m] = [(yi, yr) for (dev, yi, yr) in
                           sorted(zip(dev_bot, inj_bot, rec_bot))][min_id]

            if rollwindow:
                x = self.hinj[self.hinj != 0]
                y = d[m][self.hinj != 0]
                self.log.debug('Computing rolling statistics.')
                # split data into h-bins of width window
                xsegments, ysegments = g.bindata(x, y, rollwindow)
                # take rolling mean of y (i.e. the mean of each segment)
                # (the if condition checks segment is not empty)
                rollmean[m] = np.array([ys.mean() for ys in ysegments
                                        if any(ys)]).flatten()
                # take rolling std of y (i.e. the std of each segment)
                rollstd[m] = np.array([ys.std() for ys in ysegments
                                       if any(ys)]).flatten()
                # get avg x location to plot
                x_short[m] = np.array([(xs[0] + xs[-1])/2. for xs in xsegments
                                       if any(xs)]).flatten()
        # produce output
        return noise, slope, rmse, ymax, ymin, y1side, x_short, rollmean,\
               rollstd

    def min_h_det(self, detection_threshold, confidence=.95):
        """Returns strength of smallest detected injection, computed using the
        significance curve and the corresponding detection.

        The significance threshold is translated into a strength by means
        of the fit.
        """
        self.log.info('Obtaining min h det.')
        # obtain significance noise threshold and best--fit line slope
        noise, slope, _, _, _, y1side, _, _, _ =\
            self.quantify('s', detection_threshold, det_conf=confidence)

        if confidence == 0:
            for m in self.search_methods:
                y1side[m] = (0, 0)

        minh = {}
        for m, n in noise.iteritems():
            minh[m] = (n - (y1side[m][1] - slope[m] * y1side[m][0])) / slope[m]

        return minh

    #--------------------------------------------------------------------------
    # Open box

    def openbox(self, methods=None, det_thrsh=None, det_conf=.95,
                p_nbins=100, p_fitorder=3, band_conf=None):

        methods = methods or [self.injkind, 'Sid']

        self.log.info('Opening box.')
        # try to load search results from file; otherwise, search
        h_ob = {}
        s_ob = {}
        for m in methods:
            if self.hob[m] is None or self.sob[m] is None:
                try:
                    filename = 'ob_%s%s_%s_%s' \
                               % (self.det, self.run, self.psr, m)
                    with open(g.paths['ob'] + filename + '.p', 'rb') as f:
                        results_ob = pickle.load(f)
                    self.log.debug('OB results loaded.')
                except IOError:
                    self.pair = self.pair or g.Pair(self.psr, self.det)
                    if self.pair.data is None:
                        self.pair.load_finehet(self.run, load_vectors=True)
                    results_ob = self.pair.search(methods=methods, save=True)
                self.hob[m] = results_ob[m]['h']
                self.sob[m] = results_ob[m]['s']
            h_ob[m] = self.hob[m]
            s_ob[m] = self.sob[m]
        if det_thrsh:
            hmin = self.min_h_det(det_thrsh, det_conf)
        if band_conf:
            _, slope, _, ymax, ymin, _, _, _, _ =\
                self.quantify('s', detection_threshold=det_thrsh,
                              band_conf=band_conf, methods=methods)
        p_ob = {}
        conf_int = {}
        hmin_dist = {}
        for m in methods:
            if p_fitorder:
                self.log.debug('Fitting p-value curve.')
                s_bk = self.srec[m][self.hinj == 0] # Background (bk) values
                s_pfit = g.pvalue(s_bk, p_nbins, fitorder=p_fitorder)[2]
                p_ob[m] = s_pfit(s_ob[m])

            if band_conf:
                self.log.debug('Forming confidence interval.')
                # get fit value
                h_ob_fit = s_ob[m] / slope[m]
                # get confidence limits
                nbot = ymin[m][1] - slope[m] * ymin[m][0]
                h_ob_max = (s_ob[m] - nbot) / slope[m]
                ntop = ymax[m][1] - slope[m] * ymax[m][0]
                h_ob_min = (s_ob[m] - ntop) / slope[m]
                # pack
                conf_int[m] = (float(h_ob_min), float(h_ob_fit),
                               float(h_ob_max))
            if det_thrsh:
                self.log.debug('Comparing to h_min.')
                hmin_dist[m] = (h_ob[m] - hmin[m]) / hmin[m]

        return h_ob, s_ob, p_ob, hmin_dist, conf_int

    #--------------------------------------------------------------------------
    # Plots

    def plot(self, kind, det_thrsh=.999, band_conf=.95, det_conf=None,
             methods=None, title=True, filetype='png', alpha=.3, shade=True,
             scale=1., hide_data=False, legend=True, xlim=None, ylim=None,
             rollwindow=2e-26, rollcolor=None, band=True, bestonly=False,
             suffix='', path=g.paths['plots'], bbox='tight', legendfont=14):
        """Plots 'kind' (hrec/srec) vs hinj for methods listed in 'methods'.
        """

        methods = methods or self.search_methods

        self.log.info('Plotting.')
        # obtain data
        y, kindname = self.pickseries(kind)
        if det_thrsh:
            # obtain fit & noise threshold
            noise, slope, rmse, ytop, ybot, y1side, x_short, rollmean, rollstd\
                = self.quantify(kind, det_thrsh, det_conf=det_conf,
                                band_conf=band_conf, methods=methods,
                                rollwindow=rollwindow)
            # find "best" method
            maxslope = np.max([slope[m] for m in methods])

        # process
        fig, ax = plt.subplots(1)
        for m in methods:
            # plot
            if not hide_data:
                ax.plot(self.hinj, y[m], g.plotcolor[m]+'+', label=m)

            if det_thrsh and slope[m] == maxslope:  # BUG HERE!
                # set axes limits
                ax.set_xlim(xlim or (0., scale * max(self.hinj)))
                if kind == 'h':
                    ylim = ylim or (0., np.max(self.hinj))
                    plt.gca().set_ylim(ylim[0], ylim[1])
                else:
                    ax.set_ylim(ylim or (0., scale * np.around(y[m].max(), 1)))

            #extra features
            switch = not bestonly or slope[m] == maxslope

            if det_thrsh and switch:
                # plot noise line
                noise_line = [noise[m]] * self.ninst
                ax.plot(self.hinj, noise_line, color=g.plotcolor[m])

            if det_conf and switch:
                det_line = slope[m] * self.hinj + (y1side[m][1] - slope[m] *
                                                   y1side[m][0])
                ax.plot(self.hinj, det_line, color=g.plotcolor[m],
                        linestyle='--')

            if band_conf is not None and band and switch:
                # plot band lines
                bestfit_line = slope[m] * self.hinj
                ax.plot(self.hinj, bestfit_line, color=g.plotcolor[m],
                        alpha=alpha)
                if band_conf:
                    topband_line = slope[m] * self.hinj + \
                                   (ytop[m][1] - slope[m] * ytop[m][0])
                    botband_line = slope[m] * self.hinj +\
                                   (ybot[m][1] - slope[m] * ybot[m][0])
                    ax.plot(self.hinj, topband_line,  color=g.plotcolor[m],
                            alpha=alpha)
                    ax.plot(self.hinj, botband_line,  color=g.plotcolor[m],
                            alpha=alpha)

                    if shade:
                        # shade confidence band
                        ax.fill_between(self.hinj, botband_line, topband_line,
                                        color=g.plotcolor[m], alpha=alpha/10,
                                        where=self.hinj > 0)
                        # note the where argument is necessary to close polygon

            if det_conf and rollcolor and switch:
                ax.plot(x_short[m], rollmean[m], rollcolor,
                        linewidth=2)
                ax.plot(x_short[m], rollmean[m] + rollstd[m],
                        rollcolor)
                ax.plot(x_short[m], rollmean[m] - rollstd[m],
                        rollcolor)

                if shade:
                    ax.fill_between(x_short[m], rollmean[m]-rollstd[m],
                                    rollmean[m] + rollstd[m],
                                    color=rollcolor, alpha=.3)

        # add labels indicating noise threshold and band confidence
        if det_thrsh:
            ax.text(.02, .7, 'Detection threshold: ' + str(det_thrsh),
                    fontsize=10, transform=ax.transAxes)
        if band_conf:
            ax.text(.02, .65, 'Band confidence: ' + str(band_conf),
                    fontsize=10, transform=ax.transAxes)
        # style
        ax.set_xlabel(r'$h_{\rm inj}$')
        ax.set_ylabel(kindname)
        if legend:
            ax.legend(numpoints=1, loc=2, prop={'size': legendfont})
        if title:
            ax.set_title('%s injections on %s %s data for %s'
                         % (self.injkind, self.det, self.run, self.psr))
        # check destination directory exists
        try:
            os.makedirs(path)
            self.log.debug('Plot directory created.')
        except OSError:
            self.log.debug('Plot directory already exists.')
        # save
        filename = 'injsrch_' + self.det + self.run + '_' + self.injkind +\
                   '_' + self.psr + '_' + kind
        p = path + filename + suffix + '.' + filetype
        fig.savefig(p, bbox_inches=bbox)
        plt.close(fig)
        print 'Figure saved: %r.' % p

    def plot_p(self, kind, methods=None, nbins=100, star=None, starsize=10,
               starcolor='y', fit=True, title=True, legend=True, legendloc=3,
               legendfont=14, xlim=False, ylim=(1e-4,1), path=g.paths['plots'],
               suffix='', filetype='png', hidedata=False, manyfiles=False,
               grid=False):

        self.log.info('Plotting 1-CDF')
        methods = methods or self.search_methods

        # obtain data
        d, kindname = self.pickseries(kind)
        suffix2 = ''
        if not manyfiles:
            fig, ax = plt.subplots(1)

        # process
        for m in methods:
            # p-value of falsepositives (y)
            x, y, pfit = g.pvalue(d[m][self.hinj == 0], nbins)
            # plot
            if manyfiles:
                fig, ax = plt.subplots(1)
                suffix2 = m

            if not hidedata:
                ax.plot(x, y, g.plotcolor[m]+'+', label=m)
            if fit:
                ax.plot(x, pfit(x), g.plotcolor[m], label=m+' fit')

            # format plots only if writing to many files or if all the plots
            # were already added to the figure
            if manyfiles or m == methods[-1]:
                # plot any extra points if they were provided
                if isinstance(star, dict):
                    # assuming a dictionary indexed by method was provided
                    # elements of dictionary should be a list/array
                    # significances/hrecs
                    s = star[m]
                    ax.plot(s, 0.1, starcolor + '*', markersize=starsize)
                elif isinstance(star, list) or isinstance(star, np.ndarray):
                    # assuming a list of significances was provided
                    for s in star:
                        ax.plot(s, 0.1, starcolor + '*', markersize=starsize)
                elif isinstance(star, (int, float)):
                    ax.plot(star, 0.1, starcolor + '*', markersize=starsize)

                # style
                ax.set_yscale('log')
                ax.set_xlabel(kindname)
                ax.set_ylabel('1-CDF')
                if xlim:
                    ax.set_xlim(xlim)
                if ylim:
                    ax.set_ylim(ylim)
                if legend:
                    ax.legend(numpoints=1, loc=legendloc,
                              prop={'size': legendfont})
                if title:
                    ax.set_title(self.injkind + ' ' +
                                 'injections on ' + self.det + ' ' + self.run
                                 + ' data for ' + self.psr)

                ax.grid(grid)

                try:
                    os.makedirs(path)
                    self.log.debug('Plot directory created.')
                except OSError:
                    self.log.debug('Plot directory already exists.')
                # save
                filename = 'pvalue_' + self.det + self.run + '_' +\
                           self.injkind + '_' + self.psr + '_'\
                           + kind
                saveto = '%s%s%s%s.%s'\
                         % (path, filename, suffix, suffix2, filetype)

                fig.savefig(saveto, bbox_inches='tight')
                plt.close(fig)
                print 'Plot saved: %r' % saveto

        # close all plots just in case
        plt.close('all')

    def plot_hs(self, methods=None, rollcolor=None, rollwindow=2e-26,
                title=True, xlim=None, ylim=None, hide_data=False, grid=False,
                shade=False, legend=False, legendloc=4, path=g.paths['plots'],
                filetype='png',  suffix='', manyfiles=True, bbox='tight'):
        """Plots 'kind' (hrec/srec) vs hinj for methods listed in 'methods'.
        """

        methods = methods or self.search_methods

        self.log.info('Plotting clean s vs h.')
        if not manyfiles:
            fig, ax = plt.subplots(1)
            plotname = ''

        for m in methods:
            ## 1. Plot clean background data
            self.log.debug('Obtaining data.')
            # 1a. gather data
            x = self.hrec[m][self.hinj == 0]
            y = self.srec[m][self.hinj == 0]
            # 1b. plot data
            if manyfiles:
                fig, ax = plt.subplots(1)
                plotname = '_search' + m
            if not hide_data:
                ax.plot(x, y, g.plotcolor[m]+'+', label=m)

            ## 2. Get rolling statistics
            if rollcolor:
                self.log.debug('Cumputing rolling statistics.')
                # split data into h-bins of width rollwindow
                xsegments, ysegments = g.bindata(x, y, rollwindow)
                # take rolling mean and std (i.e. the mean of each segment)
                # (the if condition checks segment not empty)
                rollmean = np.array([ys.mean() for ys in ysegments
                                     if any(ys)]).flatten()
                rollstd = np.array([ys.std() for ys in ysegments
                                    if any(ys)]).flatten()
                # get avg x location to plot
                x_short = np.array([(xs[0] + xs[-1])/2. for xs in xsegments
                                    if any(xs)]).flatten()
                # plot
                ax.plot(x_short, rollmean, rollcolor, linewidth=2)
                ax.plot(x_short, rollmean + rollstd, rollcolor)
                ax.plot(x_short, rollmean - rollstd, rollcolor)

                if shade:
                    ax.fill_between(x_short, rollmean - rollstd,
                                    rollmean + rollstd, color=rollcolor,
                                    alpha=.3)

            if manyfiles or m == methods[-1]:
                ## 4. Style
                ax.set_xlabel(r'$h_{\rm rec}$')
                ax.set_ylabel('s')
                if xlim:
                    ax.set_xlim(xlim)
                if ylim:
                    ax.set_ylim(ylim)
                if legend:
                    ax.legend(numpoints=1, loc=legendloc)
                if title:
                    ax.set_title('%s injections on %s %s data for %s'
                                 % (self.injkind, self.det,
                                    self.run, self.psr))
                ax.grid(grid)

                ## 5. Save
                try:
                    os.makedirs(path)
                    self.log.debug('Plot directory created.')
                except OSError:
                    self.log.debug('Plot directory already exists.')

                filename = 'hs_' + self.det + self.run + '_' + self.injkind\
                           + '_' + self.psr
                if manyfiles:
                    suffix2 = m
                else:
                    suffix2 = ''
                saveto = '%s%s%s%s.%s'\
                         % (path, filename, suffix, suffix2, filetype)
                fig.savefig(saveto, bbox_inches=bbox)
                plt.close(fig)
                print 'Plot saved: %r' % saveto
        # just in case, close all plots
        plt.close('all')


class ResultsMP(object):
    def __init__(self, injkind, det='H1', run='S5', path=None):
        self.det = det
        self.run = run
        self.injkind = injkind
        self.prefix = ''
        self.suffix = ''
        self.psrlist = []
        self._psrs = []
        self._results = OrderedDict()
        self.failed = []
        self._statkinds = ['hmin', 's_slope', 'h_slope', 's_noise',
                           'h_noise', 's_rmse', 'h_rmse']
        self.asd = None
        self.ts = None
        if path is not None:
            self.load(path=path)

    def _loadpsrs(self):
        for psr in self.psrlist:
            self._psrs.append(g.Pulsar(psr))

    #--------------------------------------------------------------------------
    # IO operations
    def load(self, path=None, prefix=None, suffix='',
             listid='all', verbose=False):
        """Load PSR results for all pulsars in list. Saves results objects to
        dictionary 'self.results'.
        """
        print 'Loading PSR results.'
        ### SETUP ###
        self.prefix = prefix or self.prefix
        self.suffix = prefix or self.suffix
        # Determine source file path:
        #   if a path was provided, use it;
        #   if not, create Cluster object and use its public d
        p = path or g.Cluster().public_dir
        # Load PSR lists
        psrlist = g.read_psrlist(name=listid, det=self.det, run=self.run)
        badpsrs = g.read_psrlist(name='bad', run=self.run)
        goodpsrs = set(psrlist) - set(badpsrs)
        ### PROCESS ###
        for psr in goodpsrs:
            try:
                r = Results(self.det, self.run, psr, self.injkind,
                            prefix=self.prefix, suffix=suffix)
                r.load(path=p)
                self._results[psr] = r
                # note that if there's an error loading, PSR object is NOT
                # added to list
            except IOError:
                self.failed.append(psr)
                print 'Warning: Unable to load ' + psr + ' results.'
                if verbose:
                    print sys.exc_info()
        self.psrlist = list(goodpsrs - set(self.failed))
        self._loadpsrs()

    def load_asd(self, path='/home/misi/Documents/P1200104-v4/'):
        if self.asd is None:
            self.log.info('ASD already loaded.')
        else:
            if self.run == 'S5':
                if self.det == 'H1':
                    self.ts = 527.*86400.
                    self.asd = np.loadtxt(path + 'lho4k_070318_strain.txt',
                                          comments='%')
                elif self.det == 'H2':
                    self.asd = np.loadtxt(path + 'lho2k_070514_strain.txt',
                                          comments='%')
                    self.ts = 535.*86400.
                elif self.det == 'L1':
                    self.ts = 405.*86400.
                    self.asd = np.loadtxt(path + 'llo_060604_strain.txt',
                                          comments='%')
                else:
                    print 'ERROR: invalid detector run combination: %r and %r'\
                          % (self.det, self.run)
                    return
            elif self.run == 'S6':
                if self.det == 'H1':
                    self.ts = 238.*86400.
                    self.asd = np.loadtxt(path + 'lho4k_15May2010_05hrs17min45'
                                                 'secUTC_strain.txt',
                                          comments='%')
                elif self.det == 'L1':
                    self.ts = 225.*86400.
                    self.asd = np.loadtxt(path + 'llo4k_31May2010_09hrs05min45s'
                                                 'ecUTC_strain.txt',
                                          comments='%')
                else:
                    print 'ERROR: invalid detector run combination: %r and %r'\
                          % (self.det, self.run)
                    return
            else:
                print "ERROR: %r injections are not supported." % self.run
                return
    #--------------------------------------------------------------------------
    # Statistics
    def getstat(self, kindstat, det_thrsh=.999, det_conf=.95, verbose=False):
        """Take all efficiency statistics for all PSRs loaded.
        """
        # values in case of failure
        default = {'noise': np.Inf, 'slope': None, 'rmse': np.Inf}
        # process
        output = {}
        for m in g.SEARCHMETHODS:
            output[m] = []
        if kindstat in ['hmin', 'h_min', 'minh']:
            for psr in self.psrlist:
                try:
                    r = self._results[psr]
                    hmin_psr = r.min_h_det(det_thrsh, confidence=det_conf)
                    for m in r.search_methods:
                        output[m].append(hmin_psr[m])
                except RuntimeError:
                    self.failed.append(psr)
                    for m in output.keys():
                        output[m].append(np.Inf)
                    print 'WARNING: bad statistics for PSR %s.' % psr
                    if verbose:
                        print sys.exc_info()
            for m in g.SEARCHMETHODS:
                output[m] = np.array(output[m])
        else:
            kind, statname = kindstat.split('_')
            for psr in self.psrlist:
                try:
                    r = self._results[psr]
                    stats = r.quantify(kind, det_thrsh)[0:3]
                    statsdict = {
                        'noise': stats[0],
                        'slope': stats[1],
                        'rmse': stats[2]
                    }
                    for m in r.search_methods:
                        output[m].append(statsdict[statname][m])
                except RuntimeError:
                    self.failed.append(psr)
                    for m in output.keys():
                        output[m].append(default[statname])
                    print 'WARNING: bad statistics for PSR %s.' % psr
                    if verbose:
                        print sys.exc_info()
            for m in g.SEARCHMETHODS:
                output[m] = np.array(output[m])
        return output

    def sortby(self, target, by, methods=None, det_thrsh=.999, det_conf=.95):
        """Returns instance of list of name 'target' sorted by 'by'.
        'target' can be 'psrlist', 'fgw' (=2*FR0) or the name of a PSR
        parameter (e.g. 'RAS')
        """
        methods = methods or g.SEARCHMETHODS
        # parse name of list to be sorted
        if target == 'psrlist':
            y = self.psrlist
        elif target == 'psrs':
            y = self._psrs
        elif 'gw' in target:
            # assume fgw requested
            ylist = [psr.param['FR0'] for psr in self._psrs]
            y = 2 * np.array(ylist).astype('float')
        else:
            ylist = [psr.param[target] for psr in self._psrs]
            y = np.array(ylist).astype('float')
        # parse name of list to be sorted BY
        if by in self._statkinds:
            x = self.getstat(by, det_thrsh=det_thrsh, det_conf=det_conf)
        else:
            try:
                x = [psr.param[by] for psr in self._psrs]
            except KeyError:
                print 'ERROR: unsupported `by` argument: %r' % by
                return
        ## PROCESS
        if isinstance(x, dict):
            y_sort = {}
            x_sort = {}
            for m in methods:
                y_sort[m] = np.array([yi for (xi, yi) in sorted(zip(x[m], y))])
                x_sort[m] = sorted(x[m])
        else:
            y_sort = np.array([yi for (xi, yi) in sorted(zip(x, y))])
            x_sort = sorted(x)
        return y_sort, x_sort

    def getpsrparam(self, paramname):
        paramlist = []
        for psr in self._psrs:
            paramlist.append(psr.param[paramname])
        return np.array(paramlist)

    def scalefactor(self):
        """
        Obtain ratio (\rho) between obtained hmin and expected sensitivity.
        """

    #--------------------------------------------------------------------------
    # Open boxes
    def openboxes(self, methods=g.SEARCHMETHODS, det_thrsh=None, det_conf=.95,
                  band_conf=None, p_fitorder=None):
        """Opens boxes for all PSRs in list and saves results.
        """
        print 'Opening boxes for all PSRs.'
        h_ob = {}
        s_ob = {}
        p_ob = {}
        h_conf = {}
        for m in methods:
            h_ob[m] = []
            s_ob[m] = []
            p_ob[m] = []
            h_conf[m] = []
        for psr, r in self._results.iteritems():
            try:
                h_ob_psr, s_ob_psr, p_ob_psr, _, h_conf_psr =\
                    r.openbox(methods=methods, det_thrsh=det_thrsh,
                              det_conf=det_conf, band_conf=band_conf,
                              p_fitorder=p_fitorder)
                for m in methods:
                    h_ob[m].append(h_ob_psr[m])
                    s_ob[m].append(s_ob_psr[m])
                    if p_fitorder:
                        p_ob[m].append(p_ob_psr[m])
                    if band_conf:
                        h_conf[m].append(h_conf_psr[m])
            except RuntimeError:
                print 'WARNING: bad statistics for PSR %s.' % psr
                for m in methods:
                    h_ob[m].append(0)
                    s_ob[m].append(0)
                    p_ob[m].append(1)
                    h_conf[m].append((0, 0, 0))
        for m in methods:
            h_ob[m] = np.array(h_ob[m])
            s_ob[m] = np.array(s_ob[m])
            p_ob[m] = np.array(p_ob[m])
            h_conf[m] = np.array(h_conf[m])
        return h_ob, s_ob, p_ob, h_conf

    def sortbyob(self, target, by, methods=None, det_thrsh=.999, det_conf=.95):
        """Returns instance of list of name 'target' sorted by 'by' (OB).
        'target' can be 'psrlist', 'fgw' (=2*FR0) or the name of a PSR
        parameter (e.g. 'RAS')
        """
        methods = methods or g.SEARCHMETHODS
        # parse name of list to be sorted
        if target == 'psrlist':
            y = self.psrlist
        elif target == 'psrs':
            y = self._psrs
        elif 'gw' in target:
            # assume fgw requested
            ylist = [psr.param['FR0'] for psr in self._psrs]
            y = 2 * np.array(ylist).astype('float')
        else:
            ylist = [psr.param[target] for psr in self._psrs]
            y = np.array(ylist).astype('float')
        # parse name of list to be sorted BY
        if by == 'h_ob':
            x = self.openboxes(det_thrsh=det_thrsh, det_conf=det_conf)[0]
        elif by == 's_ob':
            x = self.openboxes(det_thrsh=det_thrsh, det_conf=det_conf)[1]
        elif by == 'p_ob':
            x = self.openboxes(det_thrsh=det_thrsh, det_conf=det_conf,
                               p_fitorder=2)[2]
        else:
            print 'ERROR: unsupported `by` argument: %r' % by
            return
        ## PROCESS
        y_sort = {}
        x_sort = {}
        for m in methods:
            y_sort[m] = np.array([yi for (xi, yi) in sorted(zip(x[m], y))])
            x_sort[m] = sorted(x[m])
        return y_sort, x_sort

    #--------------------------------------------------------------------------
    # Plots
    def plot(self, psrparam, statkinds, det_thrsh=.999, det_conf=.95,
             suffix='', methods=None, path=g.paths['plots'], filetype='pdf',
             logy=False, title=True, legend=True, legend_loc='lower right',
             logx=False, xlim=None, ylim=None, grid=True, bbox='tight',
             legendfont=14):
        """Produces plot of efficiency indicator (noise, min-hrec) vs a PSR
        parameter (e.g. FR0, DEC, 'RAS').
        """
        methods = methods or [self.injkind]
        # obtain x-axis values
        # NOTE: astype(float) needed in case data was saved as string
        if 'gw' in psrparam:
            # plot GW frequency, not rotational frequency
            x = 2 * np.array([psr.param['FR0']
                              for psr in self._psrs]).astype(float)
        else:
            x = np.array([psr.param[psrparam]
                          for psr in self._psrs]).astype(float)

        if statkinds == 'all':
            kinds = ['h_slope', 'h_rmse', 'h_noise', 's_slope', 's_rmse',
                     's_noise', 'hmin']
        elif isinstance(statkinds, str):
            kinds = [statkinds]
        elif isinstance(statkinds, (list, tuple)):
            kinds = statkinds
        else:
            print 'ERROR: unrecognized statkind %r.' % statkinds
            return

        for kind in kinds:
            print 'Plotting %s vs. PSR %s' % (kind, psrparam)
            # obtain y-axis values
            y = self.getstat(kind, det_thrsh=det_thrsh, det_conf=det_conf)

            # create new figure for each stat kind
            fig, ax = plt.subplots(1)
            # plot all methods on same axes
            for m in methods:
                ax.plot(x, y[m], g.plotcolor[m]+'+', label=m)

            if xlim:
                ax.set_xlim(xlim[0], xlim[1])
            if ylim:
                ax.set_ylim(ylim[0], ylim[1])
            if grid:
                ax.grid()

            # Style
            if legend:
                plt.legend(loc=legend_loc, numpoints=1,
                           prop={'size': legendfont})
            if logy and kind != 's_noise':
                ax.set_yscale('log')
            if logx:
                ax.set_xscale('log')
            # (logscale if requested and if not sig noise threshold)
            # (allows for better formatting with 'all')
            if 'gw' in psrparam:
                ax.set_xlabel('GW Frequency [Hz]')
            else:
                ax.set_xlabel(psrparam)

            # parse kind
            if kind.startswith('h'):
                t = 'Strength'
                yl = r'$h_{\rm rec}$'
            elif kind.startswith('s'):
                t = 'Significance'
                yl = '$s$'

            if kind[2:] == 'slope':
                t += ' best fit slope'
                ylabel = 'Slope (' + yl + r' vs. $h_{\rm inj}$)'
            elif kind[2:] == 'rmse':
                t += ' best fit RMSE'
                ylabel = 'RMSE (' + yl + r' vs. $h_{\rm inj}$)'
            elif kind[2:] == 'noise':
                t += ' of noise threshold'
                ylabel = yl
            else:
                t = 'Lowest injection strength detected'
                ylabel = r'$h_{\rm min} (%.1f\%%)$' % (det_conf * 100.)
            t += '\nwith %.3f detection threshold at %.3f confidence'\
                 % (det_thrsh, det_conf)

            ax.set_ylabel(ylabel)
            if title:
                tt = '%s injections on %s %s %s data\n%s'\
                     % (self.injkind, self.det, self.run,
                        self.prefix, t)
                ax.set_title(tt)

            # save
            filename = 'mp_' + self.det + self.run + '_' + self.injkind +\
                       suffix + '_' + kind
            plt.savefig(path + self.prefix + filename + '.' + filetype,
                        bbox_inches=bbox)
            plt.close()
            print 'Plot saved: ' + path + self.prefix + filename + suffix + \
                  '.' + filetype

    def plot_ref(self, det_thrsh=.999, det_conf=.95, methods=None,
                 linecolor='gray', linewidth=1, linealpha=.5,
                 prefix='', suffix='', marker=None, markersize=10,
                 markercolor=None, legendloc=4, ylim=(1e-27, 1e-23),
                 xlim=(30, 1500), path=g.paths['plots'],
                 lp = '/home/misi/Documents/P1200104-v4/', filetype='pdf'):

        methods = methods or ['Sid', self.injkind]

        prefix = prefix or self.prefix
        suffix = suffix or self.suffix

        self.load_asd(path=lp)

        # get a range of frequencies
        fs = np.arange(xlim[0], xlim[1], 0.1)
        # interpolate amplitude spectral densities to same freq range
        intf = interp1d(self.asd[:, 0], self.asd[:, 1], kind='linear')
        # convert to power spectra
        sint = np.square(intf(fs))

        # scale factor for sensitivity estimate (Dupuis 2005)
        sf = 10.8
        # get the harmonic mean of the S5 time weighted sensitivities
        sens = sf * np.sqrt(sint/self.ts)  # ! MIGHT BE WRONG

        freq = 2 * np.array([psr.param['FR0']
                             for psr in self._psrs]).astype(float)
        hmin = self.getstat('hmin', det_thrsh=det_thrsh, det_conf=det_conf)

        # PLOT
        fig = plt.figure(figsize=(14, 10), dpi=300)

        plt.loglog(fs, sens, color=linecolor, alpha=linealpha,
                   linewidth=linewidth, label='Expected: %s %s'
                                              % (self.det, self.run))
        plt.hold(True)

        for m in methods:
            plt.loglog(freq, hmin[m], marker or g.plotmarker[m],
                       label="%s template" % m, markersize=markersize,
                       color=markercolor or g.plotcolor[m])

        plt.xlim(xlim[0], xlim[1])
        plt.ylim(ylim[0], ylim[1])
        plt.xlabel(r'Gravitational-wave Frequency (Hz)', fontsize=20,
                   fontweight=100)
        plt.ylabel(r'Strain Sensitivity %.2f (%i pulsars)' % (det_conf,
                                                              len(freq)),
                   fontsize=20, fontweight=100)

        plt.legend(numpoints=1, loc=legendloc)

        plt.setp(plt.gca().get_legend().get_texts(), fontsize='12')
        #fig.subplots_adjust(left=0.18, bottom=0.15)
        p = path + '%ssenscurve_%s%s_%s%s.%s' \
                   % (prefix, self.det, self.run, self.injkind, suffix,
                      filetype)
        fig.savefig(p)
        plt.close()
        print "Plot saved: %r" % p

    def plot_ob(self, psrparam, statkinds, det_thrsh=None, det_conf=None,
                band_conf=None, p_fitorder=None, methods=None,
                path=g.paths['plots'], filetype='pdf', log=False, title=True,
                legend=True, legend_loc='lower right', xlim=None, ylim=None,
                grid=True, fill=False):
        """ Produces plot of ob indicator (h_ob, s_ob, p_ob) vs a PSR parameter
        (e.g. FR0, DEC, 'RAS').
        """
        methods = methods or g.SEARCHMETHODS
        # parse x-axis name
        if 'gw' in psrparam:
            # plot GW frequency, not rotational frequency
            x = 2 * np.array([psr.param['FR0']
                              for psr in self._psrs]).astype(float)
        else:
            x = np.array([psr.param[psrparam]
                          for psr in self._psrs]).astype(float)
        # parse y-axis name
        if statkinds == 'all':
            kinds = ['h_ob', 's_ob', 'p_ob', 'h_conf']
        elif isinstance(statkinds, str):
            kinds = [statkinds]
        elif isinstance(statkinds, list):
            kinds = statkinds
        else:
            print 'ERROR: unrecognized statkind "'+str(statkinds)+'".'
            return
        # open boxes
        if 'p_ob' in kinds:
            # det_thrsh = det_thrsh or .999
            # det_conf = det_conf or .95
            # band_conf = band_conf or .95
            p_fitorder = p_fitorder or 2
        elif 'h_conf' in kinds:
            det_thrsh = det_thrsh or .999
            det_conf = det_conf or .95
            band_conf = band_conf or .95

        ob = self.openboxes(det_thrsh=det_thrsh, det_conf=det_conf,
                            band_conf=band_conf, p_fitorder=p_fitorder)
        ob_dict = {'h_ob': ob[0], 's_ob': ob[1],
                   'p_ob': ob[2], 'h_conf': ob[3]}
        # process
        for kind in kinds:
            print 'Plotting ' + kind + ' vs. PSR ' + psrparam
            # obtain y-axis values
            y = ob_dict[kind]
            # create new figure for each stat kind
            fig, ax = plt.subplots(1)
            # plot all methods on same axes
            for m in methods:
                if kind == 'h_conf':
                    xs = sorted(x)
                    ys = np.array([yi for (xi, yi) in sorted(zip(x, y[m]))])
                    plt.plot(xs, ys[:, 0], g.plotcolor[m], linestyle=':')
                    plt.plot(xs, ys[:, 1], g.plotcolor[m],
                             label=m, linewidth=1)
                    plt.plot(xs, ys[:, 2], g.plotcolor[m], linestyle=':')
                    if fill:
                        plt.fill_between(xs, ys[:, 0], ys[:, 2],
                                         color=g.plotcolor[m], alpha=.03)
                else:
                    plt.plot(x, y[m], g.plotcolor[m]+'+', label=m)
            # Style
            if legend:
                plt.legend(loc=legend_loc, numpoints=1, prop={'size': 14})
            if log:
                ax.set_yscale('log')
            if 'gw' in psrparam:
                ax.set_xlabel('GW Frequency [Hz]')
            else:
                ax.set_xlabel(psrparam)
            # parse kind
            if kind == 'h_ob':
                ylabel = r'$h_{\rm rec}$'
                t = 'Recovered strength'
            elif kind == 's_ob':
                ylabel = '$s$'
                t = 'Recovered significance'
            elif kind == 'p_ob':
                ylabel = '$p$'
                t = 'Detection $p$-value'
            elif kind == 'h_conf':
                ylabel = r'$h_{\rm rec}$'
                t = 'Recovered strength with confidence'
            else:
                print "WARNING: no label for kind %r" % kind
                ylabel = ''
            ax.set_ylabel(ylabel)
            if title:
                ax.set_title('Open boxes for ' + self.det + ' ' + self.run +
                             ' data ' + self.prefix + '\n' + t)
            if xlim:
                ax.set_xlim(xlim[0], xlim[1])
            if ylim:
                ax.set_ylim(ylim[0], ylim[1])
            if grid:
                ax.grid()
            # save
            filename = 'mpOB_' + self.det + self.run + '_' + self.injkind +\
                       '_' + kind
            plt.savefig(path + self.prefix + filename + '.' + filetype,
                        bbox_inches='tight')
            plt.close()
            print 'Plot saved: ' + path + self.prefix + filename + '.' +\
                  filetype
