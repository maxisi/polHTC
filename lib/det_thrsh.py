import logging
import cPickle as pickle
import sys
import os

import numpy as np
from collections import OrderedDict

import h5py


# set up plotting backend
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = "stix"

from . import general as g

g.setuplog('results')


class Results(object):
    def __init__(self, det, run, psr, kind, pdif, methods=g.SEARCHMETHODS,
                 extra_name='', verbose=False):

        self.log = logging.getLogger('results.'+kind+'.'+psr)

        self.verbose = verbose

        self.det = det
        self.run = run
        self.psr = psr
        self.injkind = kind
        self.pdif = pdif
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
            'analysis': g.analysis_path(det, run, psr, kind, pdif),
            'export': extra_name + 'results_' + det + run + '_' + psr + '_' +
                      kind + pdif + '.hdf5'
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

    #--------------------------------------------------------------------------
    # IO operations

    def collect(self):
        """
        Collects results from multiple files in analysis_path/results as
        produced by submitting several injsrch_process jobs to Condor.
        NOT USED BY NEWEST VERSION.
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
            h = [h * ap(i, self.pdif) for _, ap in
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
        """
        Exports results to hdf5 file. If no destination path is provided, takes
        cluster  public access folder as destination. This defaults to pwd when
        cluster is not identified.
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

    def load(self, path=None):
        """
        Loads results from hdf5 file. If no origin path is provided, takes
        cluster public access folder as destination. This defaults to pwd when
        cluster is not identified.
        """
        self.log.info('Loading results.')

        # Determine origin destination.
        path = path or g.Cluster().public_dir
        export_path = path + self.paths['export']

        try:
            with h5py.File(export_path, 'r') as f:
                # load injected h
                self.hinj = f['hinj'][:]
                self.ninst = len(self.hinj)
                self.ninj = len(self.hinj[self.hinj>0])
                # load recovered h and s
                for m in self.search_methods:
                    self.hrec[m] = f[m + '/hrec'][:]
                    self.srec[m] = f[m + '/srec'][:]
        except:
            message = 'Unable to load collected results from: ' + export_path
            self.log.error(message, exc_info=self.verbose)
            return

    def pickseries(self, kind):
        """
        Returns already-loaded recovered data of 'kind', making sure it exists.
        (Small snippet of code use multiple times below.)
        """

        try:
            if kind in ['s', 'srec', 'sig']:
                y = self.srec
                name = '$s$'
            elif kind in ['h', 'hrec', 'h0']:
                y = self.hrec
                name = '$h_{rec}$'
            else:
                self.log.error('Did not recognize value "' + str(kind) + '".', exc_info=True)
                return
            return y, name

        except:
            self.log.error('No data. Try loading results with .load()', exc_info=True)
            return

    #--------------------------------------------------------------------------
    # Statistics

    def get_detection_threshold(self, kind, threshold, methods=None):
        """
        Returns the value of hrec/srec which is above a threshold (def 99%)
        of the noise only instantiations.

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

    def quantify(self, kind, det_thrsh, methods=None, band_conf=None,
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
        self.log.info('Performing linear fit of %s data.' % kind)
        methods = methods or self.search_methods
        # obtain data
        d, _ = self.pickseries(kind)
        # obtain noise levels
        noise = self.get_detection_threshold(kind, det_thrsh, methods=methods)
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
                oneside_id = int(np.floor((1-band_conf) * len(deviations)))
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

    def min_h_det(self, det_thrsh, confidence=.95):
        """Returns strength of smallest detected injection, computed using the
        significance curve and the corresponding detection.

        The significance threshold is translated into a strength by means
        of the fit.
        """

        self.log.info('Obtaining min h det.')
        # obtain significance noise threshold and best--fit line slope
        noise, slope, _, _, _, y1side, _, _, _ =\
            self.quantify('s', det_thrsh, det_conf=confidence)

        if confidence == 0:
            for m in self.search_methods:
                y1side[m] = (0, 0)

        minh = {}
        for m, n in noise.iteritems():
            minh[m] = (n - (y1side[m][1] - slope[m] * y1side[m][0])) / slope[m]

        return minh

    #--------------------------------------------------------------------------
    # Open box

    def openbox(self, methods=None, det_thrsh=.999, det_conf=.95,
                p_nbins=100, p_fitorder=3, band_conf=.95):

        methods = methods or [self.injkind, 'Sid']

        self.log.info('Opening box.')
        # try to load search results from file; otherwise, search
        for m in methods:
            try:
                filename = 'ob_%s%s_%s_%s' % (self.det, self.run, self.psr, m)
                with open(g.paths['ob'] + filename + '.p', 'rb') as f:
                    results_ob = pickle.load(f)
                self.log.debug('OB results loaded.')
            except IOError:
                self.pair = self.pair or g.Pair(self.psr, self.det)
                if self.pair.data is None:
                    self.pair.load_finehet(self.run, load_vectors=True)
                results_ob = self.pair.search(methods=methods, save=True)

        self.log.debug('Obtaining stats.')
        hmin = self.min_h_det(det_thrsh, det_conf)

        _, slope, _, ymax, ymin, _, _, _, _ =\
            self.quantify('s', det_thrsh=det_thrsh,
                          band_conf=band_conf, methods=methods)

        h_ob = {}
        s_ob = {}
        conf_int = {}
        hmin_dist = {}
        p_ob = {}
        for m in methods:
            # Open box (ob) values
            h_ob[m] = results_ob[m]['h']
            s_ob[m] = results_ob[m]['s']
            # Background (bk) values
            s_bk = self.srec[m][self.hinj==0]

            self.log.debug('Forming confidence interval.')
            # get fit value
            h_ob_fit = s_ob[m] / slope[m]
            # get confidence limits
            nbot = (ymin[m][1]- slope[m] * ymin[m][0])
            h_ob_max = (h_ob_fit - nbot) / slope[m]
            ntop = (ymax[m][1]- slope[m] * ymax[m][0])
            h_ob_min = (h_ob_fit - ntop) / slope[m]
            # pack
            conf_int[m] = [float(h_ob_min), float(h_ob_fit), float(h_ob_max)]

            self.log.debug('Comparing to h_min.')
            hmin_dist[m] = (h_ob[m] - hmin[m]) / hmin[m]

            self.log.debug('Fitting p-value curve.')
            s_pfit = g.pvalue(s_bk, p_nbins, fitorder=p_fitorder)[2]
            p_ob[m] = s_pfit(s_ob[m])

        return h_ob, s_ob, p_ob, hmin_dist, conf_int

    #--------------------------------------------------------------------------
    # Plots

    def plot(self, kind, aux='simple', detection_threshold=.999, band_conf=.95,
             methods=None, title=True, filetype='png', alpha=.3, shade=True,
             scale=1., hide_data=False, legend=True, xlim=None, ylim=None,
             rollwindow=2e-26, rollcolor=None, band=True, bestonly=False,
             extra_name='', path=g.paths['plots']):
        """Plots 'kind' (hrec/srec) vs hinj for methods listed in 'methods'.
        """

        methods = methods or self.search_methods

        self.log.info('Plotting.')
        # obtain data
        y, kindname = self.pickseries(kind)
        if detection_threshold:
            # obtain fit & noise threshold
            noise, slope, rmse, ytop, ybot, y1side, x_short, rollmean, rollstd\
                = self.quantify(kind, detection_threshold, band_conf=band_conf,
                                methods=methods, rollwindow=rollwindow)
            # find "best" method
            maxslope = np.max([slope[m] for m in methods])

        # process
        fig, ax = plt.subplots(1)
        for m in methods:
            # plot
            if not hide_data:
                ax.plot(self.hinj, y[m], g.plotcolor[m]+'+', label=m)

            #extra features
            switch = band_conf and (not bestonly or slope[m] == maxslope)

            if detection_threshold and switch:
                # plot noise line
                noise_line = [noise[m]] * self.ninst
                ax.plot(self.hinj, noise_line, color=g.plotcolor[m])

            if band and switch:
                # plot band lines
                bestfit_line = slope[m] * self.hinj
                topband_line = slope[m] * self.hinj + (ytop[m][1]-slope[m] *
                                                       ytop[m][0])
                botband_line = slope[m] * self.hinj + (ybot[m][1]-slope[m] *
                                                       ybot[m][0])
                ax.plot(self.hinj, bestfit_line, color=g.plotcolor[m],
                        alpha=alpha)
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

            if rollcolor and switch:
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

            if detection_threshold and slope[m] == maxslope: #BUG HERE!
                # set axes limits
                ax.set_xlim(xlim or (0., scale * max(self.hinj)))
                ax.set_ylim(ylim or (0., scale * np.around(y[m].max(), 1)))

        # add labels indicating noise threshold and band confidence
        if detection_threshold:
            ax.text(.02, .7, 'Detection threshold: ' + str(detection_threshold),
                    fontsize=10, transform=ax.transAxes)
        if band_conf:
            ax.text(.02, .65, 'Band confidence: ' + str(band_conf),
                    fontsize=10, transform=ax.transAxes)
        # style
        ax.set_xlabel('$h_{inj}$')
        ax.set_ylabel(kindname)
        if legend:
            ax.legend(numpoints=1, loc=2)
        if title:
            ax.set_title('%s%s injections on %s %s data for %s'
                         % (self.injkind,self.pdif,self.det,self.run,self.psr))
        # check destination directory exists
        try:
            os.makedirs(path)
            self.log.debug('Plot directory created.')
        except OSError:
            self.log.debug('Plot directory already exists.')
        # save
        filename = 'injsrch_' + self.det + self.run + '_' + self.injkind +\
                   self.pdif + '_' + self.psr + '_' + kind
        p = path + filename + extra_name + '.' + filetype
        fig.savefig(p, bbox_inches='tight')
        plt.close(fig)
        print 'Figure saved: %r.' % p

    def plot_p(self, kind, methods=None, nbins=100, star=None, starsize=6,
               starcolor='y', fit=True, title=True, legend=True, legendloc=3,
               xlim=False, ylim=(1e-4,1), path=g.paths['plots'], extra_name='',
               filetype='png', hidedata=False, manyfiles=False):

        self.log.info('Plotting 1-CDF')
        methods = methods or self.search_methods

        # obtain data
        d, kindname = self.pickseries(kind)
        if not manyfiles:
            fig, ax = plt.subplots(1)
            plotname = ''

        # process
        for m in methods:
            # p-value of falsepositives (y)
            x, y, pfit = g.pvalue(d[m][self.hinj == 0], nbins)
            # plot
            if manyfiles:
                fig, ax = plt.subplots(1)
                plotname = m
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
                    s = star[m][0]
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
                    ax.legend(numpoints=1, loc=legendloc)
                if title:
                    ax.set_title(self.injkind + self.pdif + ' ' +
                                 'injections on ' + self.det + ' ' + self.run
                                 + ' data for ' + self.psr)

                try:
                    os.makedirs(path)
                    self.log.debug('Plot directory created.')
                except OSError:
                    self.log.debug('Plot directory already exists.')
                # save
                filename = 'pvalue_' + self.det + self.run + '_' +\
                           self.injkind + self.pdif + '_' + self.psr + '_'\
                           + kind
                saveto = '%s%s%s%s.%s'\
                         % (path, filename, extra_name, plotname, filetype)

                fig.savefig(saveto, bbox_inches='tight')
                plt.close(fig)
                print 'Plot saved: %r' % saveto

        # close all plots just in case
        plt.close('all')

    def plot_hs(self, methods=None, rollcolor=None, rollwindow=2e-26,
                title=True, xlim=None, ylim=None, hide_data=False,
                shade=False, legend=False, legendloc=4, path=g.paths['plots'],
                filetype='png',  extra_name='',  manyfiles=True):
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
                ax.set_xlabel('$h_{rec}$')
                ax.set_ylabel('s')
                if xlim:
                    ax.set_xlim(xlim)
                if ylim:
                    ax.set_ylim(ylim)
                if legend:
                    ax.legend(numpoints=1, loc=legendloc)
                if title:
                    ax.set_title('%s%s injections on %s %s data for %s'
                                 % (self.injkind, self.pdif, self.det,
                                    self.run, self.psr))

                ## 5. Save
                try:
                    os.makedirs(path)
                    self.log.debug('Plot directory created.')
                except OSError:
                    self.log.debug('Plot directory already exists.')

                filename = 'hs_' + self.det + self.run + '_' + self.injkind\
                           + self.pdif + '_' + self.psr
                saveto = '%s%s%s%s.%s'\
                         % (path, filename, extra_name, plotname, filetype)
                fig.savefig(saveto, bbox_inches='tight')
                plt.close(fig)
                print 'Plot saved: %r' % saveto
        # just in case, close all plots
        plt.close('all')


class ResultsMP(object):
    def __init__(self, injkind, det='H1', run='S5', pdif='p'):
        self.det = det
        self.run = run
        self.injkind = injkind
        self.pdif = pdif
        self.detection_threshold = None
        self.extra_name = ''
        self.psrlist = []
        self._psrs = []
        self._results = OrderedDict()
        self.failed = []
        # initialize h/s_slope/rmse/noise[m] = []
        for kind in ['h', 's']:
            for stat in ['slope', 'rmse', 'noise']:
                setattr(self, kind + '_' + stat, OrderedDict())
                for m in g.SEARCHMETHODS:
                    getattr(self, kind + '_' + stat)[m] = []
        # initialize hmin
        self.hmin = OrderedDict()
        for m in g.SEARCHMETHODS:
            self.hmin[m] = []
        # initialize h/s/p_ob
        for attr in ['h_ob', 's_ob', 'p_ob', 'h_conf']:
            setattr(self, attr, OrderedDict())

    #--------------------------------------------------------------------------
    # IO operations
    def load(self, path=None, extra_name=None, listid='all', verbose=False):
        """Load PSR results for all pulsars in list. Saves results objects to
        dictionary 'self.results'.
        """
        print 'Loading PSR results.'
        ### SETUP ###
        self.extra_name = extra_name or self.extra_name
        # Determine source file path:
        #   if a path was provided, use it;
        #   if not, create Cluster object and use its public d
        p = path or g.Cluster().public_dir
        # Load PSR lists
        psrlist = g.read_psrlist(name=listid, det=self.det, run=self.run)
        badpsrs = g.read_psrlist(name='bad')
        goodpsrs = set(psrlist) - set(badpsrs)

        ### PROCESS ###
        for psr in goodpsrs:
            try:
                r = Results(self.det, self.run, psr, self.injkind, self.pdif,
                            extra_name=self.extra_name)
                r.load(path=p)
                self._results[psr] = r
                # note that if there's an error loading, PSR object is NOT
                # added to list
            except IOError:
                print 'Warning: Unable to load ' + psr + ' results.'
                if verbose:
                    print sys.exc_info()
                self.failed += [psr]

        self.psrlist = list(goodpsrs - set(self.failed))

    #--------------------------------------------------------------------------
    # Statistics
    def get_stats(self, det_thrsh=.999, det_conf=.95, verbose=False):
        """Take all efficiency statistics for all PSRs loaded.
        """
        print 'Computing result statistics.'
        ### SETUP ###
        self.detection_threshold = det_thrsh
        for psr in self.psrlist:
            try:
                # obtain results object
                r = self._results[psr]
                # get s/hrec noise and slope
                hnoise, hslope, hrmse, _, _, _, _, _, _ = \
                    r.quantify('h', det_thrsh=det_thrsh)
                snoise, sslope, srmse, _, _, _, _, _, _ = \
                    r.quantify('s', det_thrsh=det_thrsh)
                # get min h detected
                hmin = r.min_h_det(det_thrsh, confidence=det_conf)
                # append corresponding values to h/s_slope/rmse/noise[m]
                for m in r.search_methods:
                    for stat in ['slope', 'rmse', 'noise']:
                        eval('self.h_' + stat + '[m].append(h'+stat + '[m])')
                        eval('self.s_' + stat + '[m].append(s'+stat + '[m])')
                    self.hmin[m].append(hmin[m])
                # get PSR data
                self._psrs.append(g.Pulsar(psr))
            except:
                print 'Warning: unable to get stats from ' + psr + ' results.'
                if verbose:
                    print sys.exc_info()
                self.failed.append(psr)

    def sortby(self, target, by, methods=None, det_thrsh=.999, det_conf=.95,
               band_conf=.95,):
        """Returns instance of list of name 'target' sorted by 'by'.
        'target' can be 'psrlist', 'fgw' (=2*FR0) or the name of a PSR
        parameter (e.g. 'RAS')
        """
        methods = methods or self.hmin.keys()
        # check stats are loaded and noise thresholds agree
        if det_thrsh != self.detection_threshold:
            self.get_stats(det_thrsh=det_thrsh, det_conf=det_conf)
        # parse name of list to be sorted
        if target == 'psrlist':
            y = self.psrlist
        elif 'gw' in target:
            # assume fgw requested
            ylist = [psr.param['FR0'] for psr in self._psrs]
            y = 2 * np.array(ylist).astype('float')
        else:
            ylist = [psr.param[target] for psr in self._psrs]
            y = np.array(ylist).astype('float')
        # parse name of list to be sorted BY
        if by in dir(self):
            x = getattr(self, by)
        else:
            try:
                x = [psr.param[by] for psr in self._psrs]
            except KeyError:
                print 'Unsupported `by` argument: %r' % by
                return
        ## PROCESS
        if isinstance(x, dict):
            y_sort = {}
            for m in methods:
                y_sort[m] = np.array([yi for (xi, yi)
                                      in sorted(zip(x[m], y))]).astype('float')
        else:
            y_sort = np.array([yi for (xi, yi)
                               in sorted(zip(x, y))]).astype('float')
        return y_sort

    #--------------------------------------------------------------------------
    # Open boxes
    def openboxes(self, methods=g.SEARCHMETHODS, det_thrsh=.999, det_conf=.95,
                  band_conf=.95):
        """Opens boxes for all PSRs in list and saves results.
        """
        print 'Opening boxes for all PSRs.'

        for attr in ['h_ob', 's_ob', 'p_ob', 'h_conf']:
            for m in methods:
                getattr(self, attr)[m] = []

        for psr, r in self._results.iteritems():
            ob = {}
            ob['h_ob'], ob['s_ob'], ob['p_ob'], _, ob['h_conf'] =\
                r.openbox(methods=methods, det_thrsh=det_thrsh,
                          det_conf=det_conf, band_conf=band_conf)

            for attr in ['h_ob', 's_ob', 'p_ob', 'h_conf']:
                for m in methods:
                    getattr(self, attr)[m].append(ob[attr][m])

    #--------------------------------------------------------------------------
    # Plots
    def plot(self, psrparam, statkinds, det_thrsh=.999, det_conf=.95,
             band_conf=.95, extra_name='', scale=1., methods=None,
             path=g.paths['plots'], filetype='pdf', log=False,
             title=True, legend=True, legend_loc='lower right', xlim=None,
             grid=True):
        """Produces plot of efficiency indicator (noise, min-hrec) vs a PSR
        parameter (e.g. FR0, DEC, 'RAS').
        """
        methods = methods or [self.injkind]
        # check if statistics are loaded
        if det_thrsh != self.detection_threshold:
            self.get_stats(det_thrsh=det_thrsh, det_conf=det_conf,
                           band_conf=band_conf)
        # obtain x-axis values
        if 'gw' in psrparam:
            # plot GW frequency, not rotational frequency
            x = 2 * np.array([psr.param['FR0'] for psr
                              in self._psrs]).astype('float')
        else:
            x = np.array([psr.param[psrparam] for psr
                          in self._psrs]).astype('float')

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
            y = getattr(self, kind)

            # create new figure for each stat kind
            fig, ax = plt.subplots(1)
            # plot all methods on same axes
            for m in methods:
                plt.plot(x, y[m], g.plotcolor[m]+'+', label=m)

            # Style
            if legend:
                plt.legend(loc=legend_loc, numpoints=1)
            if log and kind != 's_noise':
                ax.set_yscale('log')
            # (logscale if requested and if not sig noise threshold)
            # (allows for better formatting with 'all')
            if 'gw' in psrparam:
                ax.set_xlabel('GW Frequency [Hz]')
            else:
                ax.set_xlabel(psrparam)

            # parse kind
            if kind.startswith('h'):
                t = 'Strength'
                yl = '$h_{rec}$'
            elif kind.startswith('s'):
                t = 'Significance'
                yl = '$s$'

            if kind[2:] == 'slope':
                t += ' best fit slope with %.2f detection threshold at %.2f ' \
                     'confidence' % (det_thrsh, det_conf)
                ylabel = 'Slope (' + yl + ' vs. $h_{inj}$)'
            elif kind[2:] == 'rmse':
                t += ' best fit RMSE with %.2f detection threshold at %.2f ' \
                     'confidence' % (det_thrsh, det_conf)
                ylabel = 'RMSE (' + yl + ' vs. $h_{inj}$)'
            elif kind[2:] == 'noise':
                t += ' of noise threshold with %.2f detection threshold at %.2f ' \
                     'confidence' % (det_thrsh, det_conf)
                ylabel = yl
            else:
                t = 'Lowest injection strength detected above %.2f detection ' \
                    'threshold at %.2f confidence' % (det_thrsh, det_conf)
                ylabel = '$h_{min}$'

            ax.set_ylabel(ylabel)
            if title:
                tt = '%s%s injections on %s %s %s data\n%s'\
                     % (self.injkind, self.pdif, self.det, self.run,
                        self.extra_name, t)
                ax.set_title(tt)
            if xlim:
                ax.set_xlim(xlim[0], xlim[1])
            if grid:
                ax.grid()

            # save
            filename = 'mp_' + self.det + self.run + '_' + self.injkind +\
                       self.pdif + '_' + kind
            plt.savefig(path + self.extra_name + filename + '.' + filetype,
                        bbox_inches='tight')
            plt.close()
            print 'Plot saved: ' + path + self.extra_name + filename + '.' +\
                  filetype

    def plot_ob(self, statkinds, psrparam='fgw', extra_name='', scale=1., methods=None, path='scratch/plots/', filetype='pdf', log=False, title=True, legend=True, legend_loc='lower right', xlim=None, grid=True):
        """ Produces plot of ob indicator (h_ob, s_ob, p_ob) vs a PSR parameter
        (e.g. FR0, DEC, 'RAS').
        """

        # check if OB results already loaded
        try:
            [self.h_ob[methods[m]] for m in methods]
        except:
            self.openboxes(methods=methods)

        # parse x-axis name
        if 'gw' in psrparam:
            # plot GW frequency, not rotational frequency
            x = 2 * np.array([psr.param['FR0'] for psr in self.psrs]).astype('float')
        else:
            x = np.array([psr.param[psrparam] for psr in self.psrs]).astype('float')

        # parse y-axis name
        if statkinds=='all':
            kinds = ['h_ob', 's_ob', 'p_ob', 'h_conf']

        elif isinstance(statkinds, str):
            kinds = [statkinds]

        elif isinstance(statkinds, list):
            kinds = statkinds

        else:
            print 'ERROR: unrecognized statkind "'+str(statkinds)+'".'

        # process
        for kind in kinds:

            print 'Plotting ' + kind + ' vs. PSR ' + psrparam

            # obtain y-axis values
            y = getattr(self, kind)

            # create new figure for each stat kind
            fig, ax = plt.subplots(1)

            # plot all methods on same axes
            for m in methods:
                if kind=='h_conf':
                    plt.plot(x, y[m][:][0], g.plotcolor[m]+'+') # h_min
                    plt.plot(x, y[m][:][1], g.plotcolor[m]+'+', label=m, linewidth=2) # h_fit
                    plt.plot(x, y[m][:][2], g.plotcolor[m]+'+') # h_max
                else:
                    plt.plot(x, y[m], g.plotcolor[m]+'+', label=m)

            # Style
            if legend: plt.legend(loc=legend_loc, numpoints=1)

            if log: ax.set_yscale('log')

            if 'gw' in psrparam:
                ax.set_xlabel('GW Frequency [Hz]')
            else:
                ax.set_xlabel(psrparam)

            # parse kind
            if kind=='h_ob':
                ylabel = '$h_{rec}$'
                t = 'Recovered strength'
            elif kind=='s_ob':
                ylabel = '$s$'
                t = 'Recovered significance'
            elif kind=='p_ob':
                ylabel = '$p$'
                t = 'Detection $p$-value'
            elif kind=='h_conf':
                ylabel = '$h_{rec}$'
                t = 'Recovered strength with confidence'

            ax.set_ylabel(ylabel)

            if title: ax.set_title('Open boxes for ' + self.det + ' ' + self.run + ' data ' + self.extra_name + '\n' + t)

            if xlim: ax.set_xlim(xlim[0], xlim[1])

            if grid: ax.grid()

            # save
            filename = 'mpOB_' + self.det + self.run + '_' + self.injkind + self.pdif + '_' + kind
            plt.savefig(path + self.extra_name + filename + '.' + filetype, bbox_inches='tight')
            plt.close()

            print 'Plot saved to ' + path + self.extra_name + filename + '.' + filetype
