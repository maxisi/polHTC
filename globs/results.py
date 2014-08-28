import logging
from datetime import datetime
import numpy as np
import scipy.stats
import cPickle as pickle
import h5py
import sys
import os
import math
import socket

# set up plotting backend
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = "stix"

from globs import general as g
#import general as g


g.setuplog('results')


class Results(object):
    def __init__(self, det, run, psr, kind, pdif, methods=g.search_methods, extra_name='', verbose=False):
        
        self.log = logging.getLogger('results.'+kind+'.'+psr)
            
        self.verbose = verbose

        self.det = det
        self.run = run
        self.psr = psr
        self.injkind = kind
        self.pdif = pdif

        # Defining local srchmethods to avoid relying on the global 'search_methods'
        # This also allows one to load only results of a specific method.
        if isinstance(methods, basestring):
            # assuming only one method was requested, e.g. 'GR'
            self.search_methods = [methods]
        elif isinstance(methods, list):
            # assuming a list of methods was requested, e.g. ['GR', 'G4v']
            self.search_methods = methods
        else:
            self.log.warning('Search method: "' + str(methods) + '" is not a string or a list')
        
        self.paths = {
                    'analysis' : g.analysis_path(det, run, psr, kind, pdif),
                    'export'   : extra_name + 'results_' + det + run + '_' + psr + '_' + kind + pdif + '.hdf5'
                    }
                    
        # Initialize result containers (DELIBERATELY NOT CONCISE TO SUPPORT NUMPY 2.6.6)
        self.hrec = {} # recovered strength
        self.srec = {} # recovery significance
        for m in self.search_methods:
            self.hrec[m] = []
            self.srec[m] = []
    
    #-----------------------------------------------------------------------------
    # IO operations

    def collect(self):
        '''
        Collects results from multiple files in analysis_path/results as produced by
        submitting several injsrch_process jobs to Condor. NOT USED BY NEWEST VERSION.
        '''
        
        self.log.info('Collecting results.')
        
        path = self.paths['analysis']
        
        try:
            with h5py.File(path + '/info.hdf5', 'r') as f:
                self.hinj = f['inj/h'][:]
        except:
            self.log.error('FATAL: did not find injsrch info in: ' + path, exc_info=True)

        self.ninst = len(self.hinj)
        self.ninj = len(np.flatnonzero(self.hinj)) # count nonzero elements in hinj
        
        self.log.debug('Looping over files.')
                
        for n in np.arange(0, self.ninst):
            self.log.debug('File ' + str(n))
            try:
                filename = path + '/results/r' + str(n) + '.p'
            
                with open(filename, 'rb') as f:
                
                    # load results dictionary
                    results = pickle.load(f)
                
                    # for each method retrieve h and s
                    for m in self.search_methods:
                        self.hrec[m] += [results[m]['h']]
                        self.srec[m] += [results[m]['s']]
            except:
                message = 'Unable to load result info from: ' + filename
                self.log.error(message, exc_info=True)
                print message
                print sys.exc_info()
        
        return self.hrec, self.srec
    
    def export(self, path=''):
        '''
        Exports results to hdf5 file. If no destination path is provided, takes cluster
        public access folder as destination. This defaults to pwd when cluster is not identified.
        '''

        self.log.info('Exporting results.')
        
        if path=='':        
            # Determine export destination.
            cluster = g.Cluster()
            path = cluster.public_dir
        
        try:
            export_path = '/home/max.isi/public_html/' + self.paths['export']

            with h5py.File(export_path, 'w') as f:

                # save injected h
                f.create_dataset('hinj', data = self.hinj)

                # save recovered h and s
                for m in self.search_methods:
                    grp = f.create_group(m)
                    grp.create_dataset('hrec', data = self.hrec[m])
                    grp.create_dataset('srec', data = self.srec[m])

        except IOError:
            export_path = path + self.paths['export']

            with h5py.File(export_path, 'w') as f:
                
                # save injected h
                f.create_dataset('hinj', data = self.hinj)
                
                # save recovered h and s
                for m in self.search_methods:
                    grp = f.create_group(m)
                    grp.create_dataset('hrec', data = self.hrec[m])
                    grp.create_dataset('srec', data = self.srec[m])
        except:
            message = 'Unable to save collected results to: ' + export_path
            self.log.error(message, exc_info=self.verbose)
            
    def load(self, path=''):
        '''
        Loads results from hdf5 file. If no origin path is provided, takes cluster
        public access folder as destination. This defaults to pwd when cluster is not identified.
        '''
        self.log.info('Loading results.')
        
        if path=='':        
            # Determine origin destination.
            cluster = g.Cluster()
            path = cluster.public_dir
        
        export_path = path + self.paths['export']
        
        try:
            with h5py.File(export_path, 'r') as f:
                
                # load injected h
                self.hinj = f['hinj'][:]
                
                # load recovered h and s
                for m in self.search_methods:
                    self.hrec[m] = f[m + '/hrec'][:]
                    self.srec[m] = f[m + '/srec'][:]
        except:
            message = 'Unable to load collected results from: ' + export_path
            self.log.error(message, exc_info=self.verbose)
            sys.exit(1)
            print message

    def pickseries(self, kind):
        '''
        Returns already loaded recovered data of 'kind'  and makes sure it exists.
        (Small snippet of code use multiple times below.)
        '''

        try:
            if kind in ['s', 'srec', 'sig']:
                y = self.srec
                name = '$s$'

            elif kind in ['h', 'hrec', 'h0']:
                y = self.hrec
                name = '$h_{rec}$'
                
            else:
                self.log.error('Did not recognize value "' + str(kind) + '".', exc_info=True)
                sys.exit(1)
            
            return y, name

        except:
            self.log.error('No data. Try loading results with .load()', exc_info=True)
            sys.exit(1)
    
    #-----------------------------------------------------------------------------
    # Statistics

    def get_noise_threshold(self, kind, methods=[], threshold=.9):
        '''
        Returns the value of hrec/srec which is above a threshold (def 90%) of the false
        positives. This percentage can be adapted by providing a different 'threshold' value.
        Takes one argument that indicates whether the hrec or srec stat is computed.
        '''

        if methods==[]:
            methods = self.search_methods
        
        d, _ = self.pickseries(kind)

        noise_threshold = {}
        for m in methods:
            # array of false positives
            false_pos = d[m][self.hinj==0]
            n_false_pos = len(false_pos)

            # sort by loudness
            false_pos_sorted = np.sort(false_pos)

            # get threshold location
            t_loc = int(threshold * n_false_pos)

            # issue warning if threshold is above all false-positives
            if t_loc>=(n_false_pos - 1):
                self.log.warning('Threshold placed at loudest false positive.')
                t_loc = n_false_pos - 1

            # get threshold value
            noise_threshold[m] = false_pos_sorted[t_loc]

        # return dictionary of thresholds
        return noise_threshold

    def quantify(self, kind, methods=[], noise_threshold=.9, band_conf=.9):
        '''
        Performs a linear fit ignoring values of hinj=0 and those under noise threshold.
        Returns:
        m       such that y = m*x from lstsq fit (float)
        rmse    sum of residuals = Sum[(x-y)^2], so RMSE = sqrt(s.o.r/N) (float)
        ymax    point defining line parallel to fit that encloses b_c% of pts OVER fit (inj, rec)
        ymin    point defining line parallel to fit that encloses b_c% of pts UNDER fit (inj,rec)
        noise   result of get_noise_threshold()
        '''

        self.log.info('Performing linear fit of ' + str(kind) + ' data.')
        
        if methods==[]:
            methods = self.search_methods

        # obtain data
        d, _  = self.pickseries(kind)

        # obtain noise levels
        noise = self.get_noise_threshold(kind, methods=methods, threshold=noise_threshold)

        # get linear fit
        slope = {}
        rmse  = {}
        ymax  = {}
        ymin  = {}
        for m in methods:
            
            self.log.debug('Selecting fit data.')

            # pick false postives above noise threshold
            x = self.hinj[(self.hinj!=0) & (d[m]>noise[m])]
            y = d[m][(self.hinj!=0) & (d[m]>noise[m])]

            # put vectors in proper shape for lstsq function
            x_vertical = np.reshape(x, (len(x), 1))
            y_vertical = np.reshape(y, (len(y), 1))

            self.log.debug('Performing fit.') 

            # fit using lstsq routine
            # http://docs.scipy.org/doc/numpy/reference/generated/numpy.linalg.lstsq.html
            p, residualSum, _, _ = np.linalg.lstsq(x_vertical, y_vertical)
            
            slope[m] = p[0][0]

            # the sum of residuals is the squared Euclidean 2-norm for each column in b - a*x
            # sum of residuals = Sum[(x-y)^2], so RMSE = sqrt(s.o.r./N)
            rmse[m] = np.sqrt(residualSum[0]/len(y))

            self.log.debug('Computing bands at ' + str(band_conf) + ' confidence.')

            # compute lines parallel to fit line that enclose band_conf (def 90%) of points
            
            # 1. subtract fit from data and sort
            deviations = y - slope[m]*x
            
            # 2a. find hrec/srec (y)  value that is above (band_conf %) of all datapoints
            dev_pos = deviations[deviations>0]
            
            max_ind = int(np.ceil(band_conf * len(dev_pos)))

            if max_ind>(len(dev_pos)-1):
                max_ind = -1
            
            dev_max = np.sort(dev_pos)[max_ind]
            
            ymax_rec = y[np.where(deviations==dev_max)[0][0]]
            # (note that many points might lay on the top/bottom lines, so where might return
            # multiple values; we just need one, so we take the first.)

            # 2b. find corresponding hinj value
            ymax_inj_loc = np.where(y==ymax_rec)[0]

            if len(ymax_inj_loc)!=1:
                # error caused if the value is found multiple times (len>1) or not at all (len=0)
                self.log.error('Cannot find ' + kind + ' max inj.', exc_info=True)
          
            ymax_inj = x[ymax_inj_loc[0]]

            ymax[m] = (ymax_inj, ymax_rec)

            # 3a. find value that is below (band_conf %) of all datapoints
            dev_neg = deviations[deviations<0]
            
            min_ind = int(np.floor((1-band_conf) * len(dev_pos)))
            # (note (1-band_conf) because of negative values.)
            
            if min_ind>(len(dev_pos)-1):
                min_ind = -1

            dev_min = np.sort(dev_neg)[min_ind]
            
            ymin_rec = y[np.where(deviations==dev_min)[0][0]]

            # 3b. find corresponding hinj value
            ymin_inj_loc = np.where(y==ymin_rec)[0]

            if len(ymin_inj_loc)!=1:
                # error caused if the value is found multiple times (len>1) or not at all (len=0)
                self.log.error('Cannot find ' + kind + ' min inj.', exc_info=True)
            
            ymin_inj = x[ymin_inj_loc]

            ymin[m] = (ymin_inj, ymin_rec)

        # return dictionary of thresholds and rmse's
        return slope, rmse, ymax, ymin, noise

    def min_h_det(self, confidence=.9):
        '''
        Returns strength of smallest detected injection, computed using the significance curve
        and the corresponding noise threshold. The significance threshold is translated into a
        strenght by means of the fit.
        '''

        self.log.info('Obtaining min h det.')

        # obtain significance noise threshold and best--fit line slope
        slope, _, _, _, noise = self.quantify('s', noise_threshold=confidence)
        
        minh = {}
        for m, n in noise.iteritems():
            minh[m] = n / slope[m]

        return minh
    
    #-----------------------------------------------------------------------------
    # Plots
    
    def plot(self, kind, aux='simple', noise_threshold=.99, band_conf=.95, methods=[], path='scratch/plots/', title=True, filetype='png', alpha=.3, shade=True, scale=1., extra_name='', hide_data=False, legend=True, setylim=True):
        '''
        Plots 'kind' (hrec/srec) vs hinj for methods listed in 'methods'.
        The argument 'aux' determines what extra features to include:
        -- 'full'/'all'
            For all methods adds: noise line above 'noise_threshold' (0-1) of the false positives
            best fit line, confidence band at 'band_conf' (0-1) confidence, shading if 'shade'.
        -- 'medium'
            For the "best" method: noise line above 'noise_threshold' (0-1) of false positives
            best fit line, confidence band at 'band_conf' (0-1) confidence, shading if 'shade'.
        -- 'simple'
            For the "best" method: noise line above 'noise_threshold' (0-1) of false positives
            best fit line.
        -- other
            Just the data points.
        '''
        
        if methods==[]:
            methods = self.search_methods

        self.log.info('Plotting.')
 
        # obtain data
        y, kindname = self.pickseries(kind)

        # obtain fit & noise threshold
        slope, _, ymax, ymin, noise = self.quantify(kind, noise_threshold=noise_threshold,  band_conf=band_conf, methods=methods)
 
        # find "best" method
        maxslope = np.max([slope[m] for m in methods])
        
        # process
        fig, ax = plt.subplots(1)
        for m in methods:
            # construct noise line, best fit line and confidence band around it
            noise_line = [noise[m]] * len(self.hinj)
            bestfit_line = slope[m] * self.hinj
            topband_line = slope[m] * self.hinj + (ymax[m][1]- slope[m] * ymax[m][0])
            botband_line = slope[m] * self.hinj + (ymin[m][1]- slope[m] * ymin[m][0])
            
            # plot
            if not hide_data:
                ax.plot(self.hinj, y[m], g.plotcolor[m]+'+', label=m)

            if aux in ['all', 'full', 'simple', 'medium']:
                # plot noise line
                ax.plot(self.hinj, bestfit_line, color=g.plotcolor[m])
                         
                if aux in ['all', 'full']:
                    # plot band lines
                    ax.plot(self.hinj, noise_line, color=g.plotcolor[m], alpha=alpha)
                    ax.plot(self.hinj, topband_line,  color=g.plotcolor[m], alpha=alpha)
                    ax.plot(self.hinj, botband_line,  color=g.plotcolor[m], alpha=alpha)

                    if shade:
                        # shade confidence band
                        ax.fill_between(self.hinj, botband_line, topband_line, color=g.plotcolor[m], alpha=alpha/10, where=self.hinj>0) # note the where argument is necessary to close polygon

                elif aux in ['simple', 'medium']:

                    # just plot the loudest noise threshold
                    if slope[m]==maxslope:
                        # plot noise line
                        ax.plot(self.hinj, noise_line, g.plotcolor[m], alpha=alpha)

                        if aux == 'medium':
                            # plot band lines
                            ax.plot(self.hinj, topband_line,  g.plotcolor[m], alpha=alpha)
                            ax.plot(self.hinj, botband_line,  g.plotcolor[m], alpha=alpha) 

                            if shade:
                                # shade confidence band
                                ax.fill_between(self.hinj, botband_line, topband_line, color=g.plotcolor[m], alpha=alpha/10, where=self.hinj>0)
            
            if slope[m]==maxslope: #BUG HERE!   
                # set axes limits
                ax.set_xlim(0., scale * max(self.hinj))
                if setylim: ax.set_ylim(0., scale * np.around(y[m].max(), 1)) # without the around the axes disapear when Sid
        
        # add labels indicating noise threshold and band confidence
        if aux in ['all', 'full', 'simple', 'medium']:
            ax.text(.02, .7, 'Noise threshold: ' + str(noise_threshold), fontsize=10, transform=ax.transAxes)

            if aux != 'simple':
                ax.text(.02, .65, 'Band confidence: ' + str(band_conf), fontsize=10, transform=ax.transAxes)
                        
        # style
        ax.set_xlabel('$h_{inj}$')
        ax.set_ylabel(kindname)

        if legend: ax.legend(numpoints=1, loc=2)

        if title: ax.set_title(self.injkind+self.pdif+' injections on '+ self.det+self.run+' data for '+self.psr)

        # check destination directory exists
        try:
            os.makedirs(path)
            self.log.debug('Plot directory created.')
        except:
            self.log.debug('Plot directory already exists.')

        # save
        filename = 'injsrch_'+self.det+self.run+'_'+self.injkind+self.pdif+'_'+self.psr+'_'+kind
        fig.savefig(path + filename + extra_name + '.' + filetype, bbox_inches='tight')
        plt.close(fig)
    
        
    def plot_p(self, kind, methods=[], nbins=100, star=None, starsize=6, starcolor='y', fit=True, title=True, legend=True, legendloc=3, xlim=False, ylim=(1e-4,1), path='scratch/plots/', extra_name='', filetype='png', manyfiles=False):
    
        if methods==[]: methods = self.search_methods

        self.log.info('Plotting 1-CDF')
        
        # obtain data
        y, kindname = self.pickseries(kind)
        
        if not manyfiles:
            fig, ax = plt.subplots(1)
            namelist = ''
        
        # process
        for m in methods:

            # p-value of falsepositives (y)
            x, y, pfit = g.pvalue(y[m][self.hinj==0], nbins)
            
            # plot
            if manyfiles:
                fig, ax = plt.subplots(1)
                namelist = m
            
            ax.plot(x, y, g.plotcolor[m]+'+', label=m)
            
            if fit: ax.plot(x, pfit(x), g.plotcolor[m], label = m + ' fit')
                
            # format plots only if writing to many files or if all the plots were already added to the figure
            if manyfiles or (not manyfiles and m==methods[-1]):

                # plot any extra points if they were provided
                if isinstance(star, dict):
                    # assuming a dictionary indexed by method was provided
                    # elements of dictionary should be a list/array significances/hrecs
                    for m in methods:
                        ax.plot(star[m][0], 0.1, starcolor + '*', markersize=starsize)
                        
                elif isinstance(star, list) or isinstance(star, np.ndarray):
                    # assuming a list of significances was provided
                    for s in star:
                        ax.plot(s, 0.1, starcolor + '*', markersize=starsize)
                        
                elif isinstance(star, float) or isinstance(star, int):
                    ax.plot(star, 0.1, starcolor + '*', markersize=starsize)
                
                # style
                
                ax.set_yscale('log')
                
                ax.set_xlabel(kindname)
                ax.set_ylabel('1-CDF')
                
                if xlim: ax.set_xlim(xlim)
                ax.set_ylim(ylim)

                if legend: ax.legend(numpoints=1, loc=legendloc)

                if title:
                    ax.set_title(self.injkind+self.pdif+' injections on '+ self.det+self.run+' data for '+self.psr)

                # check destination directory exists
                try:
                    os.makedirs(path)
                    self.log.debug('Plot directory created.')
                except:
                    self.log.debug('Plot directory already exists.')

                # save
                filename = 'pvalue_'+self.det+self.run+'_'+self.injkind+self.pdif+'_'+self.psr+'_'+kind
                saveto = path + filename + extra_name + namelist +'.' + filetype
                
                fig.savefig(saveto, bbox_inches='tight')
                
                plt.close(fig)
                
                self.log.info('Plot saved in ' + saveto)
            
        # close all plots just in case
        plt.close('all')
        
    def plot_hs(self, methods=[], stats=True, stats_color='y', window=2e-26, title=True, xlim=False, ylim=False, hide_data=False, shade=False, legend=False, legendloc=4, path='scratch/plots/', filetype='png',  extra_name='',  manyfiles=True):
        '''
        Plots 'kind' (hrec/srec) vs hinj for methods listed in 'methods'.
        '''
        
        if methods==[]: methods = self.search_methods

        self.log.info('Plotting clean s vs h.')
        
        if not manyfiles:
            fig, ax = plt.subplots(1)
            plotname = ''
 
        for m in methods:
            
            ## 1. Plot clean background data
            self.log.debug('Obtaining data.')
            
            # 1a. gather data
            x = self.hrec[m][self.hinj==0]
            y = self.srec[m][self.hinj==0]
            
            # 1b. plot data
            if manyfiles:
                fig, ax = plt.subplots(1)
                plotname = '_search' + m

            if not hide_data: ax.plot(x, y, g.plotcolor[m]+'+', label=m)
            
            ## 2. Get rolling statistics
            
            if stats:  
                self.log.debug('Cumputing rolling statistics.')  
                
                # splitt data into h-bins of width window
                xsegments, ysegments = g.bindata(x, y, window)
                
                # take rolling mean of y (i.e. the mean of each segment)
                # (the if condition checks segment not empty)
                rollmean = np.array([ys.mean() for ys in ysegments if any(ys)]).flatten()
                
                # take rolling std of y (i.e. the std of each segment)
                rollstd = np.array([ys.std() for ys in ysegments if any(ys)]).flatten()
                
                # get avg x location to plot
                x_short = np.array([(xs[0] + xs[-1])/2. for xs in xsegments if any(xs)]).flatten()
                
                # plot
                ax.plot(x_short, rollmean, stats_color, linewidth=2)
                ax.plot(x_short, rollmean + rollstd, stats_color)
                ax.plot(x_short, rollmean - rollstd, stats_color)
                
                if shade:
                    ax.fill_between(x_short, rollmean - rollstd, rollmean + rollstd, color=stats_color, alpha=.3)
            
            if manyfiles or (not manyfiles and m==methods[-1]):
                ## 4. Style
                
                ax.set_xlabel('$h_{rec}$')
                ax.set_ylabel('Significance')

                if xlim: ax.set_xlim(xlim)
                if ylim: ax.set_ylim(ylim)

                if legend: ax.legend(numpoints=1, loc=legendloc)

                if title:
                    ax.set_title(self.injkind+self.pdif+' injections on '+ self.det+self.run+' data for '+self.psr)
                
                ## 5. Save

                try:
                    os.makedirs(path)
                    self.log.debug('Plot directory created.')
                except:
                    self.log.debug('Plot directory already exists.')

                filename = 'hs_'+self.det+self.run+'_'+self.injkind+self.pdif+'_'+self.psr
                saveto = path + filename + extra_name + plotname +'.' + filetype
                
                fig.savefig(saveto, bbox_inches='tight')
                
                plt.close(fig)
                
                self.log.info('Plot saved in ' + saveto)
        
        # close all plots in case there was 
        plt.close('all')
    
    #-----------------------------------------------------------------------------
    # Open box
    
    def openbox(self):
        
        pair = Pair(self.psr, self.det)
        
        
    
        
class ResultsMP(object):
    def __init__(self, injkind, det='H1', run='S5', pdif='p'):
        self.det = det
        self.run = run
        self.injkind = injkind
        self.pdif = pdif
        self.noise_threshold = 0
        # (0 indicates stats haven't been taken yet)
            
    def load(self, path=None, extra_name='', listID='all', verbose=False):
        '''
        Load PSR results for all pulsars in list.
        Saves results objects to dictionary 'self.results'.
        '''
    
        print 'Loading PSR results.'
        
        ### SETUP ###
        
        # Determine source file path:
        #   if a path was provided, use it;
        #   if not, create Cluster object and use its public dir
        p = path or g.Cluster().public_dir

        self.extra_name = extra_name
        
        # Load PSR lists
        
        psrlist = g.read_psrlist(name=listID, det=self.det, run=self.run)
        badpsrs = g.read_psrlist(name='bad')
                    
        goodpsrs = set(psrlist) - set(badpsrs)
        
        ### PROCESS ###
        
        self.failed  = []
        self.results = {}
        
        for psr in goodpsrs:
            try:
                # create results object
                r = Results(self.det, self.run, psr, self.injkind, self.pdif, extra_name=extra_name)
                # load results
                r.load(path=p)
                
                self.results[psr] = r
                # note that if there's an error loading, PSR object is NOT added to list
                
            except:
                print 'Warning: Unable to load ' + psr + ' results.'
                if verbose: print sys.exc_info()

                self.failed += [psr]
                
        self.psrlist = list( goodpsrs - set(self.failed) )
    
    def get_stats(self, noise_threshold=.99, band_conf=.95, verbose=False):
        '''
        Take all efficiency statistics for all PSRs loaded.
        '''
        
        print 'Computing result statistics.'
        
        ### SETUP ###
        
        self.noise_threshold = noise_threshold
        
        # create stat containers
        
        for kind in ['h', 's']:
            for stat in ['slope', 'rmse', 'noise']:
                setattr(self, kind + '_' + stat, {})

        self.hmin = {}

        # (purposedly verbose to be compatible with python 2.6.6)
        for m in g.search_methods:
            self.hmin[m]    = []
            
            self.h_slope[m] = []
            self.h_rmse[m]  = []
            self.h_noise[m] = []

            self.s_slope[m] = []
            self.s_rmse[m]  = []
            self.s_noise[m] = []
                    
        self.psrs = []

        # load
        
        for psr in self.psrlist:
            try:
                # obtain results object
                r = self.results[psr]
                
                # get hrec noise and slope
                hslope, hrmse, _, _, hnoise = r.quantify('h', noise_threshold=noise_threshold, band_conf=band_conf)
                
                # get srec noise and slope
                sslope, srmse, _, _, snoise = r.quantify('s', noise_threshold=noise_threshold, band_conf=band_conf)
                
                # get min h detected
                hmin = r.min_h_det(confidence=noise_threshold)
                
                # save                
                for m in r.search_methods:
                    self.h_slope[m] += [hslope[m]]
                    self.h_rmse[m]  += [hrmse[m]]
                    self.h_noise[m] += [hnoise[m]]

                    self.s_slope[m] += [sslope[m]]
                    self.s_rmse[m]  += [srmse[m]]
                    self.s_noise[m] += [snoise[m]]
                    
                    self.hmin[m] += [hmin[m]]

                # get PSR data
                self.psrs += [g.Pulsar(psr)]
                
            except:
                print 'Warning: unable to get stats from ' + psr + ' results.'
                if verbose: print sys.exc_info()
                self.failed += [psr]
    
    def sortparams(self, target='psrlist', by='hmin', methods=None, noise_threshold=None):
        '''
        Returns instance of list of name 'name' sorted by hmin.
        'name' can be 'psrlist', 'fgw' (=2*FR0) or the name of a PSR parameter (e.g. 'RAS')
        '''
        
        ## SETUP
        
        nt = noise_threshold or self.noise_threshold
        
        # check stats are loaded and noise thresholds agree
        if (nt==0 or nt!=self.noise_threshold): self.get_stats()
        
        srchmethods = methods or self.hmin.keys()
        
        # parse name of list to be sorted
        if target == 'psrlist':
            y = self.psrlist
        elif 'gw' in target:
            y = 2 * np.array([psr.param['FR0'] for psr in self.psrs]).astype('float')
        else:
            y = np.array([psr.param[target] for psr in self.psrs]).astype('float')
            
        # parse name of list to be sorted BY
        x_dict = getattr(self, by)
        
        ## PROCESS
        y_sort = {}
        for m in srchmethods:
            x = self.hmin[m]
            y_sort[m] = [yi for (xi,yi) in sorted(zip(x,y))]

        return y_sort
                    
        
                
    def plot(self, statkinds, psrparam='fgw', extra_name='', scale=1., methods=g.search_methods, path='scratch/plots/', filetype='pdf', log=False, title=True, legend=True, legend_loc='lower right', xlim=None, grid=True):
        '''
        Produces plot of efficiency indicator (noise, min-hrec) vs a PSR parameter (e.g. FR0, DEC, 'RAS').
        '''
        
        try:
            self.hmin[methods[0]]
        except:
            self.get_stats()
       
        
        # obtain x-axis values
        if 'gw' in psrparam:
            # plot GW frequency, not rotational frequency
            x = 2 * np.array([psr.param['FR0'] for psr in self.psrs]).astype('float')
        else:
            x = np.array([psr.param[psrparam] for psr in self.psrs]).astype('float')
        
        if statkinds=='all':
            kinds = ['h_slope', 'h_rmse', 'h_noise', 's_slope', 's_rmse', 's_noise', 'hmin']
        elif isinstance(statkinds, str):
            kinds = [statkinds]
        elif isinstance(statkinds, list):
            kinds = statkinds
        else:
            print 'ERROR: unrecognized statkind "'+str(statkinds)+'".'
            
        for kind in kinds:
        
            print 'Plotting ' + kind + ' vs. PSR ' + psrparam
            
            # obtain y-axis values
            y = getattr(self, kind)
        
            # create new figure for each stat kind
            fig, ax = plt.subplots(1)
            
            # plot all methods on same axes
            for m in methods: plt.plot(x, y[m], g.plotcolor[m]+'+', label=m)
            
            # Style
            if legend: plt.legend(loc=legend_loc, numpoints=1)
            
            if log and kind!='s_noise': ax.set_yscale('log')
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
                t += ' best fit slope at ' + str(self.noise_threshold) + ' detection confidence'
                ylabel = 'Slope (' + yl + ' vs. $h_{inj}$)'
                
            elif kind[2:] == 'rmse':
                t += ' best fit RMSE at ' + str(self.noise_threshold) + ' detection confidence'
                ylabel = 'RMSE (' + yl + ' vs. $h_{inj}$)'
                
            elif kind[2:] == 'noise':
                t += ' of noise threshold at ' + str(self.noise_threshold) + ' detection confidence'
                ylabel = yl
                
            else:
                t = 'Lowest injection strength detected at ' + str(self.noise_threshold) + ' confidence'
                ylabel = '$h_{min}$'
                
            ax.set_ylabel(ylabel)
            
            if title: ax.set_title(self.injkind + self.pdif + ' injections on ' + self.det + ' ' + self.run + ' data ' + self.extra_name + '\n' + t)
            
            if xlim: ax.set_xlim(xlim[0], xlim[1])
            
            if grid: ax.grid()
            
            # save
            filename = 'mp_' + self.det + self.run + '_' + self.injkind + self.pdif + '_' + kind
            plt.savefig(path + self.extra_name + filename + '.' + filetype, bbox_inches='tight')
            plt.close()
            
            print 'Plot saved to ' + path + self.extra_name + filename + '.' + filetype
