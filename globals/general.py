import logging
from datetime import datetime
import numpy as np
import cPickle as pickle
import h5py
import sys
import math
import socket
# set up plotting backend
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

##########################################################################################
## LOGGING

def setuplog(logname, logpath='logs/'):
    # Set up logging (from Logging Cookbook, Python online resources).
    # Name of log will be: lognamedDATETIME.log (spaces in logname replaced by '_')

    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                        datefmt='%m-%d %H:%M',
                        filename= logpath + logname.replace(' ','_') + '.log',
                        filemode='w')
            
    # define a Handler which writes INFO messages or higher to the sys.stderr
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)

    # set a format which is simpler for console use
    formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')

    # tell the handler to use this format
    console.setFormatter(formatter)

    # add the handler to the root logger
    logging.getLogger('').addHandler(console)


##########################################################################################
# CONSTANTS

ss = 86164.0905             # Seconds in a sidereal day

w = 2*np.pi/ss              # Sidereal angular frequency of Earth

rE = 6378.137e3             # Earth radius (m)


periodLIGO = 1/16384.       # LIGO data sampling period (s), from M. Pitkin

c = 299792458.              # Speed of light (m/s)

search_methods = ['GR', 'G4v', 'Sid'] #'AP', 'Sid'] ECONOMIC VERSION WITH JUST SID



##########################################################################################
# FUNCTIONS

def detarms(lat, lon, x_east, arm_ang, t):
    '''
    Computes detector arms for given detector parameters and time vector.
    '''
    # Angle between detector and Aries (vernal equinox) at time t
    # fiducial GPS time t0=630763213 (12hUT1 1/1/2000, JD245154).
    # See http://aa.usno.navy.mil/faq/docs/GAST.php
    offset = 67310.5484088 * w   # Aries-Greenwich angle at fiducial time (GMST)
    lmst = offset + w*(t-630763213) + lon # (LMST)
    
    # Earth center to North pole
    northPole = np.array([0, 0, 1])
    
    # The zenith is obtained using the detector location
    zenith = np.array([
        np.cos(lat)*np.cos(lmst),
        np.cos(lat)*np.sin(lmst),
        np.array([np.sin(lat)]*len(t))
        ]).transpose()  # [[x0, y0, z0], ...]

    # Local vectors are obtained from North pole and zenith
    localEast = np.cross(northPole, zenith)
    localNorth = np.cross(zenith, localEast)   
    
    # Rotating local vectors yields dx & dy. They are then normalized
    xArm = np.cos(x_east)*localEast + np.sin(x_east)*localNorth
    xArm /= np.sqrt(np.sum(xArm ** 2., axis=1))[..., None]

    # yArm is created from xArm using the angle between arms
    perp_xz = np.cross(zenith, xArm)
    yArm = xArm*np.cos(arm_ang) + perp_xz*np.sin(arm_ang) # equals perp_xz when angle between arms is 90deg
    yArm /= np.sqrt(np.sum(yArm ** 2., axis=1))[..., None]
    
    # scale dz to have norm Earth radius
    dz = rE * zenith
    
    return xArm, yArm, dz

def srcarms(dec, ras, psi):
    '''
    Computes source arms for given PSR parameters (DEC, RAS).
    '''
    north = np.array([0, 0, 1])

    # take source location vector components in celestial coordinates and invert direction multiplying by -1 to get wave vector wz
    wz = np.array([-math.cos(dec)*math.cos(ras), -math.cos(dec)*math.sin(ras), -math.sin(dec)])

    wy = np.cross(wz, north) 
    wy /= np.sqrt(np.sum(wy ** 2))

    wx = np.cross(wy, wz)
    wx /= np.sqrt(np.sum(wx ** 2))
    
    # rotate vectors by polarization angle
    wxRot = -wy*np.cos(psi) + wx*np.sin(psi)
    wyRot = wx*np.cos(psi) + wy*np.sin(psi) 
    
    return wxRot, wyRot, wz

def het(f, data, t):
    '''
    Returns data reheterodyned at frequency f for time t. Takes numpy arrays.
    '''
    omega = 2 * np.pi * f
    
    template = np.exp(1j * omega * t)
    return template * data



# INJECTION INFORMATION - used???
def phase2(pinj):
    if pinj == np.pi/2.:
        return 'p'
    elif pinj == -np.pi/2.:
        return 'm'
    elif pinj == 0:
        return '0'
    else:
        return str(pinj)
        
    
##########################################################################################
# OBJECTS

class Pulsar(object):
    def __init__(self, psrname):
        
        self.log = logging.getLogger('Pulsar')
        
        self.name = psrname
                
        self.log.info('Retrieving catalogue info.')
        try:
            psrcat = pickle.load(open(paths['psrcat'],'rb')) # psrcat should be a dict
            
            self.param = psrcat[self.name]
                
            self.log.debug('Info loaded.')
            
        except IOError:
            self.log.error('Could not find PSR catalogue in: ' + paths['psrcat'], exc_info=True)
        
        self.log.info('Computing source vectors.')
        
    def vectors(self, psi=None, save=False):
        '''
        Returns source vectors. Polarization angle psi can be provided (default psi=0).
        Vectors are not saved in class unless the save=True. Note that there is almost no
        computational advantage in saving these vectors: they are inexpensive to compute.
        '''
        
        if psi == None: psi = self.param['POL']
        
        wx, wy, wz = srcarms(self.param['DEC'], self.param['RAS'], psi)
        
        if save:
            self.wx = wx
            self.wy = wy
            self.wz = wz
            self.psi = psi
        
        return wx, wy, wz       


class Detector(object):
    def __init__(self, detname):
        self.log = logging.getLogger('Detector')

        self.name = detname
        self.observatory = detnames(detname)
        
        self.param = detparams[self.observatory]
        
        self.vecpath = paths['vectors'] + '/detvec' + self.name # path to vectors
        
    def create_vectors(self, t, filename=''):
        '''
        Returns arm vectors in Cartesian sidereal coordinates.
        '''
    
        self.log.info('Creating detector vectors.')
                    
        dx, dy, dz = detarms(
                            self.param['lat'],
                            self.param['lon'],
                            self.param['x_east'],
                            self.param['arm_ang'],
                            t
                            )
        
        self.log.info('Saving detector vectors.')
        
        try:
            f = h5py.File(self.vecpath + filename + '.hdf5', 'w')
            f.create_dataset('time', data=t)
            f.create_dataset('dx', data=dx)
            f.create_dataset('dy', data=dy)
            f.create_dataset('dz', data=dz)
            f.close()        
        except IOError:
            self.log.error('Unable to write det vecs to ' + self.vecpath + filename + '.hdf5')
        except:
            self.log.error('Unable to save det vecs to ' + self.vecpath + filename + '.hdf5', exc_info=True)

    def check_vectors(self, t, filename=''):
        self.log.info('Checking health of detector vector files.')
        try:
            with h5py.File(self.vecpath + filename + '.hdf5', 'r') as f:
                # make sure time series are the same            
                
                try:
                    if any(f['/time'][:]!=np.array(t)):
                        self.log.warning('Detector vectors do not agree with time series.')
                        return False
                    else:
                        self.log.debug('Detector vectors agree with time series.')
                        return True
                except:
                        self.log.warning('Error comparing time series.')
                        return False
        except IOError:
            self.log.warning('Did not find detector vectors in: ' + self.vecpath + filename + '.hdf5', exc_info=True)
            return False
        
    def load_vectors(self, t, filename=''):
        self.log.info('Loading detector vectors.')
        
        if not self.check_vectors(t, filename=filename): self.create_vectors(t, filename=filename)

        with h5py.File(self.vecpath + filename + '.hdf5', 'r') as f:
            self.dx = f['/dx'][:]
            self.dy = f['/dy'][:]
            self.dz = f['/dz'][:]
            
        self.log.debug('Detector vectors loaded.')


class Pair(object):
    '''
    Contains information of a PSR-det pair.
    '''
    
    def __init__(self, psrname, detname):
        self.psr = Pulsar(psrname)
        self.det = Detector(detname)
        
        try:
            self.log = logging.getLogger('Pair')
        except:
            setuplog('pair')
            self.log = logging.getLogger('Pair')

    def load_finehet(self, run, p='', check_vectors=False, load_vectors=False):
        
        self.run = run
        
        self.log.info('Checking finehet data is available.')
        
        if p=='':
            datapath = 'globals/data/' + self.det.name + '/' + run
            dataname = 'finehet_' + self.psr.name + '_' + self.det.name + '.hdf5'
            p = datapath + '/' + dataname

        try:
            finehet = h5py.File(p, 'r')
            
            self.time = finehet['/time'][:]
            self.data = finehet['/data'][:]
            
            self.log.debug('Finehet data loaded from: ' + p)
            
            # If check_vectors is True, check detector vectors exist for time in finehet.
            # If they do not, they are created.
            if check_vectors:
                self.log.debug('Checking vectors.')
                if not self.det.check_vectors(self.time, filename=self.psr.name): self.det.create_vectors(self.time, filename=self.psr.name)
            
            # If load_vectors is true, load detector vectors
            if load_vectors:
                self.log.debug('Loading vectors.')
                self.det.load_vectors(self.time, filename=self.psr.name)
            self.log.debug('Done')
            
        except:
            self.log.error('FATAL: No PSR data found in: ' + p,exc_info=True)
            print sys.exc_info()
    
    def get_sigma(self):
        '''
        Takes daily standard deviation, assuming finehet data has already been loaded.
        Returns sigma array and saves in object.
        '''
        # make sure time data has been loaded
        if 'time' not in dir(self):
            self.log.error('Cannot compute std: data not loaded.')
            sys.exit()
        
        t = self.time
        
        # find number of days which the data spans
        ndays = int(np.ceil((t[-1]-t[0])/ss))
        
        # this will also be the number of bins over which the data will be split
        count = np.histogram(t, bins=ndays)[0]
        
        
        # Note from histogram manual:
        # All but the last (righthand-most) bin is half-open. In other words, if bins is:
        # [1, 2, 3, 4]
        # then the first bin is [1, 2) (including 1, but excluding 2) and the second [2, 3).
        # The last bin, however, is [3, 4], which includes 4.
        
        # Split data series based on the count. Define initial and final indices for each
        # day, i_0 and i_f.
        i_0 = 0
        self.sigma = []
        for c in count:
            if c!=0:
                i_f = i_0 + c  # create final index

                segment = self.data[i_0:i_f] # pick day-worth of data

                self.sigma += [np.std(segment)] * c # take std and add to array

                i_0 = i_f # update initial index
        
        # See: http://www.uni-klu.ac.at/tewi/downloads/masterthesis_TrampitschStefan_Final_Version.pdf
        # For info about std of complex-valued data.
        
        return np.array(self.sigma)
          
    def signal(self, kind, pdif, pol, inc):
        
        self.log.info('Creating ' + kind + pdif + ' signal.')
        
        # Retrieve detector vectors.
        try:
            dx = self.det.dx
            dy = self.det.dy
                        
        except AttributeError:
            self.log.warning('No det vectors loaded. Attempting to load.', exc_info=True)
            
            self.det.load_vectors(self.time, filename=self.psr.name)
            
            dx = self.det.dx
            dy = self.det.dy
            
        # Retrieve source vectors.
        wx, wy, wz = self.psr.vectors(psi=pol)
        
        # Build signal. (See note under TEMPLATE INFORMATION below.)
        
        signal = np.zeros(len(self.time))+1j*np.zeros(len(self.time))
        
        for A, a in templateinfo[kind].iteritems():
            signal += a(inc, pdif) * (A(dx,dy,wx,wy,wz) + 0j)
            
        return signal
    
    def design_matrix(self, kind, pol, inc):
        # Retrieve detector vectors.
        if kind=='Sid':
            theta = w * self.time
            dm = [
                    np.ones(len(theta)),
                    np.cos(theta),
                    np.cos(2.*theta),
                    np.sin(theta),
                    np.sin(2.*theta)
                    ]
        else:
            try:
                dx = self.det.dx
                dy = self.det.dy
                        
            except AttributeError:
                self.log.warning('No det vectors loaded. Attempting to load.', exc_info=True)
            
                self.det.load_vectors(self.time, filename=self.psr.name)
            
                dx = self.det.dx
                dy = self.det.dy
            
            # Retrieve source vectors.
            wx, wy, wz = self.psr.vectors(psi=pol)
        
            # Build design matrix
            # NOTE: THERE'S NO SCALING OF h AT THIS STAGE!
        
            dm = []
            for A, a in templateinfo[kind].iteritems():
                dm += [a(inc, '0') * (A(dx,dy,wx,wy,wz) + 0j)]
        
        return np.array(dm)


class Results(object):
    def __init__(self, det, run, psr, kind, pdif, methods=search_methods):
        
        try:
            self.log = logging.getLogger('Results')
        except:
            setuplog('results_%(det)s_%(run)s_%(psr)s_%(kind)s%(pdif)s' % locals())
            self.log = logging.getLogger('Results')

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
                    'analysis' : analysis_path(det, run, psr, kind, pdif),
                    'export'   : 'results_' + det + run + '_' + psr + '_' + kind + pdif + '.hdf5'
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
            cluster = Cluster()
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
            self.log.error(message, exc_info=True)
            
    def load(self, path=''):
        '''
        Loads results from hdf5 file. If no origin path is provided, takes cluster
        public access folder as destination. This defaults to pwd when cluster is not identified.
        '''
        self.log.info('Exporting results.')
        
        if path=='':        
            # Determine origin destination.
            cluster = Cluster()
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
            self.log.error(message, exc_info=True)
            print message

    def pickseries(self, kind):
        '''
        Returns already loaded recovered data of 'kind'  and makes sure it exists.
        (Small snippet of code use multiple times below.)
        '''

        try:
            if kind in ['s', 'srec', 'sig']:
                y = self.srec
                name = 'Significance'

            elif kind in ['h', 'hrec', 'h0']:
                y = self.hrec
                name = '$h_{\rm rec}$'
                
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
            rmse[m] = np.sqrt(residualSum[0]/len(d))

            self.log.debug('Computing bands at ' + str(band_conf) + ' confidence.')

            # compute lines parallel to fit line that enclose band_conf (def 90%) of points
            
            # 1. subtract fit from data and sort
            deviations = np.sort(y - slope[m]*x)
            
            # 2a. find hrec/srec (y)  value that is above (band_conf %) of all datapoints
            dev_pos   = np.sort(deviations[deviations>0]) # sorted + deviations (yrec > yfit)
            ymax_rec_loc = int(band_conf * len(dev_pos)) # pick b_c% highest dev. index

            if ymax_rec_loc>=(len(dev_pos) - 1):
                self.log.warning('Threshold placed at loudest deviation.')
                ymax_rec_loc = len(dev_pos) - 1

            ymax_rec  = dev_pos[ymax_rec_loc] # pick b_c% highest value

            # 2b. find corresponding hinj value
            ymax_inj_loc = np.where(y==ymax_rec)[0]

            if len(ymax_inj_loc)!=1:
                # error caused if the value is found multiple times (len>1) or not at all (len=0)
                self.log.error('Cannot find ' + kind + ' max inj.', exc_info=True)
            print ymax_inj_loc, ymax_rec, y, len(x)
            ymax_inj = x[ymax_inj_loc[0]]

            ymax[m] = (ymax_inj, ymax_rec)

            # 3a. find value that is below (band_conf %) of all datapoints
            dev_neg   = np.sort(deviations[deviations<0])[::-1] # ([::-1] reverses list order)
            ymin_rec_loc = int(band_conf * len(dev_neg))
 
            if ymin_rec_loc>=(len(dev_neg) - 1):
                self.log.warning('Threshold placed at loudest deviation.')
                ymin_rec_loc = len(dev_neg) - 1

            ymin_rec  = dev_pos[ymin_rec_loc]
 
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

        return noise / slope
    
    #-----------------------------------------------------------------------------
    # Plots
    
    def plot(self, kind, aux='max', noise_threshold=.95, band_conf=.95, methods=[], dir='scratch/plots/', title=True):
         
        if methods==[]:
            methods = self.search_methods

        self.log.info('Plotting.')
 
        # obtain data
        y, kindname = self.pickseries(kind)

        # obtain fit & noise threshold
        slope, _, ymax, ymin, noise = self.quantify(kind, noise_threshold=noise_threshold,  band_conf= band_conf, methods=self.search_methods)
 
        # find "best" method
        maxslope = max([slope[m] for m in methods])
        
        # process
        for m in methods:
            # construct noise line, best fit line and confidence band around it
            noise_line = [noise[m]] * len(y)
            bestfit_line = slope[m] * self.hinj
            topband_line = slope[m] * self.hinj + (ymax[m][1]- slope[m] * ymax[m][0])
            botband_line = slope[m] * self.hinj + (ymin[m][1]- slope[m] * ymin[m][0])
            
            # plot
            plt.plot(self.hinj, y[m], plotcolor[m]+'+', label=m)

            if aux in ['all', 'full', 'simple']:
                plt.plot(self.hinj, bestfit_line, plotcolor[m])
                
                if aux in ['all', 'full']:
                    plt.plot(self.hinj, noise_line, plotcolor[m]+'.')
                    plt.plot(self.hinh, topband_line,  plotcolor[m], alpha=.5)
                    plt.plot(self.hinh, botband_line,  plotcolor[m], alpha=.5)
                    
                elif aux == 'simple':
                    # just plot the loudest noise threshold
                    if slope[m]==maxslope:
                        plt.plot(self.hinj, noise_line, plotcolor[m]+'.')
                        
        # style
        plt.xlim(0, max(self.hinj))
        plt.ylim(0, max(y))

        plt.xlabel('$h_{\rm inj}$')
        plt.ylabel(kindname)

        plt.legend(numpoints=1)

        if title: plt.title('%(self.injkind)s%(self.pdif)s injections on %(self.det)s %(self.run)s data for %(self.psr)s' % locals() )

        # check destination directory exists
        try:
            os.makedirs(dir)
            self.log.debug('Plot directory created.')
        except:
            self.log.debug('Plot directory already exists.')

        # save
        filename = 'injsrch_%(self.det)s%(self.run)s_%(self.injkind)s%(self.pdif)s_%(self.psr)s_%(kind)s' % locals()
        plt.savefig(dir + filename + '.pdf', bbox_inches='tight')
        
##########################################################################################
class Cluster(object):
    
    def __init__(self, name=''):
        
        # Set identity
        if name=='':
            # get hostname to determine what server we are on
            self.hostname = socket.gethostname()
        else:
            self.hostname = name
        
        # Determine scratch and public directories
        if 'ldas' in self.hostname:
            self.scratch_dir = '/usr1/max.isi/'
            self.public_dir  = '/home/max.isi/public_html/'
            
        elif 'atlas' in self.hostname:
            # get result of hostname -s command in bash
            hs = self.hostname.split('.')[0]
            #if hs == 'atlas3': hs = 'atlas3.atlas.aei.uni-hannover.de'
            self.scratch_dir = '/atlas/user/' + hs + '/max.isi/'
            self.public_dir  =  '/home/max.isi/WWW/LSC/'
            
        elif self.hostname in ['pccdev1', 'pcdev2', 'hydra','trout']:
            # assuming Nemo cluster
            self.scratch_dir = '/home/max.isi/scratch/'
            self.public_dir = '/home/max.isi/public_html/'
            
        else:
            self.scratch_dir = 'logs/'
            self.public_dir  = ''        

        
##########################################################################################
# TEMPLATE INFORMATION

# Polarization functions:

# - tensor
def pl(dx,dy,wx,wy,wz):
    # order matters because of broadcasting in numpy
    wxdx = np.dot(dx, wx)
    wxdy = np.dot(dy, wx)
    wydx = np.dot(dx, wy)
    wydy = np.dot(dy, wy)
    # Cross checked with PRD 85, 043005 (2012), PRD 79, 082002 (2009)
    return (wxdx**2 - wxdy**2 - wydx**2 + wydy**2)/2.
    
def cr(dx,dy,wx,wy,wz):
    wxdx = np.dot(dx, wx)
    wydx = np.dot(dx, wy)
    wxdy = np.dot(dy, wx)
    wydy = np.dot(dy, wy)
    # Cross checked with PRD 85, 043005 (2012), PRD 79, 082002 (2009)
    return wxdx * wydx - wxdy * wydy

# - vector
def xz(dx,dy,wx,wy,wz):
    wxdx = np.dot(dx, wx)
    wzdx = np.dot(dx, wz)
    wxdy = np.dot(dy, wx)
    wzdy = np.dot(dy, wz)
    # Cross checked with PRD 85, 043005 (2012), PRD 79, 082002 (2009)
    return wxdx * wzdx - wxdy * wzdy

def yz(dx,dy,wx,wy,wz):
    wydx = np.dot(dx, wy)
    wzdx = np.dot(dx, wz)
    wydy = np.dot(dy, wy)
    wzdy = np.dot(dy, wz)
    # Cross checked with PRD 85, 043005 (2012), PRD 79, 082002 (2009)
    return wydx * wzdx - wydy * wzdy

# - scalar
def br(dx,dy,wx,wy,wz):
    wxdx = np.dot(dx, wx)
    wxdy = np.dot(dy, wx)
    wydx = np.dot(dx, wy)
    wydy = np.dot(dy, wy)
    # Cross checked with PRD 85, 043005 (2012), PRD 79, 082002 (2009)
    return (wxdx**2 - wxdy**2 + wydx**2 - wydy**2)/2.

def lo(dx,dy,wx,wy,wz):
    wzdx = np.dot(dx, wz)
    wzdy = np.dot(dy, wz)
    # Cross checked with PRD 85, 043005 (2012), PRD 79, 082002 (2009)
    return np.sqrt(2)*(wzdx**2 - wzdy**2)/2.

# templateinfo dictionary includes an entry for each template. In turn, each entry contains
# a dictionary indexed by AP function and containing the weight associated to the polarization.
# NOTE: both the index (e.g. cr) and the content (e.g. np.cos(iota) ) are FUNCTIONS, not
# strings, so that they can be called directly without need for 'eval' or other methods.
 
templateinfo = { # (n = norm, p = phase)
        'GR'  : {
                pl : lambda iota, pdif : (1. + np.cos(iota)**2)/2. + 0j,
                cr : lambda iota, pdif : np.cos(iota) * np.exp(1j*pcat[pdif])
               },
                
        'G4v' : {
                xz : lambda iota, pdif : np.sin(iota) + 0j,
                yz : lambda iota, pdif : np.sin(iota)*math.cos(iota) * np.exp(1j*pcat[pdif])
                },
        'AP'  : {
                pl : lambda iota, pdif : 1. + 0j,
                cr : lambda iota, pdif : 1. + 0j,
                xz : lambda iota, pdif : 1. + 0j,
                yz : lambda iota, pdif : 1. + 0j,
                br : lambda iota, pdif : 1. + 0j
                }
               }

pcat = {
        'p' : np.pi/2.,
        'm' : - np.pi/2.,
        '0' : 0
        }

plotcolor = {
            'GR'  : 'g',
            'G4v' : 'r',
            'AP'  : 'b',
            'Sid' : 'm'
            }
##########################################################################################
# DETECTOR INFORMATION

ligoruns = ('S5', 'S6')
virgoruns = ('S1', 'S2')

detruns = {
            'H1' : ligoruns,
            'H2' : ligoruns,
            'L1' : ligoruns,
            'V1' : virgoruns,
            }
            
def detnames(d):
    if d in ['H1', 'H2']:
        det = 'LHO'
    elif d == 'L1':
        det = 'LHO'
    elif d == 'V1':
        det = 'VIR'
    elif d in['LHO', 'LLO', 'VIR']:
        det = d
    else:
        log.warning('g.datenames: This is not a valid det name: ' + str(d))
        sys.exit()
    return det
    
detparams = { # source for LIGO: LIGO-T980044; else: http://arxiv.org/abs/gr-qc/9607075)
        'LHO': {
                'lat': 0.81079526, # N 46deg 27'18.528''
                'lon': -2.084056769, # W 119deg 24'27.5657''
                'x_east': 2.199, # N 35.9994deg W
                'arm_ang': np.pi/2.
                },
    
        'LLO': {
                'lat': 0.533423135, # N 30deg 33'46.4196''
                'lon': -1.584309371, # W 90deg 46'27.2654''
                'x_east': 3.4508, # S 72.2835deg W
                'arm_ang': np.pi/2.
                },
    
        'GEO': {
                'lat': 0.9119345341670372,
                'lon': 0.17121679962064373,
                'x_east': 0.37716565135597474,
                'arm_ang': 1.646369083406251
                },
    
        'VIR': {
                'lat': 0.761487152645126,
                'lon': 0.1832595714594046,
                'x_east': 1.2479104151759457,
                'arm_ang': 1.5707963267948966
                },
    
        'TAM': {
                'lat': 0.6227334771115768,
                'lon': 2.4354324382328874,
                'x_east': 3.141592653589793,
                'arm_ang': 1.5707963267948966
                }
        }


##########################################################################################
# POLARIZATIONS

pol_sym = { 'pl' : '+', 'cr' : '\\times', 'xz' : 'x', 'yz' : 'y', 'br' : 'b', 'lo' : 'l'}

pol_names = {
            'pl' : 'plus',
            'cr' : 'cross',
            'xz' : 'vector x',
            'yz' : 'vector y',
            'br' : 'breathing',
            'lo' : 'longitudinal'}

pol_kinds = {
            'vector' : ['xz', 'yz'],
            'tensor' : ['pl', 'cr'],
            'scalar' : ['br', 'lo']
            }

pols = ['pl', 'cr', 'br', 'lo', 'xz', 'yz']


##########################################################################################
# PULSAR INFORMATION

paramNames = [
                '#',
                None,
                'RAS',
                'RAS error',
                'DEC',
                'DEC error',
            ]                

extraParamNames = [None, 'POL', 'POL error', 'INC', 'INC error']

extraParamStandard = {
                        'POL' : 0,
                        'POL error' : np.pi/4,
                        'INC' : 0,
                        'INC error' : np.pi/4
                        }

# formatting guide
paramFormat = {
                'RAS' : lambda x: hms_rad(x),
                'RAS error': lambda x: hms_rad(0., 0., x),
                'DEC' : lambda x: np.radians(dms_deg(x)),
                'DEC error' : lambda x: np.radians(dms_deg(0., 0., x))
                }

# read PSR list
def read_psrlist(name=''):
    log = logging.getLogger('psrlst')
    psrs = []
    
    if name=='':
        p = paths['psrlist']
    else:
        p = 'config/psrlist_' + name + '.txt'
    
    try:
        with open(p, 'r') as f:
            for line in f.readlines():
                psrs += [line.strip()] # (.strip() removes \n character)
        return psrs

    except:
        message = 'Could not open psrlist text in: ' + p
        log.error(message, exc_info=True)

# CONVERSIONS
def hmsformat(*args):
    if len(args[0]) == 1:
        # Assume hh:mm:ss format
        if type(args) != str:
            argument = args[0][0]
        else:
            argument = args
            
        hms = argument.split(':')
        if len(hms)==3:
            h, m, s = [float(x) for x in hms]
        elif len(hms)==2:
            m, s = [float(x) for x in hms]
            h = 0.0
        elif len(hms)==1:
            s = [float(x) for x in hms]
            h = 0.0
            m = 0.0
        else:
            print 'ERROR: hmsformat cannot convert: ' + argument
        
    elif len(args[0]) == 3:
        h, m, s = [float(x) for x in args[0]]

    else:
        print 'ERROR in hmsformat: can\'t take %d arguments' % len(args)    
    return h, m, s
    

def hms_rad(*args):
    # Converts hours, minutes, seconds to radians using the sidereal angular frequency of the Earth
    h, m, s = hmsformat(args)
    sec = s + 60*(m + 60*h)
    return sec*w
    
    
def dms_deg(*args):
    # Converts degrees, minutes, seconds to decimal degrees
    d, m, s = hmsformat(args)
    return d + m/60 + s/(60**2)
    
def masyr_rads(masyr):
    # Converts milliarcseconds/yr to radians/second
    asyr = masyr * 10 ** -3                     # mas/yr to arcseconds/yr 
    radyr = asyr * np.pi / 648000.              # as/yr to rad/yr (Wikipedia)
    rads = radyr / ss                           # rad/yr to rad/s
    return rads
    
def mjd_gps(mjd):
    # Converts MJD time to GPS time (taken from LALBarycenter.c line 749)
    tgps = 86400.*(mjd - 44244.) - 51.184
    return tgps





##########################################################################################
# PATHS

paths = {
            'extrapsrparam' : 'config/psrextra.txt',
            'psrcat' : 'globals/psrcat.p',
            'psrlist' : 'config/psrlist.txt',
            'badpsrs' : 'config/badpsrs.txt',
            'originalData' : '/home/matthew/analyses/S6_all/results',
            'vectors' : 'globals/vectors'
            }
            
def analysis_path(det, run, psr, kind, pdif):
    analysis = 'injsrch_' + det + run + '_' + psr + '_' + kind + pdif
    pathname = 'analyses/' + det + '/' + run + '/' + analysis
    return pathname
    
def submit_path(det, run, psr, kind, pdif, name='injsrch'):
    p = 'subs/%(name)s_%(det)s%(run)s_%(psr)s_%(kind)s%(pdif)s.sub' % locals()
    return p
    
def dag_path(det, run, psr, name='htc'):
    p = 'subs/%(name)s_%(det)s%(run)s_%(psr)s.dag' % locals()
    return p
    
localpaths = [
            'logs/',
            'results/'
            ]
