import logging
from datetime import datetime
import numpy as np
import cPickle as pickle
import h5py
import sys
import math
import socket

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
        
        if not self.check_vectors(t, filename=''): self.create_vectors(t, filename='')

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
    def __init__(self, det, run, psr, kind, pdif):
        
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
        
        self.paths = {
                    'analysis' : analysis_path(det, run, psr, kind, pdif),
                    'export'   : 'results_' + det + run + '_' + psr + '_' + kind + pdif + '.hdf5'
                    }
                    
        # Initialize result containers (DELIBERATELY NOT CONCISE TO SUPPORT NUMPY 2.6.6)
        self.hrec = {} # recovered strength
        self.srec = {} # recovery significance
        for m in search_methods:
            self.hrec[m] = []
            self.srec[m] = []

    def collect(self):
        '''
        Collects results from multiple files in analysis_path/results as produced by
        submitting several injsrch_process jobs to Condor.
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
                    for m in search_methods:
                        self.hrec[m] += [results[m]['h']]
                        self.srec[m] += [results[m]['s']]
            except:
                message = 'Unable to load result info from: ' + filename
                self.log.error(message, exc_info=True)
                print message
                print sys.exc_info()
                # TEST: seeing if this makes 
                sys.exit(1)
        
        return self.hrec, self.srec
    
    def export(self, path=''):
        '''
        Exports collected result to hdf5 file. If no origin path is provided, takes cluster
        public access folder as destination.
        '''

        self.log.info('Exporting results.')
        
        if path=='':        
            # Determine export destination.
            cluster = Cluster()
            path = cluster.public_dir
        
        export_path = path + self.paths['export']
        
        try:
            with h5py.File(export_path, 'w') as f:
                
                # save injected h
                f.create_dataset('hinj', data = self.hinj)
                
                # save recovered h and s
                for m in search_methods:
                    grp = f.create_group(m)
                    grp.create_dataset('hrec', data = self.hrec[m])
                    grp.create_dataset('srec', data = self.srec[m])
        except:
            message = 'Unable to save collected results to: ' + export_path
            self.log.error(message, exc_info=True)
            print message
            print sys.exc_info()
            
    def load(self, path=''):
        '''
        Loads collected results from hdf5 file. If no origin path is provided, takes cluster
        public access folder as destination.
        '''
        self.log.info('Exporting results.')
        
        if path=='':        
            # Determine export destination.
            cluster = Cluster()
            path = cluster.public_dir
        
        export_path = path + self.paths['export']
        
        try:
            with h5py.File(export_path, 'r') as f:
                
                # save injected h
                self.hinj = f['hinj'][:]
                
                # save recovered h and s
                for m in search_methods:
                    self.hrec[m] = f[m + '/hrec'][:]
                    self.srec[m] = f[m + '/srec'][:]
        except:
            message = 'Unable to load collected results from: ' + export_path
            self.log.error(message, exc_info=True)
            print message
            
    def get_stats(self, methods=search_methods, types=['h','s']):
        
        self.log.info('Getting analysis statistics.')
        
        type_names = {'h' : 'hrec', 's' : 'srec'}
        
        # Initialize stats container: method/h or s/ stat kind
        self.stats = {}
        
        for m in methods:
            for type in types:
                for sk in stat_kinds:
                    self.stats[m] = {type : {sk: 0}}
                    # INCOMPLETE!
            
        
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
# STATISTICS

statkinds = [
            'slope',    # slope of linear fit to recovered injections (ignoring hinj=0)
            'noise',    # noise line
            'inter',    # hinj for which fit intersects noise line
            'rmse',     # root mean square error of lin. fit
            'min_det'   # minimum value above noise line
            ]

def lin_fit(x_in, y_in):
    '''
    Performs a linear fit: y = m x. Returns m. (Ignores values of x=0)
    '''
    
    ## Get injections
    
    # Mask zero-valued entries so that they are ignored
    # (See http://docs.scipy.org/doc/numpy/reference/maskedarray.generic.html)
    x_nozeros = np.ma.masked_equal(x_in,0)
    
    # get indices of x=0
    zero_index = np.flatnonzero(x_in)
    # get valid values of y
    y_nozeros = np.delete(y_in, zero_index)
    
    # Put vectors in proper shape
    x = np.reshape(x_nozeros, (len(x_nozeros), 1))
    y = np.reshape(y_nozeros, (len(y_nozeros), 1))
    
    # fit
    p, _, _, _ = np.linalg.lstsq(x, y)
    
    return np.poly1d([p[0][0], 0])
    
        
def noise_line(d):
    # input can be a dataframe
    # find elements with index 0 (no injection) and take their max over columns and rows
    try:
        max_n = d[0].max()
    except KeyError:
        max_n =d.T[0].max().max()
    # return corresponding constant polynomial
    return np.poly1d([max_n])

def rmse(d):
    # computes RMSE between d and its best-fit line

    # get injections
    d = d.drop([0])

    x = np.reshape(d.index, (len(d), 1)) # take injections
    y = np.reshape(d, (len(d), 1)) # take recovered values

    # fit and obtain sum of residuals (squared Euclidean 2-norm for each column in b - a*x)
    _, residuals, _, _ = np.linalg.lstsq(x, y)

    # sum of residuals = Sum[(x-y)^2], so RMSE = sqrt(s.o.r/N)
    return np.sqrt(residuals[0]/len(d))

def min_det_h(d):
    # assumes input is a series
    # get max noise
    noise = noise_line(d)(1)
    
    # find max injection below that
    d2 = d[d<noise]
    det_injs = d2[d2.index>0]
    
    try:
        print 'Get min'
        
        return det_injs.index[np.argmax(det_injs)]
    except ValueError:
	print 'ValueError'
        return 0

    
def fit_intersect_noise(d):
    # assumes input is a series
    
    noise = noise_line(d)       # noise line
    fit = lin_fit(d)            # fit line
    hdetmin = min_det_h(d)      # min detected h (probably around intersection)
    
    h_intersect = fsolve(lambda x: fit(x) - noise(x), hdetmin)
    
    return h_intersect[0]


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
