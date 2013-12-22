import warnings

def formatwarning(message, category, filename, lineno, line = 0):
    strout = "\n***********\n%s:%s: %s:\n\t%s\n***********\n" % (filename, lineno, category.__name__, message)
    return strout

warnings.formatwarning = formatwarning

warn = warnings.warn
from collections import defaultdict

from netCDF4 import Dataset
from numpy import exp, log, pi, zeros
import sys

RHOSO4  = 1.8e3 # bulk density of aerosol sulfate
RHONH4  = 1.8e3 # bulk density of aerosol ammonium
RHONO3  = 1.8e3 # bulk density of aerosol nitrate
RHOORG  = 1.3e3 # bulk density for aerosol organics following Carlton et al. 2010
RHOSOIL = 2.6e3 # bulk density for aerosol soil dust
RHOSEAS = 2.2e3 # bulk density for marine aerosol
RHOANTH = 2.2e3 # bulk density for anthropogenic aerosol
SGINIAT = 1.7   # initial sigma-G for Aitken mode
SGINIAC = 2.0   # initial sigma-G for accumulation mode
SGINICO = 2.2   # initial sigma-G for coarse mode
DGINIAT = 0.01E-6  # geometric mean diameter for Aitken mode [ m ]
DGINIAC = 0.07E-6  # geometric mean diameter for accum  mode [ m ]
DGINICO = 1.0E-6   # geometric mean diameter for coarse mode [ m ]
CONMIN  = 1.0E-30  # minimum concentration [ ug/m**3 ]
nspecies = 57      # number of aerosol species treated

#...conversion factors for number and surface area
NUMFAC = {}
NUMFAC['ATKN'] = 1.0 / ( ( DGINIAT ** 3.0 ) * exp( ( 9.0 / 2.0 ) * ( ( log( SGINIAT ) ) ** 2.0 ) ) )
NUMFAC['ACC']  = 1.0 / ( ( DGINIAC ** 3.0 ) * exp( ( 9.0 / 2.0 ) * ( ( log( SGINIAC ) ) ** 2.0 ) ) )
NUMFAC['COR']  = 1.0 / ( ( DGINICO ** 3.0 ) * exp( ( 9.0 / 2.0 ) * ( ( log( SGINICO ) ) ** 2.0 ) ) )
SRFFAC = {}
SRFFAC['ATKN'] = pi / ( DGINIAT * exp( ( 5.0 / 2.0 ) * ( ( log( SGINIAT ) ) ** 2.0 ) ) )
SRFFAC['ACC']  = pi / ( DGINIAC * exp( ( 5.0 / 2.0 ) * ( ( log( SGINIAC ) ) ** 2.0 ) ) )
SRFFAC['COR']  = pi / ( DGINICO * exp( ( 5.0 / 2.0 ) * ( ( log( SGINICO ) ) ** 2.0 ) ) )

class speciesstruct(object):
    def __init__(self, name, ind, mode, density, version, found):
        """
        name - text string identifying aerosol
        ind - variable number (still needed?)
        mode - 0: Aitken, 1: Accumulation, 2: Coarse
        density - aerosol density (kg/m^3)
        version - 0: both, 5: AERO5, 6: AERO6
        """
        self.name = name
        self.mode = mode
        self.density = density
        self.version = version
        self.found = found

def getbcspecies():         
    bcspecies = [speciesstruct (    'ACLI',   0, 'ATKN', RHOSEAS, 0, False),
                 speciesstruct (    'AECI',   0, 'ATKN', RHOANTH, 0, False),
                 speciesstruct (    'ANAI',   0, 'ATKN', RHOSEAS, 0, False),
                 speciesstruct (   'ANH4I',   0, 'ATKN',  RHONH4, 0, False),
                 speciesstruct (   'ANO3I',   0, 'ATKN',  RHONO3, 0, False),
                 speciesstruct (   'ASO4I',   0, 'ATKN',  RHOSO4, 0, False),
                 speciesstruct (    'A25I',   0, 'ATKN', RHOANTH, 5, False),
                 speciesstruct ( 'AORGPAI',   0, 'ATKN',  RHOORG, 5, False),
                 speciesstruct (  'AOTHRI',   0, 'ATKN', RHOANTH, 6, False),
                 speciesstruct ( 'APNCOMI',   0, 'ATKN',  RHOORG, 6, False),
                 speciesstruct (   'APOCI',   0, 'ATKN',  RHOORG, 6, False),
                 speciesstruct (   'AALKJ',   0,  'ACC',  RHOORG, 0, False),
                 speciesstruct (  'ABNZ1J',   0,  'ACC',  RHOORG, 0, False),
                 speciesstruct (  'ABNZ2J',   0,  'ACC',  RHOORG, 0, False),
                 speciesstruct (  'ABNZ3J',   0,  'ACC',  RHOORG, 0, False),
                 speciesstruct (    'ACLJ',   0,  'ACC', RHOSEAS, 0, False),
                 speciesstruct (    'AECJ',   0,  'ACC', RHOANTH, 0, False),
                 speciesstruct (  'AISO1J',   0,  'ACC',  RHOORG, 0, False),
                 speciesstruct (  'AISO2J',   0,  'ACC',  RHOORG, 0, False),
                 speciesstruct (  'AISO3J',   0,  'ACC',  RHOORG, 0, False),
                 speciesstruct (    'ANAJ',   0,  'ACC', RHOSEAS, 0, False),
                 speciesstruct (   'ANH4J',   0,  'ACC',  RHONH4, 0, False),
                 speciesstruct (   'ANO3J',   0,  'ACC',  RHONO3, 0, False),
                 speciesstruct (  'AOLGAJ',   0,  'ACC',  RHOORG, 0, False),
                 speciesstruct (  'AOLGBJ',   0,  'ACC',  RHOORG, 0, False),
                 speciesstruct (  'AORGCJ',   0,  'ACC',  RHOORG, 0, False),
                 speciesstruct (   'ASO4J',   0,  'ACC',  RHOSO4, 0, False),
                 speciesstruct (   'ASQTJ',   0,  'ACC',  RHOORG, 0, False),
                 speciesstruct (  'ATOL1J',   0,  'ACC',  RHOORG, 0, False),
                 speciesstruct (  'ATOL2J',   0,  'ACC',  RHOORG, 0, False),
                 speciesstruct (  'ATOL3J',   0,  'ACC',  RHOORG, 0, False),
                 speciesstruct (  'ATRP1J',   0,  'ACC',  RHOORG, 0, False),
                 speciesstruct (  'ATRP2J',   0,  'ACC',  RHOORG, 0, False),
                 speciesstruct (  'AXYL1J',   0,  'ACC',  RHOORG, 0, False),
                 speciesstruct (  'AXYL2J',   0,  'ACC',  RHOORG, 0, False),
                 speciesstruct (  'AXYL3J',   0,  'ACC',  RHOORG, 0, False),
                 speciesstruct (    'A25J',   0,  'ACC', RHOANTH, 5, False),
                 speciesstruct ( 'AORGPAJ',   0,  'ACC',  RHOORG, 5, False),
                 speciesstruct (    'AALJ',   0,  'ACC', RHOSOIL, 6, False),
                 speciesstruct (    'ACAJ',   0,  'ACC', RHOSOIL, 6, False),
                 speciesstruct (    'AFEJ',   0,  'ACC', RHOSOIL, 6, False),
                 speciesstruct (     'AKJ',   0,  'ACC', RHOSOIL, 6, False),
                 speciesstruct (    'AMGJ',   0,  'ACC', RHOSEAS, 6, False),
                 speciesstruct (    'AMNJ',   0,  'ACC', RHOSOIL, 6, False),
                 speciesstruct (  'AOTHRJ',   0,  'ACC', RHOSOIL, 6, False),
                 speciesstruct ( 'APNCOMJ',   0,  'ACC',  RHOORG, 6, False),
                 speciesstruct (   'APOCJ',   0,  'ACC',  RHOORG, 6, False),
                 speciesstruct (    'ASIJ',   0,  'ACC', RHOSOIL, 6, False),
                 speciesstruct (    'ATIJ',   0,  'ACC', RHOSOIL, 6, False),
                 speciesstruct (    'ACLK',   0,  'COR', RHOSEAS, 0, False),
                 speciesstruct (   'ACORS',   0,  'COR', RHOANTH, 0, False),
                 speciesstruct (   'ANH4K',   0,  'COR',  RHONH4, 0, False),
                 speciesstruct (   'ANO3K',   0,  'COR',  RHONO3, 0, False),
                 speciesstruct (   'ASO4K',   0,  'COR',  RHOSO4, 0, False),
                 speciesstruct (   'ASOIL',   0,  'COR', RHOSOIL, 0, False),
                 speciesstruct (    'ANAK',   0,  'COR', RHOSEAS, 5, False),
                 speciesstruct ( 'ASEACAT',   0,  'COR', RHOSEAS, 6, False)]
    return bcspecies

def repair_ae(f, myioo):
    warn = myioo.warn
    status = myioo.status
    error = myioo.error
    bcspcs = getbcspecies()
    for spc in bcspcs:
        try:
            ntimes, nlayers, nperim = f.variables[spc.name].shape
            break
        except:
            pass
    else:
        warn("There are no aerosol species")
        return
    not_found = defaultdict(lambda: [])
    for spc in bcspcs:
        try:
            spc.found = spc.name in f.variables.keys()
        except:
            spc.found = False
        if not spc.found:
            not_found[spc.version].append(spc.name)
    version_check = []
    for k, v in not_found.iteritems():
        if k != 0:
            version_check.append(k)
        warn('Some variables from %d were not found: %s' % (k, v))
    
    if len(version_check) > 1:
        warn('Some variables from aerosol versions %s were not found' % ' and '.join([str(v) for v in version_check]))
    moment3 = dict([(k, zeros((ntimes, nlayers, nperim), dtype = 'f')) for k in 'ATKN ACC COR'.split()])
    
    for spc in bcspcs:
        try:
            bcval = f.variables[spc.name]
        except KeyError:
            continue
        v = 1.0e-9*6.0/( pi*spc.density ) * bcval[:]
        moment3[spc.mode] += v
    
    for modek, modv in moment3.iteritems():
        numkey = 'NUM' + modek
        srfkey = 'SRF' + modek
        if numkey in f.variables.keys():
            numvar = f.variables[numkey]
        else:
            numvar = f.createVariable(numkey, 'f', ('TSTEP', 'LAY', 'PERIM'))
            numvar.units = '#/m**3'.ljust(16);
            numvar.long_name = numkey.ljust(16);
            numvar.var_desc = numkey.ljust(16);

        if srfkey in f.variables.keys():
            srfvar = f.variables[srfkey]
        else:
            srfvar = f.createVariable(srfkey, 'f', ('TSTEP', 'LAY', 'PERIM'))
            srfvar.units = 'm**2/m**3'.ljust(16);
            srfvar.long_name = srfkey.ljust(16);
            srfvar.var_desc = srfkey.ljust(16);
        numvar[:] = NUMFAC[modek] * modv
        srfvar[:] = SRFFAC[modek] * modv
    f.sync()

if __name__ == '__main__':
    f = Dataset(sys.argv[1], mode = 'r+')
    repair_ae(f)
    