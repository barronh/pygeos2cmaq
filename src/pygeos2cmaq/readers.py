__all__ = 'Dataset NetCDFFile METBDY3D METCRO3D ND49NC bpch cspec bcon_profile icon_profile ijavgnc getsigmafromvglvls'.split()

from fast_interp import get_interp_w
from netCDF4 import Dataset
import numpy as np
from warnings import warn
NetCDFFile = Dataset
from PseudoNetCDF.conventions.ioapi import add_cf_from_ioapi
from PseudoNetCDF import getvarpnc, PseudoNetCDFFile
from PseudoNetCDF.cmaqfiles.profile import bcon_profile as pure_bcon_profile, icon_profile as pure_icon_profile
from PseudoNetCDF.geoschemfiles import bpch as pure_bpch, cspec
from PseudoNetCDF.sci_var import extract, interpvars
def getsigmafromvglvls(ifile, shape, layeraxis):
    sigmamid = ifile.VGLVLS.repeat(2, 0)[1:-1].reshape(-1, 2).mean(1)
    return sigmamid[(None,) * layeraxis + (slice(None),) + (None,) * (len(shape) - layeraxis - 1)] * np.ones(shape)
    
class bcon_profile(pure_bcon_profile):
    def __init__(self, path):
        if isinstance(path, str):
            pure_bcon_profile.__init__(self, path)
        else:
            self.variables = path.variables
            self.dimensions = path.dimensions
            for k in path.ncattrs():
                setattr(self, k, getattr(path, k))
        
    def interptosigma(self, sigma, sources):
        if sigma.ndim == 1:
            ntimes = len(self.dimensions['TSTEP'])
            nlays = len(self.dimensions['sigma-mid'])
            nperim = len(self.dimensions['south_east_north_west'])
            shape = [ntimes, nlays, nperim]
            sigma = sigma[None, :, None] * np.ones(shape, dtype = sigma.dtype)
        sigma_in = self.getsigma()
        if np.all(sigma_in == sigma):
            return self
        else:
            weights = get_interp_w(sigma_in, sigma, axis = 1)
            outf = getvarpnc(self, sources['PROFILE'])
            return interpvars(outf, weights = weights, dimension = 'sigma-mid')
    def getsigma(self, vgtop, **kwds):
        return self.variables['sigma-mid']

class icon_profile(pure_icon_profile):
    def __init__(self, path):
        if isinstance(path, str):
            pure_icon_profile.__init__(self, path)
        else:
            self.variables = path.variables
            self.dimensions = path.dimensions
            for k in path.ncattrs():
                setattr(self, k, getattr(path, k))
    def interptosigma(self, sigma, sources):
        if sigma.ndim == 1:
            ntimes = len(self.dimensions['TSTEP'])
            nlays = len(self.dimensions['sigma-mid'])
            nperim = len(self.dimensions['south_east_north_west'])
            shape = [ntimes, nlays, nperim]
            sigma = sigma[None, :, None] * np.ones(shape, dtype = sigma.dtype)
        sigma_in = self.getsigma()
        if np.all(sigma_in == sigma):
            return self
        else:
            weights = get_interp_w(sigma_in, sigma, axis = 1)
            outf = getvarpnc(self, sources['PROFILE'])
            return interpvars(outf, weights = weights, dimension = 'sigma-mid')
    def getsigma(self, vgtop, **kwds):
        return self.variables['sigma-mid']

class METCRO3D(PseudoNetCDFFile):
    def __init__(self, path):
        if isinstance(path, PseudoNetCDFFile):
            outf = path
        else:
            outf = getvarpnc(NetCDFFile(path), ['TFLAG', 'TA', 'PRES'])
        add_cf_from_ioapi(outf)
        self.groups = dict(METCRO3D = outf)
        self.variables = outf.variables
        self.dimensions = outf.dimensions
        outf.PRES_VAR_NAME = 'PRES'
        if 'PERIM' not in self.dimensions:
            self.createDimension('PERIM', sum(map(len, [self.dimensions['ROW'], self.dimensions['COL']])) * 2 + 4)
            
        for k in outf.ncattrs():
            setattr(self, k, getattr(outf, k))

    def interptosigma(self, sigma, sources):
        if sigma.ndim == 1:
            ntimes = len(self.dimensions['TSTEP'])
            nlays = len(self.dimensions['LAY'])
            nrows = len(self.dimensions['ROW'])
            ncols = len(self.dimensions['COL'])
            shape = [ntimes, nlays, nrows, ncols]
            sigma = sigma[None, :, None, None] * np.ones(shape, dtype = sigma.dtype)
        sigma_in = self.getsigma()
        if np.all(sigma_in == sigma):
            return self
        else:
            weights = get_interp_w(sigma_in, sigma, axis = 1)
            outf = getvarpnc(self, sources['METCRO3D'])
            return interpvars(outf, weights = weights, dimension = 'LAY')
            
    def getsigma(self, vgtop, **kwds):
        assert(vgtop == self.VGTOP)
        ntimes = len(self.dimensions['TSTEP'])
        nlays = len(self.dimensions['LAY'])
        nrows = len(self.dimensions['ROW'])
        ncols = len(self.dimensions['COL'])
        shape = [ntimes, nlays, nrows, ncols]
        return getsigmafromvglvls(self, shape = shape, layeraxis = 1)

class bpch(pure_bpch):
    def __init__(self, *args, **kwds):
        if isinstance(args[0], str):
            pure_bpch.__init__(self, *args, **kwds)
        else:
            self.variables = args[0].variables
            self.dimensions = args[0].dimensions
            for k in args[0].ncattrs():
                setattr(self, k, getattr(args[0], k))
        
        
    def tooutcoords(self, coordstr, sigma, vgtop, sources):
        outf = bpch(self.extractcoords(coordstr, sources, unique = False))
        outf = bpch.interptosigma(outf, sigma, vgtop, sources)
        return bpch(outf)
        
    def extractcoords(self, coordstr, sources, unique = False):
        return extract(self, [coordstr], unique = unique)
        
    def interptosigma(self, sigma, vgtop, sources):
        if sigma.ndim == 1:
            ntimes = len(self.dimensions['time'])
            nlays = len(self.dimensions['layer'])
            if 'points' in self.dimensions:
                npoints = len(self.dimensions['points'])
                shape = [ntimes, nlays, nperim]
                sigma = sigma[None, :, None] * np.ones(shape, dtype = sigma.dtype)
            else:
                nrows = len(self.dimensions['latitude'])
                ncols = len(self.dimensions['longitude'])
                shape = [ntimes, nlays, nrows, ncols]
                sigma = sigma[None, :, None, None] * np.ones(shape, dtype = sigma.dtype)
        
        sigma_in = self.getsigma(vgtop)
        if np.all(sigma_in == sigma):
            return self
        else:
            weights = get_interp_w(sigma_in, sigma, axis = 1)
            outf = getvarpnc(self, sources['IJ-AVG-$'])
            return bpch(interpvars(outf, weights = weights, dimension = 'LAY'))
            
        

    def _getshape(self):
        ntimes = len(self.dimensions['time'])
        nlays = len(self.dimensions['layer'])
        nrows = len(self.dimensions['latitude'])
        ncols = len(self.dimensions['longitude'])
        shape = [ntimes, nlays, nrows, ncols]
        return shape

    def getsigma(self, vgtop):
        try:
            sigma = _get_geos_sigma_from_psurf(self, vgtop)
        except:
            warn('Using ETA-Pressures for sigma')
            ipress = self.variables['etai_pressure']
            apress = self.variables['etam_pressure']
            sigma = (apress - vgtop) / (ipress[0] - vgtop)
            shape = self._getshape()
            sigma = sigma[None, :, None, None] * np.ones(shape)
        return sigma
            

class ND49NC(PseudoNetCDFFile):
    def __init__(self, path):
        if isinstance(path, PseudoNetCDFFile):
            outf = path
        else:
            outf = getvarpnc(NetCDFFile(path), None)
        self.groups = dict(ND49 = outf)
        self.variables = outf.variables
        self.dimensions = outf.dimensions
        for k in outf.ncattrs():
            setattr(self, k, getattr(outf, k))
    
    def tooutcoords(self, out, outvert, vgtop, config):
        outf = self.extractcoords(coordstr, sources, unique = False)
        outf = interpextracttosigma(outf, sigma, vgtop, sources)
        return ND49NC(outf)

    def interptosigma(self, sigma, vgtop, sources):
        if sigma.ndim == 1:
            ntimes = len(self.dimensions['time'])
            nlays = len(self.dimensions['layer'])
            if 'points' in self.dimensions:
                npoints = len(self.dimensions['points'])
                shape = [ntimes, nlays, nperim]
                sigma = sigma[None, :, None] * np.ones(shape, dtype = sigma.dtype)
            else:
                nrows = len(self.dimensions['latitude'])
                ncols = len(self.dimensions['longitude'])
                shape = [ntimes, nlays, nrows, ncols]
                sigma = sigma[None, :, None, None] * np.ones(shape, dtype = sigma.dtype)
        sigma_in = self.getsigma(vgtop)
        if np.all(sigma_in == sigma):
            return self
        else:
            weights = get_interp_w(sigma_in, sigma, axis = 1)
            outf = getvarpnc(self, sources['ND49'])
            return ND49NC(interpvars(outf, weights = weights, dimension = 'LAY'))


    def getsigma(self, vgtop, **kwds):
        return _get_geos_sigma_from_psurf(self, vgtop)

def _get_geos_sigma_from_psurf(nf, vgtop):
    if 'PEDGE-$_PSURF' in nf.variables.keys():
        press = nf.variables['PEDGE-$_PSURF'][:] # hPa
    elif 'PSURF' in nf.variables.keys():
        press = nf.variables['PSURF'][:] # hPa
    press = np.concatenate([press, press[:, [0]] * 0 + 0.01], 1)
    avgpress = np.apply_along_axis(func1d = lambda x_: np.convolve(np.array([0.5, 0.5]), x_, mode = 'valid'), arr = press, axis = 1)
    sigma = (avgpress - vgtop) / (press[:, [0]] - vgtop)
    return sigma

class METBDY3D(PseudoNetCDFFile):
    def __init__(self, path):
        if isinstance(path, PseudoNetCDFFile):
            outf = path
        else:
            outf = getvarpnc(NetCDFFile(path), ['TFLAG', 'TA', 'PRES'])
            add_cf_from_ioapi(outf)
        outf.PRES_VAR_NAME = 'PRES'
        self.groups = dict(METBDY3D = outf)
        self.variables = outf.variables
        self.dimensions = outf.dimensions
        if 'PERIM' not in self.dimensions:
            self.createDimension('PERIM', sum(map(len, [self.dimensions['ROW'], self.dimensions['COL']])) * 2 + 4)
            
        for k in outf.ncattrs():
            setattr(self, k, getattr(outf, k))

    def interptosigma(self, sigma, sources):
        if sigma.ndim == 1:
            ntimes = len(self.dimensions['TSTEP'])
            nlays = len(self.dimensions['LAY'])
            nperim = len(self.dimensions['PERIM'])
            shape = [ntimes, nlays, nperim]
            sigma = sigma[None, :, None] * np.ones(shape, dtype = sigma.dtype)
        sigma_in = self.getsigma()
        if np.all(sigma_in == sigma):
            return self
        else:
            weights = get_interp_w(sigma_in, sigma, axis = 1)
            outf = getvarpnc(self, sources['METBDY3D'])
            return METBDY3D(interpvars(outf, weights = weights, dimension = 'LAY'))
    
    def getsigma(self, vgtop, **kwds):
        assert(vgtop == self.VGTOP)
        ntimes = len(self.dimensions['TSTEP'])
        nlays = len(self.dimensions['LAY'])
        nperim = len(self.dimensions['PERIM'])
        shape = [ntimes, nlays, nperim]
        return getsigmafromvglvls(self, shape = shape, layeraxis = 1)

class ijavgnc(PseudoNetCDFFile):
    def __init__(self, path):
        if isinstance(path, PseudoNetCDFFile):
            outf = path
        else:
            outf = getvarpnc(NetCDFFile(path), None)
        self.groups = {'IJ-AVG-$': outf}
        self.variables = outf.variables
        self.dimensions = outf.dimensions
        for k in outf.ncattrs():
            setattr(self, k, getattr(outf, k))
    def getsigma(self, vgtop, **kwds):
        etam_pressure = self.variables['etam_pressure'][:]
        etai_pressure = self.variables['etai_pressure'][:]
        sigma = (etam_pressure - vgtop) / (etai_pressure[0] - vgtop)
        ntimes = len(self.dimesions['time'])
        nlays = len(self.dimesions['layer'])
        nrows = len(self.dimesions['latitude'])
        ncols = len(self.dimesions['longitude'])
        return sigma[None, :, None, None] * np.ones([ntimes, nlays, nrows, ncols])
