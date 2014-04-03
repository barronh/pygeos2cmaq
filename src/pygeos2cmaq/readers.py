__all__ = 'Dataset NetCDFFile METBDY3D METCRO3D bpch cspec bcon_profile icon_profile ijavgnc'.split()

from netCDF4 import Dataset
NetCDFFile = Dataset
from PseudoNetCDF.conventions.ioapi import add_cf_from_ioapi
from PseudoNetCDF import getvarpnc, PseudoNetCDFFile
from PseudoNetCDF.cmaqfiles.profile import bcon_profile, icon_profile

class METCRO3D(PseudoNetCDFFile):
    def __init__(self, path):
        outf = getvarpnc(NetCDFFile(path), ['TFLAG', 'TA', 'PRES'])
        add_cf_from_ioapi(outf)
        self.groups = dict(METCRO3D = outf)
        self.variables = outf.variables
        self.dimensions = outf.dimensions
        if 'PERIM' not in self.dimensions:
            self.createDimension('PERIM', sum(map(len, [self.dimensions['ROW'], self.dimensions['COL']])) * 2 + 4)
            
        for k in outf.ncattrs():
            setattr(self, k, getattr(outf, k))


class METBDY3D(PseudoNetCDFFile):
    def __init__(self, path):
        outf = getvarpnc(NetCDFFile(path), ['TFLAG', 'TA', 'PRES'])
        add_cf_from_ioapi(outf)
        self.groups = dict(METBDY3D = outf)
        self.variables = outf.variables
        self.dimensions = outf.dimensions
        if 'PERIM' not in self.dimensions:
            self.createDimension('PERIM', sum(map(len, [self.dimensions['ROW'], self.dimensions['COL']])) * 2 + 4)
            
        for k in outf.ncattrs():
            setattr(self, k, getattr(outf, k))

class ijavgnc(PseudoNetCDFFile):
    def __init__(self, path):
        outf = getvarpnc(NetCDFFile(path), None)
        self.groups = {'IJ-AVG-$': outf}
        self.variables = outf.variables
        self.dimensions = outf.dimensions
        for k in outf.ncattrs():
            setattr(self, k, getattr(outf, k))

from PseudoNetCDF.geoschemfiles import bpch, cspec
