__all__ = 'Dataset NetCDFFile ioapi bpch cspec profile'.split()

from netCDF4 import Dataset
NetCDFFile = Dataset
from PseudoNetCDF.conventions.ioapi import add_cf_from_ioapi
from PseudoNetCDF import getvarpnc, PseudoNetCDFFile
from PseudoNetCDF.cmaqfiles.profile import profile

class ioapi(PseudoNetCDFFile):
    def __init__(self, path):
        outf = getvarpnc(NetCDFFile(path), None)
        add_cf_from_ioapi(outf)
        self.variables = outf.variables
        self.dimensions = outf.dimensions
        for k in outf.ncattrs():
            setattr(self, k, getattr(outf, k))

from PseudoNetCDF.geoschemfiles import bpch, cspec
