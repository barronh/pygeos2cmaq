__all__ = 'NetCDFFile ioapi bpch cspec'.split()

from netCDF4 import Dataset as NetCDFFile
from PseudoNetCDF.conventions.ioapi import add_cf_from_ioapi
from PseudoNetCDF import getvarpnc

def ioapi(path):
    outf = getvarpnc(NetCDFFile(path), None)
    add_cf_from_ioapi(outf)
    return outf

from PseudoNetCDF.geoschemfiles import bpch, cspec
