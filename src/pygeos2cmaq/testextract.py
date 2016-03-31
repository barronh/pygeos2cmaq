from collections import OrderedDict
import numpy as np
from PseudoNetCDF import interpvars, getvarpnc, extract
from netCDF4 import Dataset
from PseudoNetCDF.geoschemfiles import bpch
f = bpch('testdata/ts20100101.bpch');
metf = Dataset('testdata/METBDY3D_100101')
oldbc = Dataset('geos2cmaq.20100101.nc')
def pres_from_sigma(sigma, pref, ptop, avg = False):
    pres = sigma * (pref - ptop) + ptop
    if avg:
        pres = pres[:-1] + np.diff(pres) / 2.
    return pres

vert_out = pres_from_sigma(metf.VGLVLS, 101325., metf.VGTOP, avg = True)
vert_in = f.variables['layer'][:]
from fast_interp import get_interp_w
w = get_interp_w(vert_in, vert_out)
print -1

ft = getvarpnc(f.groups['IJ-AVG-$'], None)
xf = extract(ft, [oldbc.lonlatcoords], unique = True); 
print 0
xft = getvarpnc(xf, None)
print 1
xf1 = extract(xft, [oldbc.lonlatcoords])
print 2

def check(f):
    lonlat = f.lonlatcoords.split('/')
    lat, lon = f.variables['geos_latitude_bounds'], f.variables['geos_longitude_bounds']
    mid = lat.shape[0] // 2
    for idx in range(len(lonlat)):
        print idx,
        inlon, inlat = eval(lonlat[idx])
        checklon = inlon >= lon[idx][0] and inlon < lon[idx][1]
        checklat = inlat >= lat[idx][0] and inlat < lat[idx][1]
        assert(checklon and checklat)
    print
check(oldbc)

#check(xf)
#check(xft)
#check(xf1)

nc = oldbc.NCOLS
nr = oldbc.NROWS
ss = 0
es = ss + nc + 1
se = es
ee = se + nr + 1
sn = ee
en = ee + nc + 1
sw = en
ew = sw + nr + 1
latitude = oldbc.variables['geos_latitude'][:]
longitude = oldbc.variables['geos_longitude'][:]

for edges, edgee in [(ss, es), (se, ee), (sn, en), (sw, ew)]:
    print OrderedDict([(k, 0) for k in zip(longitude[edges:edgee], latitude[edges:edgee])]).keys()
