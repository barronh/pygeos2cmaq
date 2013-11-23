import warnings

def formatwarning(message, category, filename, lineno, line = 0):
    strout = "\n***********\n%s:%s: %s:\n\t%s\n***********\n" % (filename, lineno, category.__name__, message)
    return strout

warnings.formatwarning = formatwarning

warn = warnings.warn
from collections import defaultdict
import numpy as np
from numpy import *
from readers import *
from datetime import datetime, timedelta
from fast_interp import get_interp_w
from PseudoNetCDF import PseudoNetCDFFile

def timeit(name, new):
    from time import time
    if new:
        timeit.ts.append(time())
        timeit.names.append(name)
        print name
    else:
        print '***%s min:' % timeit.names.pop(-1), (time() - timeit.ts.pop(-1)) / 60
timeit.ts = []
timeit.names = []
def process(config, verbose = False):
    """
    Take the configuration file and create boundaries
    """
    alldates = get_dates(config)
    outpaths = defaultdict(lambda: [])
    outpathtemp, tsf = config['out_template']
    for date in alldates:
        outpaths[eval(tsf)(date, outpathtemp)].append(date)
    outpaths = [(k, v) for k, v in outpaths.iteritems()]
    outpaths.sort()
    errors = set()
    for outpath, dates in outpaths:
        print outpath
        out = make_out(config, dates)
        tflag = out.variables['TFLAG']
        curdate = 0
        for di, date in enumerate(dates):
            file_objs = get_files(config, date, out.lonlatcoords)
            jday = int(date.strftime('%Y%j'))
            itime = int(date.strftime('%H%M%S'))
            tflag[di, :, :] = np.array([[jday, itime]])
        for src, name, expr, outunit in config['mappings']:
            print '\t'+name
            try:
                grp = get_group(file_objs, src, dates)
                if verbose: timeit('EVALUNIT', True)
                evalandunit(out, curdate, name, expr, grp.variables, verbose = verbose)
                if verbose: timeit('EVALUNIT', False)
            except Exception, e:
                errors.add((src, name, expr, str(e)))
                warn("Unable to map %s to %s in %s: %s" % (name, expr, src, str(e)), stacklevel = 1)
        curdate += len(file_objs[0].dimensions['time'])
        output(out, outpath, config)
    if len(errors) > 0:
        print '******************'
        print '**** Start Errors:'
        print '******************'
    for src, name, expr, err in errors:
        print src, name, expr, err
    if len(errors) > 0:
        print '******************'
        print '**** End Errors:'
        print '******************'
def output(out, outpath, config):
    """
    """
    from netCDF4 import Dataset
    from repair_ae import repair_ae
    from PseudoNetCDF.pncgen import pncgen
    tmp = defaultdict(lambda: 1)
    for k, v in config['unitconversions']['metadefs'].iteritems():
        eval(v, None, tmp)
    for k in tmp.keys():
        tmp[k] = out.variables[k]
    for k, v in config['unitconversions']['metadefs'].iteritems():
        exec('%s = %s' % (k, v), tmp, out.variables)
    print 'Vertical profile Low -> High'
    np.set_printoptions(precision = 2)
    for vark, varo in out.variables.items():
        if hasattr(varo, 'unitnow'):
            if varo.unitnow.strip() != varo.units.strip():
                expr = config['unitconversions']['%s->%s' % (varo.unitnow.strip(), varo.units.strip())].replace('<value>', 'varo')
                exec('varo[:] = %s' % expr, dict(varo = varo), out.variables)
                varo.history += ';' + expr.replace('varo', 'RESULT')
                if (varo[:, :] < 0).any():
                    warn(vark + ' has negative values', stacklevel = 1)
            del varo.unitnow
        if 'TSTEP' in varo.dimensions and 'PERIM' in varo.dimensions:
            print vark, varo[:, :].mean(0).mean(1)
    np.set_printoptions(precision = None)
    
    f = pncgen(out, outpath, inmode = 'r', outmode = 'w', format = 'NETCDF4_CLASSIC', verbose = False)
    f.sync()
    f.close()
    if config['repair_aero']:
        f = Dataset(outpath, 'r+')
        repair_ae(f)
        f.sync()
        f.close()
    #del f, out

def pres_from_sigma(sigma, pref, ptop, avg = False):
    pres = sigma * (pref - ptop) + ptop
    if avg:
        pres = pres[:-1] + np.diff(pres) / 2.
    return pres
    
def evalandunit(out, di, name, expr, variables, verbose = False):
    pref = 101325.
    vert_out = pres_from_sigma(out.VGLVLS, pref, out.VGTOP, avg = True)
    if 'sigma-mid' in variables.keys():
        sigma = np.array(variables['sigma-mid']).repeat(2, 0)[1:-1].reshape(-1, 2).mean(1)
        vert_in = pres_from_sigma(sigma, pref, out.VGTOP)
    elif 'layer' in variables.keys():
        vert_in = variables['layer'] * 100.
    else:
        raise KeyError('No sigma-mid or layer in %s' % str(variables.keys()))
    
    x = defaultdict(lambda: 1)
    eval(expr, None, x)
    key = x.keys()[0]
    metavar = variables[key]
    origunits = metavar.units.strip()
    outvar = out.variables[name]
    outunit = outvar.units.strip()
    val = eval(expr, None, variables)
    if verbose: timeit('VINTERP', True)
    outval = vinterp(val, vert_in, vert_out)
    if verbose: timeit('VINTERP', False)
    unitnow = origunits
    outvar.history = expr
    if unitnow == "None":
        warn('No unit was provided; assuming CMAQ unit; Most likely from PROFILE', stacklevel = 1)
    else:
        if outunit in ('micrograms/m**3',):
            outval *= metavar.kgpermole
            outvar.history += '; RESULT * %s' % metavar.kgpermole
            unitnow = ('kg/mol*%s' % origunits).replace('kg/mol*ppbC', 'micrograms/mol').replace('kg/mol*ppbv', 'micrograms/mol')
        
            # Special case where profile is actually assumed to be in appropriate units
            unitnow = unitnow.replace('kg/mol*ppmV', 'micrograms/m**3')
        elif origunits in ('ppbC',):
            outval /= metavar.carbon
            unitnow = 'ppbv'
            outvar.history += '; RESULT / %s' % metavar.carbon
        outvar.unitnow = unitnow
    outvar[di:di + outvar.shape[0]] += outval.squeeze()
    outvar.origunits = origunits
    for k in metavar.ncattrs():
        if k not in outvar.ncattrs():
            setattr(outvar, k, getattr(metavar, k))

    return outval, origunits
    
outval = np.zeros((0,))
def vinterp(val, vert_in, vert_out):
    """
    Performs vertical interpolation using linear algorithms
    """
    global outval
    if (outval.shape[0], outval.shape[-1]) != (val.shape[0], val.shape[-1]):
        outval = zeros((val.shape[0], len(vert_out), val.shape[-1]), dtype = val.dtype)
    w = get_interp_w(vert_in, vert_out)
    newvals = (w[:, :, None, None] * val[:].swapaxes(0, 1)[None, :]).sum(1).swapaxes(0, 1)
    import pdb; pdb.set_trace()
    outval[:, :, :] = newvals
    if not len(vert_in) == 6:
        for ti in [0, val.shape[0] - 1]:
            for pi in [0, val.shape[-1] - 1]:
                if vert_out.max() > vert_in.max():
                    right = val[ti, 0, pi] + np.diff(val[ti, :2, pi][::-1]) / np.diff(vert_in[:2][::-1]) * (vert_out[0] - vert_in[0])
                else:
                    right = None
            
                if vert_out.min() > vert_in.min():
                    left = val[ti, -1, pi] + np.diff(val[ti, -2:, pi][::-1]) / np.diff(vert_in[-2:][::-1]) * (vert_out[-1] - vert_in[-1])
                else:
                    left = None
            
                testval = np.interp(vert_out[::-1], vert_in[:val.shape[1]][::-1], val[ti, ::-1, pi], left = left, right = right)[::-1]
                try:
                    np.testing.assert_allclose(newvals[ti, :, pi], testval, rtol=1e-05, atol=0, err_msg='', verbose=True)
                except Exception, e:
                    warn(str(e))
            

    return outval

def make_out(config, dates):
    """
    Return a file to fill with data
    """
    out = PseudoNetCDFFile()
    file_objs = get_files(config, dates[0], None)
    metf = [f for f in file_objs if 'PERIM' in f.dimensions][0]
    d = out.createDimension('TSTEP', len(dates))
    d.setunlimited(True)
    out.createDimension('LAY', len(metf.dimensions['LAY']))
    out.createDimension('PERIM', len(metf.dimensions['PERIM']))
    out.createDimension('DATE-TIME', len(metf.dimensions['DATE-TIME']))
    out.createDimension('VAR', len(config['mappings']))
    mlat = metf.variables['latitude']
    lat = out.createVariable('latitude', 'f', ('PERIM',), units = mlat.units, values = mlat[:])

    mlon = metf.variables['longitude']
    lon = out.createVariable('longitude', 'f', ('PERIM',), units = mlon.units, values = mlon[:])

    var = out.createVariable('TFLAG', 'i', ('TSTEP', 'VAR', 'DATE-TIME'))
    for pk in metf.ncattrs():
        setattr(out, pk, getattr(metf, pk))
    setattr(out, 'VAR-LIST', ''.join([name.ljust(16) for src, name, expr, unit in config['mappings']]))
    setattr(out, 'NVARS', len(config['mappings']))
    setattr(out, 'SDATE', int(dates[0].strftime('%Y%j')))
    setattr(out, 'STIME', int(dates[0].strftime('%H%M%S')))
    setattr(out, 'EDATE', int(dates[-1].strftime('%Y%j')))
    setattr(out, 'ETIME', int(dates[-1].strftime('%H%M%S')))
    var.long_name = 'TFLAG'.ljust(16);
    var.var_desc = "Timestep-valid flags:  (1) YYYYDDD or (2) HHMMSS".ljust(80)
    var.units = "<YYYYDDD,HHMMSS>"
    for src, name, expr, outunit in config['mappings']:
        var = out.createVariable(name, 'f', ('TSTEP', 'LAY', 'PERIM'))
        var.long_name = name.ljust(16);
        var.var_desc = name.ljust(80);
        var.units = outunit.ljust(16)
    out.lonlatcoords = '/'.join(['%s,%s' % (o, a) for o, a in zip(mlon, mlat)])
    return out

def get_group(file_objs, src, dates):
    """
    Return the first file with the SRC group
    """
    from PseudoNetCDF.sci_var import extract, slice_dim
    for f in file_objs:
        if hasattr(f, 'groups'):
            if src in f.groups:
                grp = f.groups[src]
                if 'TFLAG' in grp.variables.keys():
                    tflag = grp.variables['TFLAG']
                    thisdate = [date.strptime('%d %06d' % (d, t), '%Y%j %H%M%S') for d, t in tflag[:, 0]]
                elif 'time' in grp.variables.keys():
                    time = grp.variables['time']
                    unit, starttxt = grp.variables['time'].units.split(' since ')
                    starttxt = starttxt.replace('UTC', '').strip()
                    start_date = datetime.strptime(starttxt, '%Y-%m-%d %H:%M:%S')
                    thisdate = [start_date + timedelta(**{unit: t}) for t in time[:]]
                elif 'tau0' in grp.variables.keys():
                    time = grp.variables['tau0']
                    unit, starttxt = time.units.split(' since ')
                    start_date = datetime.strptime(start_date)
                    thisdate = [start_date + timedelta(**{unit: t}) for t in time]
                else:
                    return grp
                thisdate = np.array(thisdate)
                idx, = np.where(np.logical_and(thisdate >= dates[0], thisdate <= dates[-1]))
                grpt = slice_dim(grp, 'time,%d,%d' % (idx[0], idx[-1] + 1))
                return grpt
                
    else:
        raise KeyError('No file provided has the %s group' % src)

def simpledate(date, p):
    return date.strftime(p)
def minus1hour(date, p):
    return (date - timedelta(hours = 1)).strftime(p)

last_coordstr = "-1-"
def get_files(config, date, coordstr):
    """
    Put date into file path and return open files
    """
    from PseudoNetCDF.sci_var import extract, slice_dim, getvarpnc
    from PseudoNetCDF.cmaqfiles.profile import profile
    import gc
    global last_file_paths
    global file_objs
    global last_coordstr
    global last_file_objs
    file_paths = [(r, eval(tsf)(date, p)) for r, p, tsf in config['file_templates']]
    if coordstr != last_coordstr:
        last_file_paths = [''] * len(file_paths)
        last_file_objs = [PseudoNetCDFFile()] * len(file_paths)
        gc.collect()
    if file_paths == last_file_paths:
        return last_file_objs
    else:
        if coordstr is not None:
            print 
            print '-' * 40
            print "Getting files for", date
            print '-' * 40
            print 
        for fi, (r, p) in enumerate(file_paths):
            lp = last_file_paths[fi]
            nf = last_file_objs[fi]
            if p != lp:
                try:
                    nf.close()
                    onf = nf = eval(r)(p)
                    if coordstr is not None:
                        if isinstance(nf, profile):
                            pnf = nf
                            metf = [f for f in last_file_objs if 'PERIM' in f.dimensions][0]
                            nc, nr = metf.NCOLS, metf.NROWS
                            nf.createDimension('PERIM', (nc + nr + 2) * 2)
                            for k, v in nf.variables.items():
                                if k in ('sigma', 'sigma-mid'): continue
                                nv = np.concatenate([v[:, [0]].repeat(nc + 1, 1), v[:, [1]].repeat(nr + 1, 1), v[:, [2]].repeat(nc + 1, 1), v[:, [3]].repeat(nr + 1, 1)], axis = 1)
                                vard = dict([(pk, getattr(v, pk)) for pk in v.ncattrs()])
                                nf.createVariable(k, v.dtype.char, ('time', 'sigma-mid', 'PERIM'), values = nv[None,:], kgpermole = 1., **vard)
                            nf.groups = dict(PROFILE = nf)
                        else:
                            if hasattr(nf, 'groups'):
                                for grpk, grpv in nf.groups.items():
                                    nf.groups[grpk] = extract(grpv, [coordstr])
                            else:
                                nf = extract(nf, [coordstr])
                    else:
                        with warnings.catch_warnings():
                            warnings.simplefilter("ignore")
                            nf = getvarpnc(onf, 'time TFLAG latitude longitude'.split())
                    if 'PERIM' in onf.dimensions.keys() and not isinstance(onf, profile):
                        nf.createDimension('PERIM', len(onf.dimensions['PERIM']))
                        nf.createDimension('LAY', len(onf.dimensions['LAY']))
                        nf.NCOLS = onf.NCOLS
                        nf.NROWS = onf.NROWS
                except Exception, e:
                    raise Exception("Could not open %s with %s: %s" % (p, r, str(e)))
            last_file_objs[fi] = nf
        last_file_paths = file_paths
        last_coordstr = coordstr
        return last_file_objs

def get_dates(config):
    """
    Get all dates
    """
    end_date = config['end_date']
    time_incr = config['time_incr']
    alldates = []
    date = config['start_date']
    
    while date <= end_date:
        alldates.append(date)
        date = date + time_incr
    return alldates
