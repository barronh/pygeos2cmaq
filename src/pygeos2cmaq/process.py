import warnings
from collections import defaultdict, OrderedDict
from string import uppercase
from datetime import datetime, timedelta
import os

import numpy as np
from numpy import *
from readers import *
from timeformatters import *
from fast_interp import get_interp_w
from PseudoNetCDF import PseudoNetCDFFile, interpvars, Pseudo2NetCDF, PseudoNetCDFVariable
from PseudoNetCDF.core import convolve_dim

from myio import myio
myioo = myio()
warn = myioo.warn
status = myioo.status
error = myioo.error
def timeit(name, new):
    from time import time
    if new:
        timeit.ts.append(time())
        timeit.names.append(name)
        status(name)
    else:
        status('***%s min: %s' % (timeit.names.pop(-1), (time() - timeit.ts.pop(-1)) / 60))
timeit.ts = []
timeit.names = []

def process(config, verbose = 0):
    """
    Take the configuration file and create boundaries
    
    This routine:
        1. gets input files (pre vertically interpolated and variable screened; see get_files)
        2. performs algebraic mapping of species
        3. does unit conversion
    
    """
    import numpy as np
    if verbose > 1: timeit('PROCESSING', True)
    alldates = get_dates(config)
    outpaths = defaultdict(lambda: [])
    outpathtemp, tsf = config['out_template']
    for date in alldates:
        outpaths[eval(tsf)(date, outpathtemp)].append(date)
    sources = defaultdict(lambda: set('time TFLAG tau0 tau1 latitude longitude latitude_bounds longitude_bounds'.split()))
    
    for src, name, expr, unit in config['mappings']:
        co = compile(expr, 'src', 'eval')
        sources[src].update(set([k_ for k_ in co.co_names if k_ not in ('np',)]))
    outpaths = [(k, v) for k, v in outpaths.iteritems()]
    outpaths.sort()
    errors = set()
    for outpath, dates in outpaths:
        myioo.clear()
        if config['no_clobber']:
            if os.path.exists(outpath):
                status('Already exists, skipping %s ...' % outpath)
                continue
        timeit(outpath, True)
        out = make_out(config, dates)
        tflag = out.variables['TFLAG']
        curdate = 0
        if verbose > 0: timeit('GET_FILES', True)
        get_files = file_getter(config = config, out = out, sources = sources, verbose = verbose).get_files
        for di, date in enumerate(dates):
            file_objs = get_files(date)
            jday = int(date.strftime('%Y%j'))
            itime = int(date.strftime('%H%M%S'))
            tflag[di, :, :] = np.array([[jday, itime]])
        if verbose > 0: timeit('GET_FILES', False)
        for src, name, expr, outunit in config['mappings']:
            if verbose > 1: timeit('MAP %s' % name, True)
            try:
                grp = get_group(file_objs, src, dates)
                evalandunit(out, curdate, name, expr, grp.variables, verbose = verbose)
            except Exception, e:
                errors.add((src, name, expr, str(e)))
                error("Unable to map %s to %s in %s: %s" % (name, expr, src, str(e)), stacklevel = 1)
            if verbose > 1: timeit('MAP %s' % name, False)

        if config['time_incr'].total_seconds() > 0:
            try:
                curdate += len(file_objs[0].dimensions['time'])
            except:
                curdate += len(file_objs[0].dimensions['TSTEP'])
        if verbose > 2: timeit('OUTPUT', True)
        output(out, outpath, config, verbose = verbose)
        if verbose > 2: timeit('OUTPUT', False)
        timeit(outpath, False)
    if len(errors) > 0:
        status('******************')
        status('**** Start Errors:')
        status('******************')
    for src, name, expr, err in errors:
        status(' '.join([src, name, expr, err]))
    if len(errors) > 0:
        status('******************')
        status('**** End Errors:')
        status('******************')
    if verbose > 1: timeit('PROCESSING', False)

def output(out, outpath, config, verbose = 0):
    """
    Make final unit conversions, add log information
    and save output to disk
    
      out     - PseudoNetCDFFile with output data
      outpath - Path for output
      config  - configuration dictionary
    """
    from repair_ae import repair_ae
    from PseudoNetCDF.pncgen import pncgen
    tmp = defaultdict(lambda: 1)
    tmp['np'] = np
    for k, v in config['unitconversions']['metadefs'].iteritems():
        eval(v, None, tmp)
    for k in tmp.keys():
        if k != 'np':
            tmp[k] = out.variables[k]
    for k, v in config['unitconversions']['metadefs'].iteritems():
        exec('%s = %s' % (k, v), tmp, out.variables)
    if verbose > 0: status('Vertical profile Low -> High')
    np.set_printoptions(precision = 2)
    for vark, varo in out.variables.items():
        if hasattr(varo, 'unitnow'):
            if varo.unitnow.strip() != varo.units.strip():
                expr = config['unitconversions']['%s->%s' % (varo.unitnow.strip(), varo.units.strip())].replace('<value>', 'varo')
                exec('varo[:] = %s' % expr, dict(varo = varo), out.variables)
                varo.history += ';' + expr.replace('varo', 'RESULT')
            del varo.unitnow

        if vark not in ('time', 'TFLAG') and 'latitude' not in vark and 'longitude' not in vark and (varo[...] <= 0).any():
            if config['zero_negs']:
                have_zeros = varo[...] <= 0
                varo[np.where(have_zeros)] = 1.e-30
                warn(vark + ' %d (of %d; %.2f%%) negative or zero values were set to 1e-30' % (have_zeros.sum(), have_zeros.size, have_zeros.mean() * 100), stacklevel = 1)
            else:
                warn(vark + ' has negative or zero values', stacklevel = 1)
        if verbose > 0 and 'TSTEP' in varo.dimensions and 'PERIM' in varo.dimensions:
            status('%s: %s' % (vark, str(varo[:, :].mean(0).mean(1))))
    np.set_printoptions(precision = None)
    out.geos2cmaq_warnings = myioo.getwarnings()
    out.geos2cmaq_errors = myioo.geterrors()
    out.geos2cmaq_status = myioo.getstatus()
    varnames = [k.ljust(16) for k in out.variables.keys() if k[:1] in uppercase and k != 'TFLAG']
    out.createDimension('VAR', len(varnames))
    tflag = out.variables['TFLAG']
    var = out.createVariable('TFLAG', 'i', ('TSTEP', 'VAR', 'DATE-TIME'))
    var.long_name = 'TFLAG'.ljust(16);
    var.var_desc = "Timestep-valid flags:  (1) YYYYDDD or (2) HHMMSS".ljust(80)
    var.units = "<YYYYDDD,HHMMSS>"
    var[:, :, :] = tflag[:, [0], :].repeat(var.shape[1], 1)
    f = pncgen(out, outpath, inmode = 'r', outmode = 'w', format = 'NETCDF3_CLASSIC', verbose = False)
    setattr(f, 'VAR-LIST', ''.join(varnames))
    setattr(f, 'NVARS', len(varnames))
    f.sync()
    f.close()
    if config['repair_aero']:
        f = Dataset(outpath, 'r+')
        repair_ae(f, myioo)
        f.sync()
        f.close()

def evalandunit(out, di, name, expr, variables, verbose = 0):
    """
    Evaluates an expression in the context of a dictionary
    and extracts units from variables used in the evaluation
    """
    
    # Create a dictionary that will create
    # keys in response to an expression
    # add the numpy library as np
    tmpdict = defaultdict(lambda: 1, np = np)
    
    # Evaluate the expression to generate keys
    # used in the expression
    co = compile(expr, expr, 'eval')

    # Get first key to pull meta data
    for k in co.co_names:
        if k in variables:
            # Get a variable to pull meta data
            metavar = variables[k]
            break
    
    # Get the original units from the meta variable
    origunits = metavar.units.strip()
    
    # Get the output variables and store the intended
    # units 
    outvar = out.variables[name]
    outunit = outvar.units.strip()
    # Evaluate the expression to create output data
    outval = eval(expr, None, variables).reshape(*(-1,) + outvar[:].shape[1:])
    
    # Store present unit and add a history to the output variable
    unitnow = origunits
    outvar.history = expr
    
    if unitnow == "None":
        # Units are unknown, which indicates the PROFILE input
        warn('No unit was provided; assuming CMAQ unit; Most likely from PROFILE', stacklevel = 1)
    else:
        if outunit in ('micrograms/m**3',):
            # For mass concentration units, the kg/mol varies
            # per species, so this multiplier must be applied here
            outval *= metavar.kgpermole
            outvar.history += '; RESULT * %s' % metavar.kgpermole
            
            # Update present unit to reflect multiplication
            unitnow = ('kg/mol*%s' % origunits)
        
        elif origunits in ('ppbC',):
            # For units of ppbC, divide by carbon and update unit
            outval /= metavar.carbon
            unitnow = 'ppbv'
            outvar.history += '; RESULT / %s' % metavar.carbon
        # Store present unit in variable
        outvar.unitnow = unitnow
        
    
    # Store output variable in outvar starting at last
    # output
    ts = outval.shape[0]
    outvar[di:di + ts] += outval[:outvar.shape[0] - di].squeeze()

    # Update variables original unit
    # to help in final unit conversion
    outvar.origunits = origunits
    
    # Add properties from the metavariable
    # Note that this will overwrite properties 
    # from previous mapping
    for k in metavar.ncattrs():
        if k not in outvar.ncattrs():
            setattr(outvar, k, getattr(metavar, k))

    return outval, origunits
    
def make_out(config, dates):
    """
    Make an output file with appropriate dimensions and variables

      config - configuration dictionary
      dates  - iterable of dates
    """
    
    from PseudoNetCDF.sci_var import extract
    out = PseudoNetCDFFile()
    
    # Ordered entries are necessary for 
    # consistency with IOAPI
    out.dimensions = OrderedDict()
    out.variables = OrderedDict()

    # Get files for the first date
    get_files = file_getter(config = config, out = None, sources = None).get_files
    file_objs = get_files(dates[0])
    metf = [f for f in file_objs if 'PERIM' in f.dimensions][0]
    outnames = OrderedDict()
    for src, name, expr, unit in config['mappings']:
        if not [src, name, expr, unit] == ['SOURCE', 'MECHSPC', 'GEOS_EXPRESSION', 'UNIT']:
            outnames[name] = 0

    d = out.createDimension('TSTEP', len(dates))
    d.setunlimited(True)
    out.createDimension('DATE-TIME', len(metf.dimensions['DATE-TIME']))
    out.createDimension('LAY', len(metf.dimensions['LAY']))
    out.createDimension('VAR', len(outnames))
    mlatb = metf.variables['latitude_bounds']
    mlonb = metf.variables['longitude_bounds']
    if metf.FTYPE == 2:
        out.createDimension('PERIM', len(metf.dimensions['PERIM']))
        out.createDimension('nv', len(metf.dimensions['nv']))
        vardims = ('TSTEP', 'LAY', 'PERIM')
        coordbounddims = ('PERIM', 'nv')
        coorddims = ('PERIM',)
        out.FTYPE = 2 # Boundary http://www.baronams.com/products/ioapi/TUTORIAL.html
    elif metf.FTYPE == 1:
        out.createDimension('ROW', metf.NROWS)
        out.createDimension('COL', metf.NCOLS)
        out.createDimension('nv', len(metf.dimensions['nv']))
        vardims = ('TSTEP', 'LAY', 'ROW', 'COL')
        coordbounddims = ('ROW', 'COL', 'nv')
        coorddims = ('ROW', 'COL')
        out.FTYPE = 1 # Gridded http://www.baronams.com/products/ioapi/TUTORIAL.html
    else:
        raise ValueError('extract_type must be icon or bcon; got %s' % extract_type)
        

    out.createVariable('latitude_bounds', 'f', coordbounddims, units = mlatb.units, values = mlatb[:])

    out.createVariable('longitude_bounds', 'f', coordbounddims, units = mlonb.units, values = mlonb[:])

    mlat = metf.variables['latitude']
    out.createVariable('latitude', 'f', coorddims, units = mlat.units, values = mlat[:])

    mlon = metf.variables['longitude']
    out.createVariable('longitude', 'f', coorddims, units = mlon.units, values = mlon[:])
    coordstr = '/'.join(['%s,%s' % (o, a) for o, a in zip(mlon[:].ravel(), mlat[:].ravel())])
    geosfs = [f for f in file_objs if 'tau0' in f.variables.keys()]
    if len(geosfs) > 0:
        geosf = geosfs[0]
        geosf = extract(geosf, [coordstr])
        glatb = geosf.variables['latitude_bounds']
        out.createVariable('geos_latitude_bounds', 'f', coordbounddims, units = glatb.units, values = glatb[:, [0, 0, 1, 1]].reshape(mlatb[:].shape))

        glonb = geosf.variables['longitude_bounds']
        out.createVariable('geos_longitude_bounds', 'f', coordbounddims, units = glonb.units, values = glonb[:, [0, 1, 1, 0]].reshape(mlonb[:].shape))

        glat = geosf.variables['latitude']
        out.createVariable('geos_latitude', 'f', coorddims, units = glat.units, values = glat[:].reshape(mlat[:].shape))

        glon = geosf.variables['longitude']
        out.createVariable('geos_longitude', 'f', coorddims, units = glon.units, values = glon[:].reshape(mlon[:].shape))

    var = out.createVariable('TFLAG', 'i', ('TSTEP', 'VAR', 'DATE-TIME'))
    var.long_name = 'TFLAG'.ljust(16);
    var.var_desc = "Timestep-valid flags:  (1) YYYYDDD or (2) HHMMSS".ljust(80)
    var.units = "<YYYYDDD,HHMMSS>"
    for pk in metf.ncattrs():
        setattr(out, pk, getattr(metf, pk))
    setattr(out, 'VAR-LIST', ''.join([name.ljust(16) for name in outnames]))
    setattr(out, 'NVARS', len(outnames))
    setattr(out, 'SDATE', int(dates[0].strftime('%Y%j')))
    setattr(out, 'STIME', int(dates[0].strftime('%H%M%S')))
    setattr(out, 'EDATE', int(dates[-1].strftime('%Y%j')))
    setattr(out, 'ETIME', int(dates[-1].strftime('%H%M%S')))
    for src, name, expr, outunit in config['mappings']:
        var = out.createVariable(name, 'f', vardims)
        var.long_name = name.ljust(16);
        var.var_desc = name.ljust(80);
        var.units = outunit.ljust(16)
    out.lonlatcoords = coordstr
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
                    thisdate = [datetime.strptime('%d %06d' % (d, t), '%Y%j %H%M%S') for d, t in tflag[:, 0]]
                    tname = 'TSTEP'
                elif 'time' in grp.variables.keys():
                    time = grp.variables['time']
                    unit, starttxt = grp.variables['time'].units.split(' since ')
                    starttxt = starttxt.replace('UTC', '').strip()
                    start_date = datetime.strptime(starttxt, '%Y-%m-%d %H:%M:%S')
                    thisdate = [start_date + timedelta(**{unit: t}) for t in time[:]]
                    tname = 'time'
                elif 'tau0' in grp.variables.keys():
                    time = grp.variables['tau0']
                    unit, starttxt = time.units.split(' since ')
                    start_date = datetime.strptime(start_date)
                    thisdate = [start_date + timedelta(**{unit: t}) for t in time]
                    tname = 'time'
                elif isinstance(f, (icon_profile, bcon_profile, profile_to_bdy, profile_to_grid)):
                    if len(dates) > 1:
                        for k in grp.variables.keys():
                            v = grp.variables[k]
                            if v.shape[0] == 1:
                                grp.variables[k] = v.repeat(len(dates), 0)
                    return grp
                else:
                    return grp
                thisdate = np.array(thisdate)
                idx, = np.where(np.logical_and(thisdate >= dates[0], thisdate <= dates[-1]))
                grpt = slice_dim(grp, '%s,%d,%d' % (tname, idx[0], idx[-1] + 1))
                return grpt
                
    else:
        raise KeyError('No file provided has the %s group' % src)

class file_getter(object):
    def __init__(self, config, out, sources, verbose = 0):
        """
        Initialize a file getting object that caches old files
        for reuse when necessary.
          config   - configuration dictionary
          out      - output file that has lonlatcoords, VGLVLS, VGTOP
          sources  - dictionary (key = groups, values = variables) that
                     identify which groups and variables will be used 
                    from GEOS-Chem will be used
          verbose  - True to see more details
        """
        self._config = config
        self._out = out
        self._sources = sources
        self._verbose = verbose
        self.last_file_paths = None
        self.last_file_objs = None
        if out is None:
            self._coordstr = None
        else:
            self._coordstr = out.lonlatcoords
        self.last_file_paths = [''] * len(config['file_templates'])
        self.last_file_objs = [PseudoNetCDFFile()] * len(config['file_templates'])
        if out is None:
            self.vert_out = None
        else:
            self._pref = 101325.
            self.vert_out = get_vert_in(out, pref = self._pref, vgtop = out.VGTOP, calcgeospress = config['interpolation']['calcgeospress'])

    def get_files(self, date):
        """
        Put date into file_templates and return open files,
        where all files have been vertically interpolated
        and only variables that will be used are present
    
          date     - date to use for files
        """
    
        from PseudoNetCDF.sci_var import extract, slice_dim, getvarpnc
        from PseudoNetCDF.cmaqfiles.profile import bcon_profile, icon_profile
        import gc
        
        # make quick references to instance variables
        sources = self._sources
        verbose = self._verbose
        coordstr = self._coordstr
        met_vert_out = vert_out = self.vert_out
        met_vaxis = vaxis = 0
        
        # Fill file path templates with date using the time
        # function provided
        file_paths = [(r, eval(tsf)(date, p)) for r, p, tsf in self._config['file_templates']]
        
        # Return cached files if appropriate
        if file_paths == self.last_file_paths:
            return self.last_file_objs


        # If coordstr is not none, this is
        # a data extraction call and the status should be
        # updated
        if coordstr is not None:
            status('')
            status('-' * 40)
            status("Getting files for " + str(date))
            status('-' * 40)
            status('')
            
        # For each file, use the reader (r) to open the
        # path (p)
        for fi, (r, p) in enumerate(file_paths):
            # If last path is this path
            # no need to update
            lp = self.last_file_paths[fi]
            nf = self.last_file_objs[fi]
            if p != lp:
                if verbose > 0: timeit('GET_FILE %s' % p, True)

                # Close old file to prevent memory leaks
                nf.close()
                
                if verbose > 1: status('Opening %s with %s' % (p, r), show = False)
                onf = nf = eval(r)(p)
                
                if coordstr is None:
                    # Coordstr is None, so this call is just for some
                    # meta-data
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        nf = getvarpnc(onf, 'time TFLAG tau0 tau1 latitude longitude latitude_bounds longitude_bounds'.split())
                    if 'PERIM' in onf.dimensions.keys() and not isinstance(onf, bcon_profile):
                        nf.createDimension('PERIM', len(onf.dimensions['PERIM']))
                        nf.createDimension('LAY', len(onf.dimensions['LAY']))
                        nf.NCOLS = onf.NCOLS
                        nf.NROWS = onf.NROWS
                else:
                    # If coordstr is not None, this is a real data
                    # call
                    #
                    # Real data calls require vertical interpolation
                    
                    ## Calculate the vertical coordinate of
                    ## the input file
                    vert_in = get_vert_in(nf, pref = self._pref, vgtop = self._out.VGTOP)

                    needsvinterp = not np.all([vert_out == vert_in])
                    weights = get_interp_w(vert_in, vert_out, axis = vaxis)
                    if not self._config['interpolation']['extrapolate']:
                        weights = np.ma.min([weights, np.ones_like(weights)], axis = 0)
                        weights = np.ma.max([weights, np.zeros_like(weights)], axis = 0)
                    if isinstance(nf, (bcon_profile, icon_profile)):
                        if needsvinterp:
                            nf = interpvars(nf, weights, dimension = 'sigma-mid')
                        metf = [f for f in self.last_file_objs if 'PERIM' in f.dimensions][0]
                        ## profile files need to be converted
                        ## to METBDY coordinates by repeating
                        ## boundaries
                        nf = profile_to_ftype(nf, ncols = metf.NCOLS, nrows = metf.NROWS, ftype = metf.FTYPE)
                    elif isinstance(nf, bpch):
                        ## Only extract groups that are used
                        ## in mappings, and only extract variables
                        ## in those groups that are used
                        grpkeys = set(nf.groups.keys())
                        srcgrps = set(sources.keys()).intersection(grpkeys)
                        nonsrcgrps = grpkeys.difference(sources)
                        if needsvinterp and self._config['interpolation']['calcgeospress']:
                            vert_in = get_vert_in(nf, pref = self._pref, vgtop = self._out.VGTOP, calcgeospress = self._config['interpolation']['calcgeospress'])
                            tmpnf = getvarpnc(nf, vert_in.dimensions)
                            var = tmpnf.createVariable('PRES', vert_in.dtype.char, vert_in.dimensions)
                            var[:] = vert_in
                            tmpnf = extract(tmpnf, [coordstr], unique = False)
                            vert_in = tmpnf.variables['PRES']
                            weights = get_interp_w(vert_in, met_vert_out, axis = met_vaxis)
                            if not self._config['interpolation']['extrapolate']:
                                weights = np.ma.min([weights, np.ones_like(weights)], axis = 0)
                                weights = np.ma.max([weights, np.zeros_like(weights)], axis = 0)
                        for grpk in srcgrps:
                            grp = nf.groups[grpk]
                            grp = nf.groups[grpk] = getvarpnc(grp, sources[grpk])
                            
                            # If needs interpolation, first extract unique points
                            # then interpolate, then repeat unique points for
                            # all boundary points
                            if self._config['interpolation']['calcgeospress']:
                                grp = nf.groups[grpk] = extract(grp, [coordstr])
                            else:
                                grp = nf.groups[grpk] = extract(grp, [coordstr], unique = True)
                            if needsvinterp:
                                # GEOS-Chem files can have multiple layer
                                # dimensions that share pressure coordinates
                                # but may have fewer output points
                                #
                                # for each dimensions, interpolation must 
                                # be done separately
                                for dimk in grp.dimensions:
                                    if dimk[:5] == 'layer':
                                        nlays = len(grp.dimensions[dimk])
                                        grp = nf.groups[grpk] = interpvars(grp, weights = weights.take(range(nlays), met_vaxis + 1), dimension = dimk, loginterp = (grp.variables.keys() if self._config['interpolation']['log'] else []))
                            
                            # Extract values that are appropriate
                            # for each boundary grid cell
                            if not self._config['interpolation']['calcgeospress']:
                                nf.groups[grpk] = extract(grp, [coordstr], gridded = False)

                        for grpk in nonsrcgrps:
                            # Unused groups are explicitly
                            # removed
                            del nf.groups[grpk]
                        
                    elif isinstance(nf, (Dataset, NetCDFFile, METBDY3D, METCRO3D, PseudoNetCDFFile)):
                        ## Dataset and NetCDFFile could have groups,
                        ## but that is not treated at this time.
                        
                        # If needs interpolation, first extract unique points
                        # then interpolate, then repeat unique points for
                        # all boundary points
                        if isinstance(nf, (METBDY3D, METCRO3D)):
                            tmpnf = getvarpnc(nf, ['PRES'])
                            tmpnf = convolve_dim(tmpnf, 'TSTEP,valid,0.5,0.5')
                            self.met_vert_out = met_vert_out = tmpnf.variables['PRES'][:]
                            met_vaxis = list(tmpnf.variables['PRES'].dimensions).index('LAY')
                        else:
                            nf = extract(nf, [coordstr], unique = True)
                            if needsvinterp:
                                if 'LAY' in nf.dimensions:
                                    nf = interpvars(nf, weights = weights, dimension = 'LAY')
                                else:
                                    nf = interpvars(nf, weights = weights, dimension = 'layer47')
                            nf = extract(nf, [coordstr], gridded = False)
                    else:
                        raise IOError('Unknown type %s; add type to readers' % type(nf))
                
                if verbose > 0: timeit('GET_FILE %s' % p, False)

            self.last_file_objs[fi] = nf
        self.last_file_paths = file_paths
        return self.last_file_objs

def get_dates(config):
    """
    Get all dates using, start_date, end_date and time_incr
    
    start_date - first date object
    end_date - last date object
    time_incr - timedelta object to increment
    """
    end_date = config['end_date']
    time_incr = config['time_incr']
    alldates = []
    date = config['start_date']
    if time_incr.total_seconds() == 0:
        return [date]
    
    while date <= end_date:
        alldates.append(date)
        date = date + time_incr
    return alldates

def pres_from_sigma(sigma, pref, ptop, avg = False):
    """
    Calculates pressure from sigma coordinates
        pres = sigma * (pref - ptop) + ptop
    where
        sigma = a sigma coordinate
        pref  = surface reference pressure    
        ptop  = the top of the model
        avg   = True if should be an average from edges
    """
    pres = sigma * (pref - ptop) + ptop
    if avg:
        pres = pres[:-1] + np.diff(pres) / 2.
    return pres
    
def get_vert_in(nf, pref, vgtop, calcgeospress = False):
    """
    Calculates a pressure coordinate
    from sigma coordinates in a netcdf-like
    file
    """
    #if hasattr(nf, 'IOAPI_VERSION'):
    #    if 'PRES' in nf.variables.keys():
    #        vert_in = nf.variables['PRES'][:]
    #    elif hasattr(nf, 'VGLVLS'):
    #        vert_in = pres_from_sigma(nf.VGLVLS, pref, vgtop, avg = True)
    if hasattr(nf, 'VGLVLS'):
        vert_in = pres_from_sigma(nf.VGLVLS, pref, vgtop, avg = True)
    elif 'sigma-mid' in nf.variables.keys():
        sigma = np.array(nf.variables['sigma-mid'])
        # assuming VGTOP of model is consistent with
        # VGTOP of profile
        vert_in = pres_from_sigma(sigma, pref, vgtop)
    elif calcgeospress:
        if 'PEDGE-\$_PSURF' in nf.variables.keys():
            vert_in = nf.variables['PEDGE-\$_PSURF'][:]
        elif 'TIME-SER_AIRDEN' in nf.variables.keys() and 'DAO-3D-$_TMPU' in nf.variables.keys():
            warn('Calculating pressure from Air Density and Temperature')
            airden = nf.variables['TIME-SER_AIRDEN'][:] * 1e6
            temperature = nf.variables['DAO-3D-$_TMPU'][:]
            from scipy.constants import Avogadro, R
            vert_in = (airden / Avogadro * temperature * R)
        else:
            raise ValueError('Could not find geos pressure; disable calculation of pressure.')
    elif 'layer' in nf.variables.keys():
        warn('Using pressure from website')
        if 'PEDGE-\$_PSURF' in nf.variables.keys():
            warn('Pressure is available; enable calcgeospress to utilize.')
        elif 'TIME-SER_AIRDEN' in nf.variables.keys() and 'DAO-3D-$_TMPU' in nf.variables.keys():
            warn('AIRDEN and TMPU are available; enable calcgeospress to utilize.')
        layerv = nf.variables['layer']
        punit = layerv.units.strip()
        if punit in ('millibar', 'hPa'):
            vert_in = layerv[:] * 100.
        elif punit == 'Pa':
            vert_in = layerv[:]
        elif punit == 'model layer':
            warn('Interpreting bpch model layer unit as hPa; deprecated and improved in updated PseudoNetCDF.')
            vert_in = layerv[:] * 100.
        else:
            raise ValueError('Not sure how to process unit %s; need millibar, hPa or Pa' % punit)
            
    else:
        raise KeyError('No sigma-mid or layer in %s' % str(nf.variables.keys()))
    return vert_in

def profile_to_ftype(nf, ncols, nrows, ftype):
    if ftype == 2:
        return profile_to_bdy(nf, ncols, nrows)
    elif ftype == 1:
        return profile_to_grid(nf, ncols, nrows)
    else:
        raise ValueError('ftype must be 2 (BNDY) or 1 (CRO3D); got %d' % ftype)

class profile_to_bdy(PseudoNetCDFFile):
    def __init__(self, nf, ncols, nrows):
        """
        nf - profile file (netcdf-like)
        ncols - number of columns
        nrows - number of rows
        """
        outf = self
        p2n = Pseudo2NetCDF()
        p2n.addGlobalProperties(nf, outf)
        outf.createDimension('sigma-mid', len(nf.dimensions['sigma-mid']))
        outf.createDimension('TSTEP', 1)
        outf.createDimension('PERIM', (ncols + nrows + 2) * 2)
        for k, v in nf.variables.items():
            if k in ('sigma', 'sigma-mid'): continue
            nv = np.concatenate([v[:, [0]].repeat(ncols + 1, 1), v[:, [1]].repeat(nrows + 1, 1), v[:, [2]].repeat(ncols + 1, 1), v[:, [3]].repeat(nrows + 1, 1)], axis = 1)
            vard = dict([(pk, getattr(v, pk)) for pk in v.ncattrs()])
            outf.createVariable(k, v.dtype.char, ('time', 'sigma-mid', 'PERIM'), values = nv[None,:], kgpermole = 1., **vard)
        outf.groups = dict(PROFILE = outf)

class profile_to_grid(PseudoNetCDFFile):
    def __init__(self, nf, ncols, nrows):
        """
        nf - profile file (netcdf-like)
        ncols - number of columns
        nrows - number of rows
        """
    
        outf = self
        p2n = Pseudo2NetCDF()
        p2n.addGlobalProperties(nf, outf)
        outf.createDimension('sigma-mid', len(nf.dimensions['sigma-mid']))
        outf.createDimension('ROW', nrows)
        outf.createDimension('COL', ncols)
        outf.createDimension('TSTEP', 1)
    
        j, i = np.indices([nrows, ncols])
        north_fraction = j.astype('d') / j.max()
        south_fraction = 1 - north_fraction
        east_fraction = i.astype('d') / i.max()
        west_fraction = 1 - east_fraction
        outval = np.zeros((1, len(nf.dimensions['sigma-mid']), nrows, ncols), dtype = 'f')
        for k, v in nf.variables.items():
            if k in ('sigma', 'sigma-mid'): continue
            vard = dict([(pk, getattr(v, pk)) for pk in v.ncattrs()])
            # If this is an icon file
            if v.ndim == 2 and v.shape[1] == 1:
                v = v.repeat(4, 1)[:, :, None, None]
            else:
                v = v[:, :, None, None]
            outval[0, :] = (north_fraction * v[:, 2] + south_fraction * v[:, 0] + west_fraction * v[:, 3] + east_fraction * v[:, 1]) / 2
            outf.createVariable(k, v.dtype.char, ('time', 'sigma-mid', 'ROW', 'COL'), values = outval.astype(v.dtype), kgpermole = 1., **vard)
        outf.groups = dict(PROFILE = outf)
