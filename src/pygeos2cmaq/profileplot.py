import sys
import pylab as pl
import numpy as np
from warnings import warn
from netCDF4 import MFDataset
from collections import defaultdict
from datetime import datetime, timedelta

unitconvert = {('ppmV', 'ppb'): lambda x: x * 1000.}

def pres_from_sigma(sigma, pref, ptop, avg = False):
    pres = sigma * (pref - ptop) + ptop
    if avg:
        pres = pres[:-1] + np.diff(pres) / 2.
    return pres

def plottes(ax, lon_bnds, lat_bnds, tespaths):
    from netCDF4 import Dataset
    from glob import glob
    allx = []
    ally = []
    lon_bnds = lon_bnds[:]
    lat_bnds = lat_bnds[:]
    tespaths = reduce(list.__add__, [glob(i) for i in tespaths])
    for path in tespaths:
        f = Dataset(path)
        lats = f.variables['latitude'][:][:, None]
        lons = f.variables['longitude'][:][:, None]
        pressure = f.variables['pressure']
        species = f.variables['species']
        x = []
        y = []
        inlon = np.logical_and(lons >= lon_bnds[None, :, 0], lons <= lon_bnds[None, :, 1])
        inlat = np.logical_and(lats >= lat_bnds[None, :, 0], lats <= lat_bnds[None, :, 1])
        inboth = np.logical_and(inlon, inlat).any(1)
        if inboth.sum(0) > 0:            
            print '******** FOUND ******', path
            x = pressure[inboth]
            y = species[inboth]
            allx.append(x)
            ally.append(y)
        else:
            warn('No data found for %s' % path)
    
    if len(allx) == 0:
        return            
    var = np.ma.masked_values(ally, -999.) * 1e9
    var = var.reshape(-1, var.shape[-1])
    vertcrd = np.ma.masked_values(allx, -999.).mean(0).mean(0)
    minval = var.swapaxes(0, 1).reshape(var.shape[1], -1).min(1)
    meanval = var[:].swapaxes(0, 1).reshape(var.shape[1], -1).mean(1)
    maxval = var[:].swapaxes(0, 1).reshape(var.shape[1], -1).max(1)
    tesline = ax.plot(meanval, vertcrd, ls = '-', lw = 2, color = 'r', label = r'TES', zorder = 3)

    x = np.ma.concatenate([minval[:vertcrd.size], maxval[:vertcrd.size][::-1]])
    y = np.ma.concatenate([vertcrd[:], vertcrd[::-1]])
    mask = x.mask | y.mask
    x = np.ma.masked_where(mask, x).compressed()
    y = np.ma.masked_where(mask, y).compressed()
    tesrange = ax.fill(x, y, facecolor = 'r', edgecolor = 'r', alpha = .7, zorder = 1, ls = 'solid', label = 'TES min/max')

def plot(paths, keys = ['O3'], prefix = 'BC', scale = 'log', minmax = (None, None), minmaxq = (0, 100), outunit = 'ppb', sigma = False, maskzeros = False, tespaths = [], edges = True):
    from pylab import figure, NullFormatter, close, rcParams
    rcParams['text.usetex'] = False
    from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm, LogNorm
    f = MFDataset(paths)
    try:
        if sigma:
            vertcrd = f.VGLVLS[:-1] + np.diff(f.VGLVLS)
        else:
            vertcrd = pres_from_sigma(f.VGLVLS, pref = 101325., ptop = f.VGTOP, avg = True) / 100.
    except:
        lay = f.variables['layer_edges'][:]
        if not sigma:
            vertcrd = lay[:-1] + np.diff(lay) / 2.
        else:
            vertcrd = (lay[:] - lay[-1]) / (lay[0] - lay[-1])
            vertcrd = vertcrd[:-1] + np.diff(vertcrd) / 2

    geosls = '-'
    geoscolor = 'k'
    geosrangecolor = 'k'
    geosrangeecolor = 'k'
    geosrangels = 'solid'
    alpha = .7
    lonb = f.variables['longitude_bounds']
    latb = f.variables['latitude_bounds']
    for var_name in keys:
        temp = defaultdict(lambda: 1)
        try:
            eval(var_name, None, temp)
            var = eval(var_name, None, f.variables)[:]
        except:
            temp[var_name]
            var = f.variables[var_name][:]
        if maskzeros: var = np.ma.masked_values(var, 0)
        unit = f.variables[temp.keys()[0]].units.strip()
        var = unitconvert.get((unit, outunit), lambda x: x)(var)
        bmap = None
        vmin, vmax = np.percentile(np.ma.compressed(var).ravel(), list(minmaxq))
        if minmax[0] is not None:
            vmin = minmax[0]
        if minmax[1] is not None:
            vmax = minmax[1]
        if edges:
            fig = pl.figure(figsize = (16, 4))
            offset = 0.05
            ax = fig.add_axes([.1 - offset, .15, .225, .725])
            ax = fig.add_axes([.325 - offset, .15, .225, .725])
            ax = fig.add_axes([.55 - offset, .15, .225, .725])
            ax = fig.add_axes([.775 - offset, .15, .225, .725])
            ss = 0
            se = ss + f.NCOLS + 1
            es = se
            ee = se + f.NROWS + 1
            ns = ee
            ne = ee + f.NCOLS + 1
            ws = ne
            we = ws + f.NROWS + 1
            axs = fig.axes
            for ax in fig.axes[1:]:
                ax.yaxis.set_major_formatter(pl.NullFormatter())
            
            vars = [var[:, :, ss:se], var[:, :, es:ee], var[:, :, ns:ne][:, :, ::-1], var[:, :, ws:we][:, :, ::-1]]
            lonbss = [lonb[ss:se], lonb[es:ee], lonb[ns:ne][::-1], lonb[ws:we][::-1]]
            latbss = [latb[ss:se], latb[es:ee], latb[ns:ne][::-1], latb[ws:we][::-1]]
            
        else:
            fig = pl.figure(figsize = (8, 4))
            ax = fig.add_axes([.1, .15, .8, .725])
            axs = fig.axes
            vars = [var]
            lonbss = [lonb[:]]
            latbss = [latb[:]]
        for ax, var, lonbs, latbs in zip(axs, vars, lonbss, latbss):
            minval = var.swapaxes(0, 1).reshape(var.shape[1], -1).min(1)
            meanval = var[:].swapaxes(0, 1).reshape(var.shape[1], -1).mean(1)
            maxval = var[:].swapaxes(0, 1).reshape(var.shape[1], -1).max(1)
            modline = ax.plot(meanval, vertcrd, ls = geosls, lw = 2, color = geoscolor, label = r'GC', zorder = 4)

            x = np.ma.concatenate([minval[:vertcrd.size], maxval[:vertcrd.size][::-1]])
            y = np.ma.concatenate([vertcrd[:], vertcrd[::-1]])
            mask = x.mask | y.mask
            x = np.ma.masked_where(mask, x).compressed()
            y = np.ma.masked_where(mask, y).compressed()
            modrange = ax.fill(x, y, facecolor = geosrangecolor, edgecolor = geosrangeecolor, alpha = alpha, zorder = 1, ls = geosrangels, label = 'GC min/max')
            ymin, ymax = vertcrd.min(), vertcrd.max()
            ax.set_ylim(ymax, ymin)
            ax.set_xscale(scale)
            ax.set_xlim(vmin, vmax)
            #if scale == 'log':
            #    ax.set_xticklabels(['%.1f' % (10**x) for x in ax.get_xticks()])
            try:
                sdate = datetime.strptime('%d %06d' % (f.SDATE, f.STIME), '%Y%j %H%M%S')
                edate = datetime.strptime('%d %06d' % (f.EDATE, f.ETIME), '%Y%j %H%M%S')
            
            except Exception, e:
                print str(e)
                sdate = datetime(1985, 1, 1, 0) + timedelta(hours = f.variables['tau0'][0])
                edate = datetime(1985, 1, 1, 0) + timedelta(hours = f.variables['tau1'][-1])
            if len(tespaths) > 0:
                plottes(ax, lonbs, latbs, tespaths)

        try:
            title = '%s to %s' % (sdate.strftime('%Y-%m-%d'), edate.strftime('%Y-%m-%d'))
        except:
            title = var_name
        if sigma:
            axs[0].set_ylabel('sigma')
        else:
            axs[0].set_ylabel('pressure')

        xmax = -np.inf
        xmin = np.inf
        for ax in fig.axes:
            tmp_xmin, tmp_xmax = ax.get_xlim()
            xmax = max(tmp_xmax, xmax)
            xmin = min(tmp_xmin, xmin)
        for ax in fig.axes:
            ax.set_xlim(xmin, xmax)
            
        if len(axs) == 1:
            axs[0].set_xlabel('%s %s' % (var_name, outunit))
        else:
            axs[0].set_xlabel('South')
            axs[1].set_xlabel('East')
            axs[2].set_xlabel('North')
            axs[3].set_xlabel('West')
        pl.legend(bbox_to_anchor = (.5, 1), loc = 'upper center', bbox_transform = fig.transFigure, ncol = 4)
        if edges:
            fig.text(0.95, 0.975, title, horizontalalignment = 'right', verticalalignment = "top", fontsize = 16)
        else:
            fig.text(0.95, 0.025, title, horizontalalignment = 'right', verticalalignment = "bottom", fontsize = 16)
        fig.savefig('%s_%s_%s.png' % (prefix, var_name, 'profile'))
        pl.close(fig)
    return fig

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser("Usage: %prog [options] ifile\n\n\tifile - netcdf input file with variables; layer dimension must be second")

    parser.add_option("-v", "--variables", dest = "variables", action = "append", default = [],
                        help = "Variable name or expression; multiple names can be provided by specifying -v multiple times")
    
    parser.add_option("-p", "--prefix", dest = "prefix", type = "string", default = None,
                        help = "Prefix for figures")

    parser.add_option("-n", "--no-map", dest = "nomap", action = "store_true", default = False,
                        help = "Try to plot with map")

    parser.add_option("", "--sigma", dest = "sigma", action = "store_true", default = False,
                        help = "Plot on sigma coordinate instead of pressure")

    parser.add_option("-s", "--scale", dest = "scale", type = "string", default = 'log',
                        help = "Defaults to deciles (i.e., 10 equal probability bins), but linear and log are also options.")

    parser.add_option("", "--minmax", dest = "minmax", type = "string", default = "None,None",
                        help = "Defaults None, None.")

    parser.add_option("", "--mask-zeros", dest = "maskzeros", action = "store_true", default = False,
                        help = "Defaults False.")

    parser.add_option("", "--minmaxq", dest = "minmaxq", type = "string", default = '0,100',
                        help = "Defaults 0,100.")

    parser.add_option("", "--out-unit", dest = "outunit", type = "string", default = 'ppb',
                        help = "Defaults ppb.")

    parser.add_option("", "--tes-paths", dest = "tespaths", type = "string", default = [], action = "append",
                        help = "Plot tes on top of boundary from paths; defaults to []")

    parser.add_option("-f", "--time-func", dest = "timefunc", default = "mean",
                        help = "Use time-func to reduce the time dimension (mean, min, max, std, var, ndarray.__iter__, etc.")
                        
    parser.add_option("-e", "--edges", dest = "edges", default = False, action = "store_true",
                        help = "Plot S,E,N,W")

    (options, args) = parser.parse_args()
    
    if not len(args) > 0:
        parser.print_help()
        exit()
    if options.prefix is None:
        options.prefix = args[0]
    if len(options.variables) == 0:
        options.variables = ['O3']
    if len(options.tespaths) > 0:
        options.tespaths = reduce(list.__add__, [tp.split(',') for tp in options.tespaths])
    fig = plot(args, keys = options.variables, prefix = options.prefix, scale = options.scale, minmax = eval(options.minmax), minmaxq = eval(options.minmaxq), sigma = options.sigma, maskzeros = options.maskzeros, outunit = options.outunit, tespaths = options.tespaths, edges = options.edges)
