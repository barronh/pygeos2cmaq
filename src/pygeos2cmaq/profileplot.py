import sys
import pylab as pl
import numpy as np
from warnings import warn
from netCDF4 import MFDataset
from collections import defaultdict

unitconvert = {('ppmV', 'ppb'): lambda x: x * 1000.}

def plot(paths, keys = ['O3'], prefix = 'BC', scale = 'log', minmax = (None, None), minmaxq = (0, 100), outunit = 'ppb'):
    from pylab import figure, NullFormatter, close, rcParams
    rcParams['text.usetex'] = False
    from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm, LogNorm
    f = MFDataset(paths)
    try:
        vertcrd = f.VGLVLS[:-1] + np.diff(f.VGLVLS)
    except:
        lay = f.variables['layer_edges'][:]
        vertcrd = (lay[:] - lay[-1]) / (lay[0] - lay[-1])
        vertcrd = vertcrd[:-1] + np.diff(vertcrd) / 2

    geosls = '-'
    geoscolor = 'k'
    geosrangecolor = 'k'
    geosrangeecolor = 'k'
    geosrangels = 'solid'
    alpha = .7
    print keys
    for var_name in keys:
            temp = defaultdict(lambda: 1)
            try:
                eval(var_name, None, temp)
                var = eval(var_name, None, f.variables)[:]
            except:
                temp[var_name]
                var = f.variables[var_name][:]
            
            unit = f.variables[temp.keys()[0]].units.strip()
            var = unitconvert.get((unit, outunit), lambda x: x)(var)
            bmap = None
            vmin, vmax = np.percentile(np.ma.compressed(var).ravel(), list(minmaxq))
            if minmax[0] is not None:
                vmin = minmax[0]
            if minmax[1] is not None:
                vmax = minmax[1]
            fig = pl.figure(figsize = (16, 4))
            ax = fig.add_subplot(1,1,1)
            minval = var.swapaxes(0, 1).reshape(var.shape[1], -1).min(1)
            meanval = var[:].swapaxes(0, 1).reshape(var.shape[1], -1).mean(1)
            maxval = var[:].swapaxes(0, 1).reshape(var.shape[1], -1).max(1)
            modline = ax.plot(meanval, vertcrd, ls = geosls, lw = 2, color = geoscolor, label = r'GC $\mathbf{\hat{y}}^{i,m}_t$', zorder = 4)

            x = np.ma.concatenate([minval[:vertcrd.size], maxval[:vertcrd.size][::-1]])
            y = np.ma.concatenate([vertcrd[:], vertcrd[::-1]])
            mask = x.mask | y.mask
            x = np.ma.masked_where(mask, x).compressed()
            y = np.ma.masked_where(mask, y).compressed()
            modrange = ax.fill(x, y, facecolor = geosrangecolor, edgecolor = geosrangeecolor, alpha = alpha, zorder = 1, ls = geosrangels, label = 'GC min/max')
            ax.set_ylim(*ax.get_ylim()[::-1])
            ax.set_xscale(scale)
            ax.set_xlim(vmin, vmax)
            #if scale == 'log':
            #    ax.set_xticklabels(['%.1f' % (10**x) for x in ax.get_xticks()])
            fig.savefig('%s_%s_%s.png' % (prefix, var_name, 'profile'))
            pl.close(fig)
    return fig
    
if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    parser.set_usage("""Usage: python -m geos2cmaq.plot [-v VAR1,VAR2] [-p prefix] ifile

    ifile - path to a file formatted as type -f
    
    """)

    parser.add_option("-v", "--variables", dest = "variables", action = "append", default = [],
                        help = "Variable names separated by ','")

    parser.add_option("-p", "--prefix", dest = "prefix", type = "string", default = None,
                        help = "Prefix for figures")

    parser.add_option("-n", "--no-map", dest = "nomap", action = "store_true", default = False,
                        help = "Try to plot with map")

    parser.add_option("-s", "--scale", dest = "scale", type = "string", default = 'log',
                        help = "Defaults to deciles (i.e., 10 equal probability bins), but linear and log are also options.")

    parser.add_option("", "--minmax", dest = "minmax", type = "string", default = "None,None",
                        help = "Defaults None, None.")

    parser.add_option("", "--minmaxq", dest = "minmaxq", type = "string", default = '0,100',
                        help = "Defaults 0,100.")

    parser.add_option("-f", "--time-func", dest = "timefunc", default = "mean",
                        help = "Use time-func to reduce the time dimension (mean, min, max, std, var, ndarray.__iter__, etc.")

    (options, args) = parser.parse_args()
    
    if not len(args) > 0:
        parser.print_help()
        exit()
    if options.prefix is None:
        options.prefix = args[0]
    if len(options.variables) == 0:
        options.variables = ['O3']
    fig = plot(args, keys = reduce(list.__add__, [v.split(',') for v in options.variables]), prefix = options.prefix, scale = options.scale, minmax = eval(options.minmax), minmaxq = eval(options.minmaxq))
