import sys
import pylab as pl
import numpy as np
from warnings import warn
from netCDF4 import MFDataset, Dataset
from datetime import datetime, timedelta


def aqmeiidomain():
    from mpl_toolkits.basemap import Basemap
    # From griddesc
    aqmeii_proj = {}
    aqmeii_proj['llcrnrx'] = -2556000.
    aqmeii_proj['llcrnry'] = -1728000.
    aqmeii_proj['dx'] = 12000.0
    aqmeii_proj['dy'] = 12000.0
    aqmeii_proj['nx'] = 459
    aqmeii_proj['ny'] = 299

    # Derived
    exec('width = nx * dx', None, aqmeii_proj)
    exec('height = ny * dy', None, aqmeii_proj)
    exec('urcrnrx = llcrnrx + width', None, aqmeii_proj)
    exec('urcrnry = llcrnry + height', None, aqmeii_proj)

    cmaqmap = Basemap(rsphere = (6370000., 6370000.),\
                        resolution = 'c', projection = 'lcc',\
                        lat_1 = 33., lat_2 = 45., lat_0 = 40., lon_0 = -97.,\
                        llcrnrx = aqmeii_proj['llcrnrx'], llcrnry = aqmeii_proj['llcrnry'],\
                        urcrnrx = aqmeii_proj['urcrnrx'], urcrnry = aqmeii_proj['urcrnry'])
    return cmaqmap

def domainfromcmaq(f):
    from mpl_toolkits.basemap import Basemap
    # From griddesc
    aqmeii_proj = {}
    aqmeii_proj['llcrnrx'] = f.XORIG
    aqmeii_proj['llcrnry'] = f.YORIG
    aqmeii_proj['dx'] = f.XCELL
    aqmeii_proj['dy'] = f.YCELL
    aqmeii_proj['nx'] = f.NCOLS
    aqmeii_proj['ny'] = f.NROWS

    # Derived
    exec('width = nx * dx', None, aqmeii_proj)
    exec('height = ny * dy', None, aqmeii_proj)
    exec('urcrnrx = llcrnrx + width', None, aqmeii_proj)
    exec('urcrnry = llcrnry + height', None, aqmeii_proj)

    cmaqmap = Basemap(rsphere = (6370000., 6370000.),\
                        resolution = 'c', projection = 'lcc',\
                        lat_1 = f.P_ALP, lat_2 = f.P_BET, lat_0 = f.YCENT, lon_0 = f.XCENT,\
                        llcrnrx = aqmeii_proj['llcrnrx'], llcrnry = aqmeii_proj['llcrnry'],\
                        urcrnrx = aqmeii_proj['urcrnrx'], urcrnry = aqmeii_proj['urcrnry'])
    return cmaqmap


def pres_from_sigma(sigma, pref, ptop, avg = False):
    pres = sigma * (pref - ptop) + ptop
    if avg:
        pres = pres[:-1] + np.diff(pres) / 2.
    return pres
    
def plot(paths, keys = ['O3'], func = 'mean', map = True, prefix = 'BC', scale = 'deciles', minmax = (None, None), minmaxq = (0, 100), sigma = False, maskzeros = False):
    from pylab import figure, NullFormatter, close, rcParams
    rcParams['text.usetex'] = False
    from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm, LogNorm
    if len(paths) == 1:
        f = Dataset(paths[0])
    else:
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
    reversevert = not (np.diff(vertcrd) > 0).all()
    for var_name in keys:
        var = eval(var_name, None, f.variables)[:]
        if func == 'each':
            vars = [('time%02d' % vi, v) for vi, v in enumerate(var)]
        else:
            vars = [(func, getattr(np, func)(var, axis = 0))]
        for func, var in vars:
            bmap = None
            if maskzeros: var = np.ma.masked_values(var, 0)
            vmin, vmax = np.percentile(np.ma.compressed(var).ravel(), list(minmaxq))
            if minmax[0] is not None:
                vmin = minmax[0]
            if minmax[1] is not None:
                vmax = minmax[1]
            if scale == 'log':
                bins = np.logspace(np.log10(vmin), np.log10(vmax), 11)
            elif scale == 'linear':
                bins = np.linspace(vmin, vmax, 11)
            elif scale == 'deciles':
                bins = np.percentile(np.ma.compressed(np.ma.masked_greater(np.ma.masked_less(var, vmin), vmax)).ravel(), [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
                bins[0] = vmin; bins[-1] = vmax
            norm = BoundaryNorm(bins, ncolors = 256)
            if map:
                fig = pl.figure(figsize = (8, 8))
                fig.subplots_adjust(hspace = .3, wspace = .3)
                axmap = fig.add_subplot(3,3,5)
                try:
                    cmaqmap = domainfromcmaq(f)
                    cmaqmap.drawcoastlines()
                    cmaqmap.drawcountries()
                    cmaqmap.drawstates()
                    cmaqmap.drawparallels(np.arange(-90, 100, 10), labels = [True, True, False, False])
                    cmaqmap.drawmeridians(np.arange(-180, 190, 20), labels = [False, False, True, True])
                except Exception, e:
                    warn('An error occurred and no map will be shown:\n%s' % str(e))
                axn = fig.add_subplot(3,3,2, sharex = axmap)
                axw = fig.add_subplot(3,3,4, sharey = axmap)
                axe = fig.add_subplot(3,3,6, sharey = axmap)
                axs = fig.add_subplot(3,3,8, sharex = axmap)
                cax = fig.add_axes([.8, .7, .05, .25])
                for ax in [axmap, axe]:
                    ax.yaxis.set_major_formatter(NullFormatter())
                for ax in [axmap, axn]:
                    ax.xaxis.set_major_formatter(NullFormatter())
                for ax in [axn, axs]:
                    if sigma:
                        ax.set_ylabel('sigma')
                    else:
                        ax.set_ylabel('pressure')
                for ax in [axe, axw]:
                    if sigma:
                        ax.set_xlabel('sigma')
                    else:
                        ax.set_xlabel('pressure')
                xyfactor = 1
            else:
                fig = pl.figure(figsize = (16, 4))
                fig.subplots_adjust(bottom=0.15)
                axw = fig.add_subplot(1,4,1)
                axn = fig.add_subplot(1,4,2)
                axe = fig.add_subplot(1,4,3)
                axs = fig.add_subplot(1,4,4)
                cax = fig.add_axes([.91, .1, .025, .8])
                if sigma:
                    axw.set_ylabel('sigma')
                else:
                    axw.set_ylabel('pressure')
            
                xyfactor = 1e-3 # m -> km
                     
            x = f.NCOLS + 1
            y = f.NROWS + 1
            start_south = 0
            end_south = start_south + x
            start_east = end_south
            end_east = start_east + y
            start_north = end_east
            end_north = start_north + x
            start_west = end_north
            end_west = start_west + y
            X, Y = np.meshgrid(np.arange(x) * f.XCELL * xyfactor, vertcrd)
            patchess = axs.pcolor(X, Y, var[:, start_south:end_south], cmap = bmap, vmin = vmin, vmax = vmax, norm = norm)
            if not map:
                if reversevert: axs.set_ylim(*axs.get_ylim()[::-1])
                axs.set_title('South')
                axs.set_xlabel('E to W km')
                axs.set_xlim(*axs.get_xlim()[::-1])
                
            X, Y = np.meshgrid(np.arange(-1, x - 1) * f.XCELL * xyfactor, vertcrd)
            patchesn = axn.pcolor(X, Y, var[:, start_north:end_north], cmap = bmap, vmin = vmin, vmax = vmax, norm = norm)
            if reversevert: axn.set_ylim(*axn.get_ylim()[::-1])
            if not map:
                axn.set_title('North')
                axn.set_xlabel('W to E km')

            if map:
                X, Y = np.meshgrid(vertcrd, np.arange(y) * f.YCELL)
                patchese = axe.pcolor(X, Y, var[:, start_east:end_east].swapaxes(0,1), cmap = bmap, vmin = vmin, vmax = vmax, norm = norm)
                if reversevert: axe.set_xlim(*axe.get_xlim()[::-1])
            else:
                X, Y = np.meshgrid(np.arange(y) * f.YCELL * xyfactor, vertcrd)
                patchese = axe.pcolor(X, Y, var[:, start_east:end_east], cmap = bmap, vmin = vmin, vmax = vmax, norm = norm)
                if reversevert: axe.set_ylim(*axe.get_ylim()[::-1])
                axe.set_title('East')
                axe.set_xlabel('N to S km')
                axe.set_xlim(*axe.get_xlim()[::-1])
            if map:
                X, Y = np.meshgrid(vertcrd, np.arange(-1, y - 1) * f.YCELL)
                patchesw = axw.pcolor(X, Y, var[:, start_west:end_west].swapaxes(0,1), cmap = bmap, vmin = vmin, vmax = vmax, norm = norm)
            else:
                X, Y = np.meshgrid(np.arange(-1, y - 1) * f.YCELL * xyfactor, vertcrd)
                patchesw = axw.pcolor(X, Y, var[:, start_west:end_west], cmap = bmap, vmin = vmin, vmax = vmax, norm = norm)
                if reversevert: axw.set_ylim(*axw.get_ylim()[::-1])
                axw.set_title('West')
                axw.set_xlabel('S to N km')
            if map:
                for ax in [axe, axw]:
                    ax.axis('tight', axis = 'x')
                    pl.setp( ax.xaxis.get_majorticklabels(), rotation=90 )
                for ax in [axs, axn]:
                    ax.axis('tight', axis = 'y')
            else:
                for ax in [axe, axn, axw, axs] + ([axmap] if map else []):
                    ax.axis('tight')
            
            if 'TFLAG' in f.variables.keys():
                SDATE = f.variables['TFLAG'][:][0, 0, 0]
                EDATE = f.variables['TFLAG'][:][-1, 0, 0]
                STIME = f.variables['TFLAG'][:][0, 0, 1]
                ETIME = f.variables['TFLAG'][:][-1, 0, 1]
                if SDATE == 0:
                    SDATE = 1900001
                    EDATE = 1900001
                sdate = datetime.strptime('%07d %06d' % (SDATE, STIME), '%Y%j %H%M%S')
                edate = datetime.strptime('%07d %06d' % (EDATE, ETIME), '%Y%j %H%M%S')
            elif 'tau0' in f.variables.keys():
                sdate = datetime(1985, 1, 1, 0) + timedelta(hours = f.variables['tau0'][0])
                edate = datetime(1985, 1, 1, 0) + timedelta(hours = f.variables['tau1'][-1])
            else:
                sdate = datetime(1900, 1, 1, 0)
                edate = datetime(1900, 1, 1, 0)
            try:
                title = '%s %s to %s' % (var_name, sdate.strftime('%Y-%m-%d'), edate.strftime('%Y-%m-%d'))
            except:
                title = var_name
            fig.suptitle(title.replace('O3', 'Ozone at Regional Boundaries'))
            fig.colorbar(patchesw, cax = cax, boundaries = bins)
            cax.set_xlabel('ppm')
            fig.savefig('%s_%s_%s.png' % (prefix, var_name, func))
            pl.close(fig)
    
if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser("Usage: %prog [options] ifile\n\n\tifile - netcdf input file with variables; time dimension must be first, and layer dimension must be second")

    parser.add_option("-v", "--variables", dest = "variables", action = "append", default = [],
                        help = "Variable name or expression; multiple names can be provided by specifying -v multiple times")

    parser.add_option("-p", "--prefix", dest = "prefix", type = "string", default = None,
                        help = "Prefix for figures")

    parser.add_option("-n", "--no-map", dest = "nomap", action = "store_true", default = False,
                        help = "Try to plot with map")

    parser.add_option("", "--sigma", dest = "sigma", action = "store_true", default = False,
                        help = "Plot data on sigma coordinate.")

    parser.add_option("-s", "--scale", dest = "scale", type = "string", default = 'deciles',
                        help = "Defaults to deciles (i.e., 10 equal probability bins), but linear and log are also options.")

    parser.add_option("", "--mask-zeros", dest = "maskzeros", action = "store_true", default = False,
                        help = "Defaults False.")

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
    plot(args, keys = options.variables, map = not options.nomap, prefix = options.prefix, func = options.timefunc, scale = options.scale, minmax = eval(options.minmax), minmaxq = eval(options.minmaxq), sigma = options.sigma, maskzeros = options.maskzeros)
