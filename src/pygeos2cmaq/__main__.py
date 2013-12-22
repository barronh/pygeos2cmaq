from optparse import OptionParser
from glob import glob
import os
from warnings import warn

def run():
    defaultmech = "%s/mapping/cb05cl_ae6_aq.csv" % os.path.dirname(__file__)
    parser = OptionParser()
    parser.set_usage("Usage: %prog [-tq] \n"+(" "*16)+" [-i <init name>] [-f <final name>] <yamlfile>")
    parser.add_option("-t", "--template", dest = "template", action = "store_true", default = False, help="Output template on standard out (configurable with -m and -c", metavar="Template")

    parser.add_option("-v", "--verbose", dest = "verbose", action = "count", default = 0, help = "extra output for debugging", metavar = "VERBOSE")
    
    paths = glob(os.path.join(os.path.dirname(__file__), 'mapping', '*_*.csv'))
    mechanisms = ', '.join(['_'.join(path.split('/')[-1].split('_')[:])[:-4] for path in paths])
    parser.add_option("-c", "--configuration", dest="configuration", default = None,
                        help = "Chemical mechanisms: %s (for use with -t)" % mechanisms, metavar="CONFIG")

    (options, args) = parser.parse_args()
    if options.template:
        from template import template
        if options.configuration is None:
            warn("Using default mechanism: %s" % defaultmech)
        else:
            if os.path.exists(options.configuration):
                pass
            else:
                options.configuration = "%s/mapping/%s.csv" % (os.path.dirname(__file__), options.configuration)
                if not os.path.exists(options.configuration):
                    raise IOError('Cannot find file %s; must be either you own file or in %s' % (options.configuration, mechanisms))
        print template(options.configuration)
        parser.exit()
    if len(args)<1:
        parser.error(msg="Requires a yaml file as an argument.  For a template use the -t option.  The template will be output to the stdout.")
    else:
        yamlpath=args[0]
        from load import loader
        from process import process
        outf = process(config = loader(yamlpath), verbose = options.verbose)

if __name__ == '__main__':
    run()