from optparse import OptionParser
from glob import glob
import os

def run():
    parser = OptionParser()
    parser.set_usage("Usage: %prog [-tq] \n"+(" "*16)+" [-i <init name>] [-f <final name>] <yamlfile>")
    parser.add_option("-t", "--template", dest = "template", action = "store_true", default = False, help="Output template on standard out (configurable with -m and -c", metavar="Template")

    parser.add_option("-v", "--verbose", dest = "verbose", action = "store_true", default = False, help = "extra output for debugging", metavar = "VERBOSE")
    
    paths = glob(os.path.join(os.path.dirname(__file__), 'mapping', '*_*.yaml'))
    mechanisms = ', '.join(['_'.join(path.split('/')[-1].split('_')[1:])[:-5] for path in paths])
    parser.add_option("-c", "--configuration", dest="configuration", default="%s/mapping/cb05cl_ae6_aq.csv" % os.path.dirname(__file__),
                        help="Chemical mechanisms: %s (for use with -t)" % mechanisms, metavar="CONFIG")

    (options, args) = parser.parse_args()
    if options.template:
        from template import template
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