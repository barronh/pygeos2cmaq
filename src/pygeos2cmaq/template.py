template_string = """
# Start and End date, inclusive, will be output
# in the boundary condition file using the time
# increment descibed by time_incr.
#
# start_date and end_date in the form YYYY-MM-DD HH:MM:SS
# time_incr in the form val units (units must be plural, e.g., hours, seconds)
start_date: '2010-01-01 23:00:00'
end_date:   '2010-01-02 01:00:00'
time_incr: 1 hours

# Output files are stored in files matching
# the out_template
out_template: 'geos2cmaq.%%Y%%m%%d.nc'

# Files used by pygeos2cmaq will be found based on the patterns
# below after converting current date using a function and
# then the file will be read using an interpreter.
#
# [interpreter, pathtemplate, function]
# - interpreter is a file reading function (e.g., NetCDFFile, ioapi, bpch, cspec)
# - pathtemplate will be processed by a strftime like function to create a path
# - function is a function that uses strftime and the current date to process the path template
file_templates:
  - [bpch, 'testdata/ts%%Y%%m%%d.bpch', minus1hour]
  - [ioapi, 'testdata/METBDY3D_100101', simpledate]
# Not currently supporting PROFILES
#  - [profile, 'testdata/profile.dat', simpledate]


# GEOS-Chem species will automatically convert from
#    ppbC to ppbv
#    ppbv to micrograms/mol
#    ppbC to micrograms/mol
#
# Some non-automatic, but default unit conversions
# are supplied below
# 
# Each unit conversion has a key and an expression
# - key tells the from and to unit (from->to)
# - expression will be evaluated with <value> replaced
#   with the variable value
unitconversions:
    ppbv->ppmV: <value>/1000.
    micrograms/mol->micrograms/m**3: <value> * AIRMOLPERM3
    molec/cm3->ppmV: <value> * CM3PERMMOLEC
    metadefs:
        AIRMOLPERM3: 1e6 * AIRDEN / 6.022E+23
        CM3PERMMOLEC: 1e6 /  AIRDEN

# Mappings provide algebraic processing of GEOS-Chem
# variables to make CMAQ variables, with unit processing
# as described above.
mappings:
%s

"""

def template(config = 'mapping/cb05cl_ae6_aq.csv'):
    """
    Create a template by loading mappings from config
    and adding them to a template string
    """
    import numpy as np
    
    # Get mappings from CSV file
    mappings = np.recfromtxt(config, delimiter = ',', names = True, comments = '+|$')
    
    # Initialize lines list
    lines = []
    
    # Add yaml formatting
    for maps in mappings:
        line = '  - %s' % map(lambda x: x.strip(), maps)
        lines.append(line)
    
    # Put mappings into template and return
    out = template_string % ('\n'.join(lines))
    return out

if __name__ == '__main__':
    print template()