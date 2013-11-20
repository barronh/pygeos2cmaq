template_string = """
start_date: '2010-01-01 23:00:00'
end_date:   '2010-01-02 01:00:00'
time_incr: 1 hours
out_template: 'geos2cmaq.%%Y%%m%%d.nc'

unitconversions:
    ppbv->ppmV: <value>/1000.
    ppbv->micrograms/m**3: <value> * AIRMOLPERM3
    molec/cm3->ppmV: <value> * CM3PERMMOLEC
    metadefs:
        AIRMOLPERM3: 1e6 * AIRDEN / 6.022E+23
        CM3PERMMOLEC: 1e6 /  AIRDEN
# PPB2NMOL [=] mole/m3; nmole/mole * mole/m3 = nmole/m3
# MCC2PPM [=] cm3/Mmolecule; molecule/cm3 * cm3/Mmolecule = mole/Mmol = umol / mol

file_templates:
  - [bpch, 'testdata/ts%%Y%%m%%d.bpch', minus1hour]
  - [ioapi, 'testdata/METBDY3D_100101', simpledate]
#  - [profile, 'testdata/profile.dat', simpledate]
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