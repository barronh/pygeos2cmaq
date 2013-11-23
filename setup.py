from distutils.core import setup

setup(name = 'pygeos2cmaq',
      version = '0.9',
      author = 'Barron Henderson',
      author_email = 'barronh@gmail.com',
      maintainer = 'Barron Henderson',
      maintainer_email = 'barronh@gmail.com',
      url='https://github.com/barronh/pygeos2cmaq/',
      download_url='https://github.com/barronh/pygeos2cmaq/archive/master.zip',
      long_description="pygeos2cmaq creates boundary conditions from GEOS-Chem for CMAQ.",
      packages = ['pygeos2cmaq'],
      package_dir = {'pygeos2cmaq': 'src/pygeos2cmaq'},
      package_data = {'pygeos2cmaq': ['mapping/*.csv']},
      scripts = ['scripts/pygeos2cmaq'],
      install_requires = ['numpy', 'scipy', 'PseudoNetCDF', 'pyyaml']
      )
