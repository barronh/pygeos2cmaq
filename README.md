Overview
========

This is intended to drastically simplify GEOS2CMAQ by relying on PseudoNetCDF

User Steps:
-----------

1. Create a template mapping file (python -m pygeos2cmaq -t > config.yaml
2. Edit paths and start/stop date
3. Run pygeos2cmaq (python -m pygeos2cmaq config.yaml

Tutorial:

1. Install the software (see Installation)
2. ``python -m pygeos2cmaq -t > config.yaml``
3. Copy the testdata folder to the same folder where config.yaml is
4. python -m pygeos2cmaq config.yaml


Program Steps:
--------------

1. Open a METBDY file
    * METBDY provides coordinates for horizontal and vertical interpolation
    * latitude/longitude are calculated by PseudoNetCDF
    * VGLVLS stores Sigma Levels by default
2. Open a GEOS-Chem file
    * latitude/longitude (and bounds) are calculated by PseudoNetCDF
    * sigma-eta levels are provided by PseudoNetCDF
3. Find GEOS-Chem cells that belong to METBDY cells
    * Only cells from GEOS-Chem that overlap with METBDY are selected
    * Only variables that will be used are selected (including AIRDEN)
4. Algebra is applied to map GEOS-Chem species to CMAQ
5. A "profile.dat" file fills in the gaps.
6. Vertical interpolation is applied


Install Instructions
--------------------

1. Windows
  1. Install python-xy using Windows executable
  2. Open DOS command prompt (start run -> type "cmd" -> hit enter
  3. type "pip install https://github.com/barronh/pygeos2cmaq/archive/master.zip"
2. Mac OSX
  1. Install the Anaconda Python Distribution
  2. Install following "Get Pip" instructions at http://www.pip-installer.org/en/latest/installing.html
  3. Open a Terminal
  4. type "sudo pip install https://github.com/barronh/pygeos2cmaq/archive/master.zip"
5. Ubuntu Linux
  1. Open Terminal
  2. sudo apt-get install python-scipy
  3. sudo apt-get install python-pip
  4. sudo pip install http://cheeseshop.python.org/packages/source/p/pyparsing/pyparsing-1.5.5.tar.gz
  5. type "sudo pip install https://github.com/barronh/pygeos2cmaq/archive/master.zip"
