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
3. Ubuntu Linux
  1. Open Terminal
  2. sudo apt-get install python-scipy
  3. sudo apt-get install python-pip
  4. sudo pip install http://cheeseshop.python.org/packages/source/p/pyparsing/pyparsing-1.5.5.tar.gz
  5. type "sudo pip install https://github.com/barronh/pygeos2cmaq/archive/master.zip"
4. Linux without Permissions
  1. Navigate in terminal to the folder you want to install/work in.
  2. type `curl -kLO http://repo.continuum.io/miniconda/Miniconda-3.5.5-Linux-x86_64.sh`
    1. path may change in the future
  3. type `bash Miniconda-3.5.5-Linux-x86_64.sh`
    1. confirm prompts as necessary until you get to the install path (at the time of this writing: hit enter; hit space; type yes; hit enter;)
    2. Choose `./miniconda` as the install path
    3. Decide whether or not all future sessions should use this python
      1. Choose yes to enable this miniconda for all future sessions;
      2. Choose no to require miniconda to be enabled as needed;
  4. Enable miniconda (this step is required for at least this session): `export PATH=${PWD}/miniconda/bin:${PATH}`
  5. Install packages (type code in each step and confirm prompts)
    1. `conda install numpy`
    2. `conda install scipy`
    3. `conda install matplotlib`
    4. `conda install basemap`
    5. `conda install netcdf4`
    6. `conda install pip`
    7. `pip install --no-deps --upgrade --install-option="--prefix=${PWD}/miniconda/" git+https://code.google.com/p/pseudonetcdf/`
    8. `pip install --no-deps --upgrade --install-option="--prefix=${PWD}/miniconda/" https://github.com/barronh/pygeos2cmaq/archive/master.zip`
  6. You're ready to go!
  



