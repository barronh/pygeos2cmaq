Overview
========

This is intended to drastically simplify GEOS2CMAQ by relying on PseudoNetCDF

User Steps:

1. Create a template mapping file (python -m pygeos2cmaq -t > config.yaml
2. Edit paths and start/stop date
3. Run pygeos2cmaq (python -m pygeos2cmaq config.yaml

Program Steps:

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


5. Vertical interpolation is applied
