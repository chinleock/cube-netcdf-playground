# Cube/ NetCDF Playground

Space Time Cube, or Cube, is developed by Esri Spatial Statistics Team, trying to extend the analytics to space and time dimension. Cube is built based on NetCDF4-HDF5.

This Github repository is about some experiments to discover what can Cube do, and what's the capability of NetCDF4-Python can achieve.

## Required Packages
Numpy with MKL

[NetCDF4-Python API](http://netcdf4-python.googlecode.com/svn/trunk/docs/netCDF4-module.html)

## Optional Packages
...Depands on each script files

[ArcGIS Desktop](http://www.esri.com/software/arcgis/arcgis-for-desktop)

[Pymc](https://pypi.python.org/pypi/pymc)

## Scripts
* Experiment_cube_size_by_dimension.py: 

   An experiment to test how dimension size affect netcdf file size, with the same amount of data

* Utilities_netcdf_cube_convertor.py: 

   An utility tool to help convert an netCDF file with 3-Dim data to a Cube format, which can be used in Esri [Space Time Pattern Mining Toolbox](http://desktop.arcgis.com/en/desktop/latest/tools/space-time-pattern-mining-toolbox/an-overview-of-the-space-time-pattern-mining-toolbox.htm)
Still onging, an Esri PYT toolbox will be create when finalize