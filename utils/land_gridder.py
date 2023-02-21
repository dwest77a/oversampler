## Graph Creator for L2 files
## Version 1.0 - 10/05/2021 09:15
## Daniel Westwood - daniel.westwood@stfc.ac.uk

import os
import sys

import shapely.affinity
from shapely.geometry import Point, Polygon
import geopy
import geopy.distance

import time
import numpy as np
import numpy.ma as ma

from netCDF4 import Dataset

#from mpl_toolkits import mplot3d
#from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap
#from matplotlib.patches import Ellipse, Rectangle

from getopt import getopt

import math

# Standard dw package tools
STD_PY_LOC = '/home/users/dwest77/Documents/std_py'
try:
    sys.path.append(STD_PY_LOC)
    import file_import as fm
    import find_files as ff
    import pmath as pm
    import output_data as od
except: 
    print('ImportError: Std_py library missing: requires file_import.py, find_files.py and pmath.py as minimum')
    sys.exit()
g = 9.81
map_bounds = [0,0,0,0]
bound = False
isfilter1=True
isfilter2=False
isfilter3=True
isfilter4=True
istpwsec = True

mmair = 28.964001
mmh2o = 18
g = 9.80665
dw = 1000
PI = 3.1415926


map_example = Basemap(projection='cyl',llcrnrlat=-89.5, urcrnrlat=89.5, 
									   llcrnrlon=-179.5, urcrnrlon=179.5, lat_ts=20, resolution='l')

def determine_land(lat, lon):
	xpt, ypt = map_example(lon/1000, lat/1000)
	return int(map_example.is_land(xpt, ypt))
	
def main():
	land_flag = []
	latlist = [-89.5 + i/10 for i in range(0,1790)]
	lonlist = [-179.5 + i/10 for i in range(0,3580)]
	for ilat in range(0,1790):
		print(ilat)
		row = []
		for ilon in range(0,3580):
			
			lat = -89.5 + ilat/10
			lon = -179.5 + ilon/10
			
			flag = determine_land(lat, lon)
			row.append(flag)
		land_flag.append(row)
	print(land_flag)
	
	ncf_new = Dataset('land_grid.nc', 'w', format='NETCDF4')
	
	
	lat = ncf_new.createDimension('lat',len(latlist))
	lon = ncf_new.createDimension('lon',len(lonlist))
		
	lat = ncf_new.createVariable('lat', np.float32, ('lat',))
	lat.long_name = 'latitude'
	lat.units = 'deg'
	lat[:] = latlist
		
	lon = ncf_new.createVariable('lon', np.float32, ('lon',))
	lon.long_name = 'longitude'
	lon.units = 'deg'
	lon[:] = lonlist
	
	land_flag_var = ncf_new.createVariable('land_flag', np.float32, ('lat','lon',))
	land_flag_var.long_name = 'Land Flag (0.1x0.1 ref grid)'
	land_flag_var[:] = land_flag
	
	ncf_new.close()
		
main()
