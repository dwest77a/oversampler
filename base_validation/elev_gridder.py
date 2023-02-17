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

res = 0.125
map_bounds = [-89.5, 89.5, -179.5, 179.5]

class do_elev_grid:
	def __init__(self, res, map_bounds):
		self.res = res
		self.map_bounds = map_bounds
		
		self.create_grid()
		
	def create_grid(self):
		latrange = int((self.map_bounds[1] - self.map_bounds[0])/self.res)
		lonrange = int((self.map_bounds[3] - self.map_bounds[2])/self.res)
	
		self.data = np.repeat(0.0, latrange*lonrange)
		self.data = np.reshape(self.data, (latrange, lonrange))
	
		self.count = np.repeat(0.0, latrange*lonrange)
		self.count = np.reshape(self.count, (latrange, lonrange))
		
		self.latlist = np.linspace(self.map_bounds[0], self.map_bounds[1], num = int((self.map_bounds[1] - self.map_bounds[0])/res) )
		self.lonlist = np.linspace(self.map_bounds[2], self.map_bounds[3], num = int((self.map_bounds[3] - self.map_bounds[2])/res) )
	
	def do_file(self, fileN):
		ncf = Dataset(fileN, 'r', format='NETCDF4')
		elev = np.array(ncf['ecmwf_altitude'])
		lat = np.array(ncf['latitude'])
		lon = np.array(ncf['longitude'])
	
		abslat = lat -self.map_bounds[0]
		abslon = lon -self.map_bounds[2]
	
		dlat = abslat / self.res#*(self.map_bounds[1] - self.map_bounds[0])
		dlon = abslon / self.res#*(self.map_bounds[3] - self.map_bounds[2])
		
		dlat[dlat > len(self.data)-1] = len(self.data)-1
		dlon[dlon > len(self.data[0])-1] = len(self.data[0])-1
	
		dlat = dlat.astype('int')
		dlon = dlon.astype('int')
	
		self.data[dlat,dlon] += elev
		self.count[dlat,dlon] += 1
		
	def find_files(self):
		filepath = '/gws/pw/j05/rsg_share/public/transfer/barry/6ad76b04-b3df-11eb-a4d3-024ad8e814ad/cris_l2/valid_l2a/v00.00c/2016/'
		files_array = ff.list_files(filepath, starts='ral',ends='.nc')
		print('found {} files'.format(len(files_array)))
		x=input()
		index = 0
		is_20 = False
		while index < 2000:#en(files_array) and not is_20:
			print(index+1,'/',len(files_array))
			self.do_file(files_array[index])
			assess = (self.data.size - np.count_nonzero(self.data)) / self.data.size
			if assess < 0.2:
				is_20 = True
			print(assess, end='	')
			index += 1
			
			
		self.average_grid()
		self.save_to_nc()
		
	def average_grid(self):
		self.data = self.data / self.count
		
	def save_to_nc(self):
		
		ncf_new = Dataset('elev_grid_125.nc','w',format='NETCDF4')
		
		lat = ncf_new.createDimension('lat',len(self.latlist))
		lon = ncf_new.createDimension('lon',len(self.lonlist))
		
		lat = ncf_new.createVariable('lat', np.float32, ('lat',))
		lat.long_name = 'latitude'
		lat.units = 'deg'
		lat[:] = self.latlist
		
		lon = ncf_new.createVariable('lon', np.float32, ('lon',))
		lon.long_name = 'longitude'
		lon.units = 'deg'
		lon[:] = self.lonlist
		
		elev = ncf_new.createVariable('elev',np.float32, ('lat','lon',))
		elev.long_name = 'Estimated Surface Elevation'
		elev.units='m'
		elev[:,:] = self.data
		
		ncf_new.close()
		
EG = do_elev_grid(res, map_bounds)
EG.find_files()
