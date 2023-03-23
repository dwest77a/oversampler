"""
L3 oversampled gridding script
--------------
Takes oversampled daily files and averages to monthly files on a regular grid

"""
__author__	= "Daniel Westwood"
__date__	  = "23-03-2023"
__copyright__ = "Copyright 2020 United Kingdom Research and Innovation"

import os
import sys
from datetime import datetime
from getopt import getopt

import numpy as np
from netCDF4 import Dataset

# Standard pyanalysis package tools
from pyanalysis import fileIO as ff
from pyanalysis import datasetMath as dm
    
    
# Combine daily files into monthly gridded files

def create_grid(res, map_bounds):
	mbs = np.array(map_bounds) / res
	data = np.zeros((int(mbs[1]-mbs[0]), int(mbs[3]-mbs[2])))
	return data
	
def create_grid_3d(res, map_bounds):
	mbs = np.array(map_bounds) / res
	grid = np.zeros((int(mbs[1]-mbs[0]), int(mbs[3]-mbs[2]), 30))
	grid[:] = np.nan
	return grid
	
class var_grid:
	def __init__(self, map_bounds, res, vnames, dn):
		self.map_bounds = map_bounds
		self.res = res
		
		self.dn = dn
		
		self.data = []
		self.count = []
		self.grid_squares = []
		
		self.vnames = vnames
		
		for vname in vnames:
			self.data.append(create_grid(res, map_bounds))
			self.grid_squares.append(create_grid_3d(res, map_bounds))
		self.count = create_grid(res, map_bounds)
		
		self.latlist = [map_bounds[0] + i*res for i in range(len(self.data[0]))]
		self.lonlist = [map_bounds[2] + i*res for i in range(len(self.data[0][0]))]
		
		self.data_bounds = [0, len(self.data[0]), 0, len(self.data[0][0])]
		
		lat_range = map_bounds[1] - map_bounds[0]
		lon_range = map_bounds[3] - map_bounds[2]
		
		data_lat_range = self.data_bounds[1] - self.data_bounds[0] 
		data_lon_range = self.data_bounds[3] - self.data_bounds[2]
		
		self.latratio = data_lat_range/lat_range
		self.lonratio = data_lon_range/lon_range
	
	def add_file(self, filename):
		ncf = Dataset(filename,'r',format='NETCDF4')
		
		tpw_filter = np.array(ncf['tpw']) < 40
		nh3 = np.array(ncf['nh3'])
		nh3_filter = np.logical_not(np.isnan(nh3))
		#nh3_f1 = nh3 < 0.8
		#h3_f2 = nh3 > -0.5
		
		full_filter = tpw_filter & nh3_filter #& nh3_f1 & nh3_f2
		
		def pr(arr):
			print(arr.size, np.count_nonzero(arr.astype('int')))
		
		latlist = np.array(ncf['lat'])[full_filter]
		lonlist = np.array(ncf['lon'])[full_filter]
		
		varlist = []
		for vname in self.vnames:
			varlist.append(np.array(ncf[vname])[full_filter])
		
		if len(latlist) > 0:
		
		
			dlatlist = np.array(self.data_bounds[0] + (latlist-self.map_bounds[0])*self.latratio)
			dlonlist = np.array(self.data_bounds[2] + (lonlist-self.map_bounds[2])*self.lonratio)
			
			dlatlist[dlatlist > len(self.data[0])-1] = len(self.data[0])-1
			dlonlist[dlonlist > len(self.data[0][0])-1] = len(self.data[0][0])-1
			
			
			dlatlist = dlatlist.astype('int')
			dlonlist = dlonlist.astype('int')
			for v in range(len(self.vnames)):
				self.data[v][dlatlist, dlonlist] += varlist[v]
				
			self.count[dlatlist,dlonlist] += 1
			return True
	
		else:
			print('Date skipped')
			return False
		
	def average_grid(self):
		
		for v in range(len(self.vnames)):
			self.data[v][:,:] = self.data[v][:,:]/self.count[:,:]
		
	def save_data(self, fname):
		print(np.count_nonzero(self.data[0]))
		os.system('touch {}'.format(fname))
		ncf_new = Dataset(fname,'w',format='NETCDF4')
		
		lat = ncf_new.createDimension('lat',len(self.latlist))
		lon = ncf_new.createDimension('lon',len(self.lonlist))
		count = 0
		
		lat = ncf_new.createVariable('lat', np.float32, ('lat',))
		lat.long_name = 'latitude'
		lat.units = 'deg'
		lat[:] = self.latlist
		
		lon = ncf_new.createVariable('lon', np.float32, ('lon',))
		lon.long_name = 'longitude'
		lon.units = 'deg'
		lon[:] = self.lonlist
		
		count = ncf_new.createVariable('count',np.float32, ('lat','lon',))
		count[:,:] = self.count
		
		
		nh3 = ncf_new.createVariable('nh3',np.float32, ('lat','lon',))
		nh3.long_name = 'Ammonia Concentration Mean'
		nh3.units = 'ppbv'
		nh3[:,:] = self.data[0]
		
		print(np.nanmax(self.data[0]), np.nanmin(self.data[0]))
		
		ncf_new.close()



options, operands = getopt(sys.argv[1:], "", ["year=","month=","dn=","instrt=","version=","subversion="])

year = operands[0]
month = operands[1]
dn = operands[2]
instrt = operands[3]
version = operands[4]
subversion = operands[5]

instrt_long = instrt

fpath = '/gws/pw/j05/rsg_share/public/transfer/barry/6ad76b04-b3df-11eb-a4d3-024ad8e814ad/{}_l2/oversampled_l2/{}/{}/{}/'.format(instrt_long, version, year, month)
fname = '{}_os_reference_l2_{}{}{}_{}_{}.nc'.format(instrt, year, month, dn, version, subversion) 
CLOBBER = True
if CLOBBER or not os.path.isfile(fpath+fname):

	filename = '/gws/pw/j05/rsg_share/public/transfer/barry/6ad76b04-b3df-11eb-a4d3-024ad8e814ad/{}_l2/valid_l2a/{}/{}/{}/{}_valid_l2a_{}{}{}_{}.nc'
	map_bounds = [-89.5,89.5,-179.5,179.5]
	res = 0.25
	nh3_map = var_grid(map_bounds, res, ['nh3'], dn)

	not_empty = False
	for d in range(1,31):
		day = format(d,'02d')
		print(day,'/31')
		try:
			e = nh3_map.add_file(filename.format(instrt_long, version, year, month, instrt, year, month, day, version))
			not_empty = not_empty or e
		except FileNotFoundError:
			print('skipped',year, month, day)
	if not_empty:
		if not os.path.isdir(fpath):
			os.makedirs(fpath)
		if not os.path.isfile(fpath+fname):
			os.system('touch {}{}'.format(fpath, fname))
		nh3_map.average_grid()
		nh3_map.save_data(fpath+fname)
		print('saved to ',fpath+fname)
else:
	print('Skipping existing file')

