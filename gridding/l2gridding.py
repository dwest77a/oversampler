import os
import sys
from datetime import datetime
from getopt import getopt

import numpy as np
from netCDF4 import Dataset
    
    
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
		
		spatial1 = np.array(ncf['lat']) > map_bounds[0]
		spatial2 = np.array(ncf['lat']) < map_bounds[1]
		spatial3 = np.array(ncf['lon']) > map_bounds[2]
		spatial4 = np.array(ncf['lon']) < map_bounds[3]
		tempfilter = np.array(ncf['dt1000']) > 10
		tpwfilter = np.array(ncf['tpw']) < 40
		
		print(np.array(ncf['tpw_base']))
		
		# day_night filter
		if self.dn == 'day':
			dnfilter = np.array(ncf['day_flag']) == 1
		else:
			dnfilter = np.array(ncf['day_flag']) != 1
		
		
		# surface elevation filter? - apply later
		
		
		# Additional filters
		# tpw > 40 (with sec)
		# surface elevation < 450m
		# dt1000 > 15
		
		all_space = spatial1 & spatial2 & spatial3 & spatial4 & tpwfilter & dnfilter #& tempfilter
		
		def pr(arr):
			print(arr.size, np.count_nonzero(arr.astype('int')))
		
		latlist = np.array(ncf['lat'])[all_space]
		lonlist = np.array(ncf['lon'])[all_space]
		varlist = []
		for vname in self.vnames:
			
			if vname == 'nh3_adj':
				nh3_adj = np.array(ncf['nh3']) - 0.003*np.array(ncf['tpw']) + 0.06 # IASI L2 night sea
				varlist.append(nh3_adj[all_space])
			else:
				varlist.append(np.array(ncf[vname])[all_space])
		
		if len(latlist) > 0:
		
		
			dlatlist = np.array(self.data_bounds[0] + (latlist-self.map_bounds[0])*self.latratio)
			dlonlist = np.array(self.data_bounds[2] + (lonlist-self.map_bounds[2])*self.lonratio)
			
			
			dlatlist = dlatlist.astype('int')
			
			dlonlist = dlonlist.astype('int')
			for v in range(len(self.vnames)):
				self.data[v][dlatlist, dlonlist] += varlist[v]
				
				for i in range(len(dlatlist)):
					ilat = dlatlist[i]
					jlon = dlonlist[i]
					index_next = 0
					while not np.isnan(self.grid_squares[v][ilat, jlon, index_next]) and index_next < len(self.grid_squares[v][ilat, jlon])-1:
						index_next += 1
					self.grid_squares[v][ilat,jlon, index_next] = varlist[v][i]
				
			self.count[dlatlist,dlonlist] += 1
			return True
	
		else:
			print('Date skipped')
			return False
		
	def average_grid(self):
		
		for v in range(len(self.vnames)):
			self.data[v][:,:] = self.data[v][:,:]/self.count[:,:]
			
		# Standard deviation:
			self.grid_squares[v] = np.nanstd(self.grid_squares[v], axis=2)
			
		
	def save_data(self, fname):
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
		nh3.long_name = 'Ammonia Concentration'
		nh3.units = 'ppbv'
		nh3[:,:] = self.data[0]
		
		nh3_sd = ncf_new.createVariable('nh3_sd',np.float32, ('lat','lon',))
		nh3_sd.long_name = 'Ammonia Concentration Standard Deviation'
		nh3_sd.units = 'ppbv'
		nh3_sd[:,:] = self.grid_squares[0]
		
		tpw_base = ncf_new.createVariable('tpw_base',np.float32, ('lat','lon',))
		tpw_base.long_name = 'Total Precipitable Water Column LOS'
		tpw_base.units = 'mm'
		tpw_base[:,:] = self.data[1]
		
		tpw_base_sd = ncf_new.createVariable('tpw_base_sd',np.float32, ('lat','lon',))
		tpw_base_sd.long_name = 'Total Precipitable Water Column LOS Standard Deviation'
		tpw_base_sd.units = 'mm'
		tpw_base_sd[:,:] = self.grid_squares[1]
		
		tpw = ncf_new.createVariable('tpw',np.float32, ('lat','lon',))
		tpw.long_name = 'Total Precipitable Water with adjusted zenith angle'
		tpw.units = 'mm'
		tpw[:,:] = self.data[2]
		
		tpw_sd = ncf_new.createVariable('tpw_sd',np.float32, ('lat','lon',))
		tpw_sd.long_name = 'Total Precipitable Water with adjusted zenith angle Standard Deviation'
		tpw_sd.units = 'mm'
		tpw_sd[:,:] = self.grid_squares[2]
		
		dt1000 = ncf_new.createVariable('dt1000',np.float32, ('lat','lon',))
		dt1000.long_name = 'Surface Temperature to 1km difference'
		dt1000.units = 'K'
		dt1000[:,:] = self.data[3]
		
		dt1000_sd = ncf_new.createVariable('dt1000_sd',np.float32, ('lat','lon',))
		dt1000_sd.long_name = 'Surface Temperature to 1km difference Standard Deviation'
		dt1000_sd.units = 'K'
		dt1000_sd[:,:] = self.grid_squares[3]
		
		nh3_adj = ncf_new.createVariable('nh3_tpwvza', np.float32, ('lat','lon',))
		nh3_adj.long_name = 'Ammonia Concentration with TPW VZA adjustment'
		nh3_adj.units = 'ppbv'
		nh3_adj[:,:] = self.data[4]
		
		nh3_adj_sd = ncf_new.createVariable('nh3_tpwvza_sd', np.float32, ('lat','lon',))
		nh3_adj_sd.long_name = 'Ammonia Concentration with TPW VZA adjustment Standard Deviation'
		nh3_adj_sd.units = 'ppbv'
		nh3_adj_sd[:,:] = self.grid_squares[4]
		
		ncf_new.close()



options, operands = getopt(sys.argv[1:], "", ["year=","month=","dn=","instrt=","version=","subversion="])

year = operands[0]
month = operands[1]
dn = operands[2]
instrt = operands[3]
version = operands[4]
subversion = operands[5]

instrt_long = instrt

fpath = '/gws/pw/j05/rsg_share/public/transfer/barry/6ad76b04-b3df-11eb-a4d3-024ad8e814ad/{}_l2/gridded_l2/{}/{}/{}/'.format(instrt_long, version, year, month)
fname = '{}_gridded_l2_{}{}{}_{}_{}.nc'.format(instrt, year, month, dn, version, subversion) 
if True:#not os.path.isfile(fpath+fname):

	filename = '/gws/pw/j05/rsg_share/public/transfer/barry/6ad76b04-b3df-11eb-a4d3-024ad8e814ad/{}_l2/valid_l2a/{}/{}/{}/{}_valid_l2a_{}{}{}_{}.nc'
	map_bounds = [-89.5,89.5,-179.5,179.5]
	res = 0.25
	nh3_map = var_grid(map_bounds, res, ['nh3','tpw_base','tpw','dt1000', 'nh3_adj'], dn)

	not_empty = False
	for d in range(1,31):
		day = format(d,'02d')
		if True:
			e = nh3_map.add_file(filename.format(instrt_long, version, year, month, instrt, year, month, day, version))
			not_empty = not_empty or e
		#except:
			#print('skipped',year, month, day)
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

