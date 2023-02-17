import os
import sys
from datetime import datetime
from getopt import getopt

import numpy as np
from netCDF4 import Dataset

# Standard dw package tools
STD_PY_LOC = '/home/users/dwest77/Documents/std_py'
try:
    sys.path.append(STD_PY_LOC)
    import file_import as fm
    import find_files as ff
    import pmath as pm
except: 
    print('ImportError: Std_py library missing: requires file_import.py, find_files.py and pmath.py as minimum')
    sys.exit()
    
    
# Combine daily files into monthly gridded files

def create_grid(res, map_bounds):
	mbs = np.array(map_bounds) / res
	data = np.zeros((int(mbs[1]-mbs[0]), int(mbs[3]-mbs[2])))
	return data
	
class var_grid:
	def __init__(self, map_bounds, res, vnames, dn):
		self.map_bounds = map_bounds
		self.res = res
		
		self.dn = dn
		
		self.data = []
		self.count = []
		
		self.vnames = vnames
		
		self.elev = create_grid(res, map_bounds)
		
		self.count = create_grid(res, map_bounds)
		
		self.latlist = [map_bounds[0] + i*res for i in range(len(self.elev))]
		self.lonlist = [map_bounds[2] + i*res for i in range(len(self.elev[0]))]
		
		self.data_bounds = [0, len(self.elev), 0, len(self.elev[0])]
		
		lat_range = map_bounds[1] - map_bounds[0]
		lon_range = map_bounds[3] - map_bounds[2]
		
		data_lat_range = self.data_bounds[1]
		data_lon_range = self.data_bounds[3]
		
		self.latratio = data_lat_range/lat_range
		self.lonratio = data_lon_range/lon_range
	
	def add_file(self, filename):
		ncf = Dataset(filename,'r',format='NETCDF4')
		print(ncf.variables)
		for var in ncf.variables:
			print(var, ncf.variables[var].long_name)
		x=input()
		
		latlist = np.array(ncf['lat'])
		lonlist = np.array(ncf['lon'])
		elev = np.array(ncf[''])
		
		dlatlist = np.array(self.data_bounds[0] + latlist*self.latratio)
		dlonlist = np.array(self.data_bounds[2] + lonlist*self.lonratio)
			
		dlatlist = dlatlist.astype('int')
		dlonlist = dlonlist.astype('int')
		
		for v in range(len(self.vnames)):
			self.elev[dlatlist, dlonlist] += elev
		self.count[dlatlist,dlonlist] += 1
		return True
		
	def average_grid(self):
		
		self.elev[:,:] = self.elev[:,:]/self.count[:,:]
		
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
		count[:,:] = self.count[0]
		
		elev = ncf_new.createVariable('surf_elev',np.float32, ('lat','lon',))
		elev.long_name = 'Surface Elevation from IASI...'
		elev.units = 'm'
		elev[:,:] = self.elev
		
		
		
		ncf_new.close()


options, operands = getopt(sys.argv[1:], "", ["year=","month=","dn=","version=","subversion="])

year = operands[0]
month = operands[1]
dn = operands[2]
version = operands[3]
subversion = operands[4]

instrt = 'iasi'
instrt_long = 'iasi_mhs'

filename = '/gws/pw/j05/rsg_share/public/transfer/barry/6ad76b04-b3df-11eb-a4d3-024ad8e814ad/{}_l2/valid_l2a/v00.03/{}/{}/{}_valid_l2a_{}{}{}_{}.nc'
map_bounds = [-89.5,89.5,-179.5,179.5]
res = 0.25
nh3_map = var_grid(map_bounds, res, ['nh3','tpw_base','tpw'], dn)
not_empty = False
for d in range(1,31):
	day = format(d,'02d')
	if True:
		e = nh3_map.add_file(filename.format(instrt_long, year, month, instrt, year, month, day, version))
		print(d, e)
		not_empty = not_empty or e
	#except:
		#print('skipped',year, month, day)
if not_empty:
	fpath = '/gws/pw/j05/rsg_share/public/transfer/barry/6ad76b04-b3df-11eb-a4d3-024ad8e814ad/{}_l2/gridded_l2/{}/{}/{}/'.format(instrt_long, version, year, month)
	fname = '{}_gridded_l2_{}{}{}_{}_{}.nc'.format(instrt, year, month, dn, version, subversion)
	if not os.path.isdir(fpath):
		os.makedirs(fpath)
	if not os.path.isfile(fpath+fname):
		os.system('touch {}{}'.format(fpath, fname))
	nh3_map.average_grid()
	nh3_map.save_data(fpath+fname)
	print('saved to ',fpath+fname)

