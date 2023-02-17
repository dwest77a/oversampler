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

from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Ellipse, Rectangle

from getopt import getopt

import math

sys.path.append('/home/users/dwest77/Documents/std_py')
import find_files as ff
import pmath as pm
import output_data as od

#BASE_INPUT = '/gws/pw/j05/rsg_share/public/projects/ims/data/lv2/output_ghg_cv9_v12_lam_fgsnwp_nobc_newbc_rbc_ram2_rem20-ncam_rcl5_raer6st_rmg_oap_cnv2_o3f113_csyn_swir2/npp/{}/{}/{}'
BASE_INPUT = '/gws/pw/j05/rsg_share/public/transfer/barry/6ad76b04-b3df-11eb-a4d3-024ad8e814ad/{}_l2/valid_l2a/{}/{}/{}'
BASE_OUTPUT = '/gws/pw/j05/rsg_share/public/transfer/barry/6ad76b04-b3df-11eb-a4d3-024ad8e814ad/{}_l2/oversampled_l2/{}/{}/{}'
PI = 3.1415926
km2deg = 9e-3

# CRIS L2 NH3 TPW Adjustment Values
SCALE = 0.0012
OFFSET = 0.0024

class l3_file: # All files going into a single end file
	def __init__(self, date, instrt, version, subversion):
		
		self.grid_res = 0.1
		self.map_bounds = [-89.5, 89.5, -179.5, 179.5]
		
		self.version = version
		self.subversion = subversion
		
		self.dlc=[]
		self.create_grids()
		
		self.instrt = instrt.upper()
		instrt_full = instrt
		
		self.date = [date[0:4],date[4:6],format(int(date[6:len(date)-1]),'02d')]
		self.daynight = date[len(date)-1:len(date)] # Last character
		self.timeformat = 'Monthly'
		
		self.inpath = BASE_INPUT.format(instrt_full, self.version, self.date[0],self.date[1])
		self.outpath = BASE_OUTPUT.format(instrt_full, self.version, self.date[0], self.date[1])
		print(self.outpath)
		
		print('Processing files for {}/{}/{}'.format(self.date[2],self.date[1],self.date[0]))
		self.main()
		
	def create_grids(self):
		self.lonlist = []
		self.latlist = []
		
		nh3calcdata, nh3tpwvza_calcdata, tpwcalcdata, tpwbcalcdata = [], [], [], []
		nh3basedata, nh3tpwvza_basedata, tpwbasedata, tpwbbasedata = [], [], [], []
		
		print('Creating Grids')
		
		num_lats = int((self.map_bounds[1] - self.map_bounds[0])/self.grid_res)
		num_lons = int((self.map_bounds[3] - self.map_bounds[2])/self.grid_res)
		
		print(num_lats, num_lons, self.grid_res)
		
		self.latlist = (np.arange(0,num_lats)*self.grid_res) + self.map_bounds[0]
		self.lonlist = (np.arange(0,num_lons)*self.grid_res) + self.map_bounds[2]
		
		base = np.reshape( np.repeat(0, num_lats*num_lons), (num_lats, num_lons)).astype('float')
		calc = np.reshape( np.repeat(0, 3*num_lats*num_lons), (3, num_lats, num_lons)).astype('float')
		
		self.nh3calcdata = calc.copy()
		self.nh3tpwvza_calcdata = calc.copy()
		self.tpwcalcdata = calc.copy()
		self.tpwbcalcdata = calc.copy()
		
		self.nh3basedata = base.copy()
		self.nh3tpwvza_basedata = base.copy()
		self.tpwbasedata = base.copy()
		self.tpwbbasedata = base.copy()
			
		print('Grids Created')
		
	def oversample_pixel(self, pixel_lat, pixel_lon, pixel_snum, nh3, tpw, tpwb, nh3_unc, tpw_unc, tpwb_unc, pixel_theta):
		# Takes single pixel with dimensions
		# Takes resolution of constant grid and initial coordinates
		# Determines overlaps of single pixel
		# Adds to Sum [0] and Denom [1] of 2d lat/lon array
		semi_major = 12 + (29/60)*abs(pixel_snum-60)
		semi_minor = 12 + (8/60)*abs(60-pixel_snum)
		
		ll_major = semi_major*km2deg
		ll_minor = semi_minor*km2deg
        
		circle = Point(pixel_lon, pixel_lat).buffer(1)
		ellipse = shapely.affinity.scale(circle, ll_major, ll_minor)
		ellipse = shapely.affinity.rotate(ellipse, pixel_theta)
		
		
		uncertainty_factor = PI*semi_major*semi_minor
        
		relative_lat = pixel_lat - self.map_bounds[0] - ll_major*1.5
		delta_lat    = self.map_bounds[1] - self.map_bounds[0]
		data_lat     = int(round( len(self.nh3calcdata[0]) * (relative_lat/delta_lat) ))
		
		relative_lon = pixel_lon - self.map_bounds[2]
		delta_lon    = self.map_bounds[3] - self.map_bounds[2]
		data_lon     = int(round( len(self.nh3calcdata[0][0]) * (relative_lon/delta_lon) ))
		
		is_complete = False
		square_count = 0
		while not is_complete:
			data_lat_check = int(round(data_lat + square_count))
			lat_check = self.map_bounds[0] + (data_lat_check/len(self.nh3calcdata[0]))*(self.map_bounds[1]-self.map_bounds[0])
			if data_lat_check < len(self.nh3calcdata[0]):
				for dl in [[-1,0],[1,1]]:
			
					is_zero = False
					side_count = dl[1]
					while not is_zero:
						data_lon_check = int(round(data_lon + side_count*dl[0]))
						
						lon_check = self.map_bounds[2] + (data_lon_check/len(self.nh3calcdata[0][0]))*(self.map_bounds[3]-self.map_bounds[2])
						if data_lon_check < len(self.nh3calcdata[0][0]):
							
							rectangle = Polygon([(lon_check    , lat_check    ),
												(lon_check+self.grid_res, lat_check    ),
												(lon_check+self.grid_res    , lat_check+self.grid_res),
												(lon_check, lat_check+self.grid_res) ])
							area_A = ellipse.intersection(rectangle).area
							if area_A > 0:
								
								if not np.isnan(nh3) and not np.isnan(nh3_unc):
									self.nh3calcdata[0][data_lat_check][data_lon_check] += area_A/(nh3_unc*uncertainty_factor)
									self.nh3calcdata[1][data_lat_check][data_lon_check] += (nh3*area_A)/(nh3_unc*uncertainty_factor)
									self.nh3calcdata[2][data_lat_check][data_lon_check] += 1
								
								if not np.isnan(tpw) and not np.isnan(tpw_unc):
									self.tpwcalcdata[0][data_lat_check][data_lon_check] += area_A/(tpw_unc*uncertainty_factor)
									self.tpwcalcdata[1][data_lat_check][data_lon_check] += (tpw*area_A)/(tpw_unc*uncertainty_factor)
									self.tpwcalcdata[2][data_lat_check][data_lon_check] += 1
								
								if not np.isnan(tpw) and not np.isnan(tpw_unc) and not np.isnan(nh3) and not np.isnan(nh3_unc):
									adjusted_unc = math.sqrt( (nh3_unc/nh3)**2 + (tpw_unc/tpw)**2 )
									nh3tpwvza = nh3 - (tpw*SCALE + OFFSET)
									nh3tpwvza_unc = nh3tpwvza*adjusted_unc
									
									self.nh3tpwvza_calcdata[0][data_lat_check][data_lon_check] += area_A/(nh3tpwvza_unc*uncertainty_factor)
									self.nh3tpwvza_calcdata[1][data_lat_check][data_lon_check] += (nh3tpwvza*area_A)/(nh3tpwvza_unc*uncertainty_factor)
									self.nh3tpwvza_calcdata[2][data_lat_check][data_lon_check] += 1
									
								if not np.isnan(tpwb) and not np.isnan(tpwb_unc):
									self.tpwbcalcdata[0][data_lat_check][data_lon_check] += area_A/(tpwb_unc*uncertainty_factor)
									self.tpwbcalcdata[1][data_lat_check][data_lon_check] += (tpwb*area_A)/(tpwb_unc*uncertainty_factor)
									self.tpwbcalcdata[2][data_lat_check][data_lon_check] += 1
								self.dlc.append([data_lat_check, data_lon_check])
								side_count += 1
							else:
								is_zero = True
								is_complete = (side_count == 0 and lat_check > pixel_lat)
								if is_complete:
									return None
						else:
							is_zero = True
			else:
				is_complete = True
						    
			square_count += 1
		return None
		
	def average_oversample(self, calcdata):
	    # Takes gridded data averaged
	    # Divide [0] by [1] for all lat/lon grid squares
	    # Returns basedata
	   data = np.take(calcdata,1, axis=0)
	   
	   area = np.take(calcdata,0, axis=0)
	   basedata = data/area
	   return basedata
		
	def add_frames_to_oversample(self,os_file):
	    # Takes single frame basedata
	    # Iterate through frame for each pixel
	    # Pass pixel to oversample pixel
	    # Receive gridded data
	    
		dt1000_f1 = True
		dt1000_f2 = False
		dt1000_f3 = False
	    
		os_ncf = Dataset(os_file,'r',format="NETCDF4")
		
		lat       = np.array(os_ncf['lat'])
		lon       = np.array(os_ncf['lon'])
		dflag     = np.array(os_ncf['day_flag'])
		
		nh3       = np.array(os_ncf['nh3'])
		nh3_err   = np.array(os_ncf['nh3_err'])
		tpw       = np.array(os_ncf['tpw'])
		tpw_err   = np.array(os_ncf['tpw'])
		tpw_base  = np.array(os_ncf['tpw_base'])
		tpw_base_err = tpw_base * (tpw_err/tpw)
		
		scl       = np.array(os_ncf['scanlinenum'])
		dt1km     = np.array(os_ncf['dt1000'])
		theta     = np.array(os_ncf['theta'])
		
		scl[scl is ma.masked] = 1
		
		
		dt1000_filter = [True for d in range(len(dt1km))]
		
		# Application of filters
		
		if dt1000_f1:
			dt1000_filter = dt1km < 15
		elif dt1000_f2:
			dt1000_filter = dt1km < 10
		elif dt1000_f3:
			dt1000_filter = dt1km < 20
		
		if self.daynight == 'd':
			zen_filter = dflag == 1
		elif self.daynight == 'n':
			zen_filter = dflag == 0
		else:
			zen_filter = [True for f in range(len(dflag))]
			
		full_filter = zen_filter #&dt1000_filter # Ignore for no filtering
		# All filters applied
		
		index_array = np.array([i for i in range(len(lat))])
		print('Total Number of pixels: {}'.format(len(index_array)))
		index_array = index_array[full_filter]
			
		print('{} pixels after filtering'.format(len(index_array)))
		idx = 0
		while idx < len(index_array):
			index = index_array[idx]
			nh3i = nh3[index]
			nh3ei = nh3_err[index]
			
			tpwi = tpw[index]
			tpwei = tpw_err[index]
			
			tpwbi = tpw_base[index]
			tpwbei = tpw_base_err[index]
			self.oversample_pixel(lat[index], lon[index], scl[index], nh3i, tpwi, tpwbi,
			                                                          nh3ei, tpwei, tpwbei,
			                                                          theta[index])
			idx += 1
		return None
		
	def save_to_nc(self):
		if self.daynight == 'd':
			daynight = 'day'
		elif self.daynight == 'n':
			daynight = 'night'
		else:
			daynight = 'daynight'
			
		if self.timeformat == 'Daily':
			base_file_format = '{}_oversampled_l2_{}{}{}_{}_{}_{}.nc'.format(self.instrt.lower(), self.date[0], self.date[1], self.date[2], daynight, self.version, self.subversion)
			out_dir = self.outpath + '/'
		else:
			base_file_format = '{}_oversampled_l2_{}{}_{}_{}_{}.nc'.format(self.instrt.lower(), self.date[0], self.date[1], daynight, self.version, self.subversion)
			
			out_dir = self.outpath + '/'
		if not(os.access(out_dir, os.F_OK)):
			os.makedirs(out_dir)
		if not os.path.isfile(out_dir+base_file_format):
			os.system('touch {}{}'.format(out_dir, base_file_format))
		print('Written to ',out_dir+base_file_format)
		ncfile = Dataset(out_dir+base_file_format, 'w', format = 'NETCDF4')
		
		lat_dim = ncfile.createDimension('lat',len(self.latlist))
		lon_dim = ncfile.createDimension('lon',len(self.lonlist))
		
		lat_var = ncfile.createVariable('lat', np.float32, ('lat',))
		lat_var.long_name = 'latitude'
		lat_var[:] = self.latlist
		
		lon_var = ncfile.createVariable('lon', np.float32, ('lon',))
		lon_var.long_name = 'longitude'
		lon_var[:] = self.lonlist
		
		nh3 = ncfile.createVariable('nh3', np.float32, ('lat','lon',))
		nh3.long_name = 'Oversampled Ammonia'
		nh3.units = 'ppbv'
		nh3[:,:] = self.nh3basedata
		
		nh3tpwvza = ncfile.createVariable('nh3tpwvza', np.float32, ('lat','lon',))
		nh3tpwvza.long_name = 'Oversampled Adjusted Ammonia'
		nh3tpwvza.units = 'ppbv'
		nh3tpwvza[:,:] = self.nh3tpwvza_basedata
		
		
		tpw = ncfile.createVariable('tpw', np.float32, ('lat','lon',))
		tpw.long_name = 'Oversampled Total Precipitable Water with angle adjustment'
		tpw.units = 'mm'
		tpw[:,:] = self.tpwbasedata
		
		tpw_base = ncfile.createVariable('tpw_base', np.float32, ('lat','lon',))
		tpw_base.long_name = 'Oversampled TPW Base'
		tpw_base.units = 'mm'
		tpw_base[:,:] = self.tpwbbasedata
		
		ncfile.close()
		 
	def main(self):
		print('Main')
		if self.timeformat == 'Daily':
			# Do single file
			filein = self.inpath + '/{}_valid_l2a_{}{}{}_{}.nc'.format(self.instrt.lower(), self.date[0], self.date[1], self.date[2], self.version)
			self.add_frames_to_oversample(filein)
			self.nh3basedata = self.average_oversample(self.nh3calcdata)
			self.nh3tpwvza_basedata = self.average_oversample(self.nh3tpwvza_calcdata)
			self.tpwbasedata = self.average_oversample(self.tpwcalcdata)
			self.tpwbbasedata = self.average_oversample(self.tpwbcalcdata)
			self.save_to_nc()
		elif self.timeformat == 'Monthly':
			for day in range(0,30):
				filein = self.inpath + '/{}_valid_l2a_{}{}{}_{}.nc'.format(self.instrt.lower(), self.date[0], self.date[1], format(int(self.date[2])+day, '02d'), self.version)
				self.add_frames_to_oversample(filein)

			self.nh3basedata = self.average_oversample(self.nh3calcdata)
			self.nh3tpwvza_basedata = self.average_oversample(self.nh3tpwvza_calcdata)
			self.tpwbasedata = self.average_oversample(self.tpwcalcdata)
			self.tpwbbasedata = self.average_oversample(self.tpwbcalcdata)
			self.save_to_nc()
				
def main(date):

	l3f = l3_file(date)	

options, operands = getopt(sys.argv[1:], "", ["date=","instrt=","version=","subversion="])
# Determine options and processing method

date = operands[0]
version = operands[2]
subversion = operands[3]
instrt = operands[1]
l3f = l3_file(date, instrt, version, subversion)
