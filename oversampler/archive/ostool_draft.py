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
import file_import as fm
import find_files as ff

BASE_INPUT = '/gws/pw/j05/rsg_share/public/projects/ims/data/lv2/output_ghg_cv9_v12_lam_fgsnwp_nobc_newbc_rbc_ram2_rem20-ncam_rcl5_raer6st_rmg_oap_cnv2_o3f113_csyn_swir2/npp/{}/{}/{}'
PI = 3.1415926
km2deg = 9e-3

class l3_file: # All files going into a single end file
	def __init__(self, date):
		
		self.grid_res = 0.1
		self.map_bounds = [40, 60, -20 ,5]
		
		self.create_grids()
		
		self.variable = 'mgf'
		self.vshort = 'nh3'
		self.instrt = 'IASI'
		
		# Date = '20200402'
		
		self.date = [date[0:4],date[4:6],format(int(date[6:len(date)-1]),'02d')]
		self.daynight = date[len(date)-1:len(date)] # Last character
		print(self.date)
		
		self.timeformat = 'Daily'
		self.inpath = BASE_INPUT.format(self.date[0],self.date[1],self.date[2])
		self.outpath = '/home/users/dwest77/Documents/IASI_os_l3/'
		print('Processing files for {}/{}/{}'.format(date[2],date[1],date[0]))
		self.main()
		
	def create_grids(self):
		self.calcdata = []
		self.basedata = []
		self.lonlist = []
		self.latlist = []
		
		print('Creating Grids')
		lat_counter = self.map_bounds[0]
		while lat_counter < self.map_bounds[1]:
			lon_counter = self.map_bounds[2]
			data_row = []
			base_row = []
			while lon_counter < self.map_bounds[3]:
				if lat_counter == self.map_bounds[0]:
					self.lonlist.append(lon_counter)
				data_row.append([0,0,0]) # Denominator, Numerator, Count
				base_row.append(0)
				
				lon_counter += self.grid_res
			self.latlist.append(lat_counter)
			lat_counter += self.grid_res
			self.calcdata.append(data_row)
			self.basedata.append(base_row)
			
		print(len(self.calcdata), len(self.calcdata[0]))
			
		print('Grids Created')
		
	def oversample_pixel(self, pixel_lat, pixel_lon, pixel_snum, pixel_val, pixel_unc, pixel_theta):
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
		
		area_unc = pixel_unc*PI*semi_major*semi_minor #pi*a*b
        
		relative_lat = pixel_lat - self.map_bounds[0] - ll_major*1.5
		delta_lat    = self.map_bounds[1] - self.map_bounds[0]
		data_lat     = int(round( len(self.calcdata) * (relative_lat/delta_lat) ))
		
		relative_lon = pixel_lon - self.map_bounds[2]
		delta_lon    = self.map_bounds[3] - self.map_bounds[2]
		data_lon     = int(round( len(self.calcdata[0]) * (relative_lon/delta_lon) ))
		
		is_complete = False
		square_count = 0
		while not is_complete:
			data_lat_check = int(round(data_lat + square_count))
			lat_check = self.map_bounds[0] + (data_lat_check/len(self.calcdata))*(self.map_bounds[1]-self.map_bounds[0])
			if data_lat_check < len(self.calcdata):
				for dl in [[-1,0],[1,1]]:
			
					is_zero = False
					side_count = dl[1]
					while not is_zero:
						data_lon_check = int(round(data_lon + side_count*dl[0]))
						lon_check = self.map_bounds[2] + (data_lon_check/len(self.calcdata[0]))*(self.map_bounds[3]-self.map_bounds[2])
						if data_lon_check < len(self.calcdata[0]):
							rectangle = Polygon([(lon_check    , lat_check    ),
												(lon_check+self.grid_res, lat_check    ),
												(lon_check+self.grid_res    , lat_check+self.grid_res),
												(lon_check, lat_check+self.grid_res) ])
							area_A = ellipse.intersection(rectangle).area
							if area_A > 0:
				
								common = area_A/area_unc
								self.calcdata[data_lat_check][data_lon_check][0] += common
								self.calcdata[data_lat_check][data_lon_check][1] += common*pixel_val
								self.calcdata[data_lat_check][data_lon_check][2] += 1
				
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
		
	def average_oversample(self):
	    # Takes gridded data averaged
	    # Divide [0] by [1] for all lat/lon grid squares
	    # Returns basedata
		self.max = 0
		self.min = 999
	    
	    
		print('Averaging')
		for irow in range(len(self.calcdata)):
			for icol in range(len(self.calcdata[0])):
				if self.calcdata[irow][icol][0] == 0:
					self.basedata[irow][icol] = 0
				else:
					self.basedata[irow][icol] = (self.calcdata[irow][icol][1]/self.calcdata[irow][icol][0])
					if self.basedata[irow][icol] > self.max:
						self.max = self.basedata[irow][icol]
					if self.basedata[irow][icol] < self.min:
						self.min = self.basedata[irow][icol]
												
		return None
		
	def add_frames_to_oversample(self,os_file):
	    # Takes single frame basedata
	    # Iterate through frame for each pixel
	    # Pass pixel to oversample pixel
	    # Receive gridded data
		os_ncf = Dataset(os_file,'r',format="NETCDF4")
		# Determine approximate angle
		valid_angle = False
		index = 0
		while not valid_angle:
			npixel_lat = os_ncf['latitude'][index]
			npixel_lon = os_ncf['longitude'][index]
		
			pixel_lat = os_ncf['latitude'][index+1]
			pixel_lon = os_ncf['longitude'][index+1]
			frac = (npixel_lat - pixel_lat) / (npixel_lon - pixel_lon)
			pixel_theta = math.atan(frac)*(180/PI)
			if not pixel_theta is ma.masked and not np.isnan(pixel_theta):
				valid_angle = True
			index += 1
		is_day = False
		pixels = []
		mcount = 0
		for index in range(0, len(os_ncf['latitude'])):
			if True:
				lat = os_ncf['latitude'][index]
				lon = os_ncf['longitude'][index]
				zen = os_ncf['solzen'][index]
				if lon < self.map_bounds[2] or lon > self.map_bounds[3]:
					pass
				elif lat < self.map_bounds[0] or lat > self.map_bounds[1]:
					pass
				else:
					val = os_ncf[self.variable][index][1]
					val_err = os_ncf[self.variable+'_err'][index][1]
					if not val is ma.masked and not val_err is ma.masked and not np.isnan(val) and not np.isnan(val_err):
						pixels.append(index)
				if zen < 90:
					is_day = True
			if not os_ncf[self.variable][index][1] is ma.masked:
				mcount += 1
		print('{} pixels within range'.format(len(pixels)))
		if (is_day and self.daynight == 'd') or (not is_day and self.daynight == 'n'):
			for pdex, index in enumerate(pixels):
				if not os_ncf[self.variable][index][1] is ma.masked:
					scanlinenum = os_ncf['ixt'][index]
					if scanlinenum is ma.masked:
						scanlinenum = 1
					self.oversample_pixel(os_ncf['latitude'][index],
									os_ncf['longitude'][index],
									scanlinenum,
									(os_ncf[self.variable][index][1]),
									abs(os_ncf[self.variable+'_err'][index][1]),
									pixel_theta)

			if len(pixels) == 0:
				print('File Contains No relevant data')
		else:
			print('File contains no data matching {} requirements'.format(self.daynight))
		return None
		
	def save_to_nc(self):
		if self.timeformat == 'Daily':
			if self.daynight == 'd':
				daynight = 'day'
			else:
				daynight = 'night'
			base_file_format = '{}_{}_{}_cv9_v12_{}_{}_{}_{}.nc'.format(self.instrt, self.vshort, self.timeformat, self.date[0], self.date[1], self.date[2], daynight)
			out_dir = self.outpath + '{}/{}/{}/'.format(self.vshort.upper(), self.date[0], self.date[1])
		else:
			base_file_format = '{}_{}_{}_cv9_v12_{}_{}.nc'.format(self.instrt, self.vshort, self.timeformat, self.date[0], self.date[1])
			out_dir = self.outpath + '{}/{}/'.format(self.vshort.upper(), self.date[0])
		if not(os.access(out_dir, os.F_OK)):
			os.makedirs(out_dir)
		
		ncfile = Dataset(out_dir+base_file_format, 'w', format = 'NETCDF4')
		
		lat_dim = ncfile.createDimension('lat',len(self.latlist))
		lon_dim = ncfile.createDimension('lon',len(self.lonlist))
		
		lat_var = ncfile.createVariable('lat', np.float32, ('lat',))
		lat_var.long_name = 'latitude'
		lat_var[:] = self.latlist
		
		lon_var = ncfile.createVariable('lon', np.float32, ('lon',))
		lon_var.long_name = 'longitude'
		lon_var[:] = self.lonlist
		
		data_var = ncfile.createVariable(self.vshort, np.float32, ('lat','lon',))
		data_var.long_name = self.variable
		data_var[:,:] = self.basedata
		
		ncfile.close()
		 
	def main(self):
		print('Main')
		if self.timeformat == 'Daily':
			print(self.inpath)
			files_array = ff.list_files(self.inpath, starts='ral',ends='.nc')
			print('found {} files'.format(len(files_array)))
			if len(files_array) > 0:
			
				for ix in range(0,len(files_array)):
					print('File {} of {} '.format(ix+1, len(files_array)))
					self.add_frames_to_oversample(files_array[ix])
				self.average_oversample()
				self.save_to_nc()
			else:
				print('{}/{}/{} skipped'.format(self.date[2],self.date[1],self.date[0]))
				
def main(date):

	l3f = l3_file(date)	

options, operands = getopt(sys.argv[1:], "", ["date="])
# Determine options and processing method

date = operands[0]
main(date)
