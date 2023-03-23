"""
L2 oversampling tool
--------------
Takes L2 pixel-by-pixel filtered files and performs oversampling to create daily files.

"""
__author__	= "Daniel Westwood"
__date__	  = "23-03-2023"
__copyright__ = "Copyright 2020 United Kingdom Research and Innovation"
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

from matplotlib.patches import Ellipse, Rectangle

from getopt import getopt

import math

# Input from l2 base-filtered files, output to oversampled l2 files
BASE_INPUT  = '/gws/pw/j05/rsg_share/public/transfer/barry/6ad76b04-b3df-11eb-a4d3-024ad8e814ad/{}_l2/valid_l2a/{}/{}/{}'
BASE_OUTPUT = '/gws/pw/j05/rsg_share/public/transfer/barry/6ad76b04-b3df-11eb-a4d3-024ad8e814ad/{}_l2/oversampled_l2/{}/{}/{}'
PI = 3.1415926
km2deg = 9e-3

def remap_to_new_grid(lat, lon, bounds, data):
	# Convenient function for mapping to differently sized grids
	# Used for SD filtering with reference grid of a different size to the os grid
	abs_lat = lat - bounds[0]
	abs_lon = lon - bounds[2]
	ratio_lat = (data[1]-data[0])/(bounds[1]-bounds[0])
	ratio_lon = (data[3]-data[2])/(bounds[3]-bounds[2])
	
	data_lat = int(abs_lat * ratio_lat)
	data_lon = int(abs_lon * ratio_lon)
	return data_lat, data_lon

class var_save:
	# Class for simplifying saving data to netcdf files (partial oversampling)
	def __init__(self, ncfile, vshort, vlong, vunits, vdata):
		self.den = ncfile.createVariable(vshort+'_den', np.float32, ('lat','lon',))
		self.den.long_name = vlong + ' Denominator'
		self.den.units = vunits
		self.den[:,:] = vdata[0]
		
		self.den = ncfile.createVariable(vshort+'_num', np.float32, ('lat','lon',))
		self.den.long_name = vlong + ' Numerator'
		self.den.units = vunits
		self.den[:,:] = vdata[1]
		
		self.ncfile = ncfile
		
class single_var:
	# Class for simplifying saving data to netcdf files (final oversampling)
	def __init__(self,ncfile,vshort,vlong, vunits, vdata):
		self.var = ncfile.createVariable(vshort,np.float32,('lat','lon'))
		self.var.long_name = vlong
		self.var.units = vunits
		self.var[:,:] = vdata
		
		self.ncfile = ncfile

class l3_file: # All files going into a single end file
	def __init__(self, date, instrt, version, subversion, time, variables, longs, units):
		
		self.grid_res = 0.1 # Higher resolution oversampled grid
		self.map_bounds = [-89.5, 89.5, -179.5, 179.5]
		self.timeformat = time
			
		self.variables  = variables
		self.longs      = longs
		self.units      = units	
		
		self.version    = version
		self.subversion = subversion
		
		self.create_grids()
		
		self.instrt     = instrt.upper()
		instrt_full     = instrt
		
		self.date = [date[0:4],date[4:6],format(int(date[6:len(date)-1]),'02d')]
		self.daynight = date[len(date)-1:len(date)] # Last character
		if self.daynight == 'd':
			daynight = 'day'
		elif self.daynight == 'n':
			daynight = 'night'
		else:
			daynight = 'daynight'
		
		#Extract Known Means for SD filtering
		try:
			ref_filepath = '/gws/pw/j05/rsg_share/public/transfer/barry/6ad76b04-b3df-11eb-a4d3-024ad8e814ad/{}_l2/'.format(instrt)
			ref_filename = 'oversampled_l2/{}/{}/{}/{}_os_reference_l2_{}{}{}_{}_ref.nc'.format(version, self.date[0], self.date[1], instrt, self.date[0], self.date[1], daynight, version, subversion)
			print(ref_filepath+ref_filename)
			ncf = Dataset(ref_filepath+ref_filename, 'r', format='NETCDF4')
			self.nh3_ref = np.array(ncf['nh3'])
			self.nh3_ref[np.isnan(self.nh3_ref)]=0
			self.is_sd = True
			
			lats = np.array(ncf['lat'])
			lons = np.array(ncf['lon'])
			self.ref_bounds = [lats[0], lats[len(lats)-1], lons[0], lons[len(lons)-1]]
			self.ref_data = [0, len(lats), 0, len(lons)]
			ncf.close()
			
		except:
			self.is_sd = False
			print('SDError: Standard Dev. Filter failed - Reference file not found')
		# Force no SD checking for now
		self.is_sd = False
		self.timeformat = time
		
		self.inpath = BASE_INPUT.format(instrt_full, self.version, self.date[0],self.date[1])
		self.outpath = BASE_OUTPUT.format(instrt_full, self.version, self.date[0], self.date[1])
		
		self.main()
		
	def create_grids(self):
		self.lonlist = []
		self.latlist = []
		
		print('Creating Grids')
		# Create static grids for the oversampling routines
		num_lats = int((self.map_bounds[1] - self.map_bounds[0])/self.grid_res)
		num_lons = int((self.map_bounds[3] - self.map_bounds[2])/self.grid_res)
		lvs = len(self.variables)
		self.basedatae = np.reshape( np.repeat(0, num_lats*num_lons*lvs), (lvs, num_lats, num_lons)).astype('float')
		# Save three sections of partial oversampled data
		self.calcdatae = np.reshape( np.repeat(0, 3*num_lats*num_lons*lvs), (lvs, 3, num_lats, num_lons)).astype('float')
		self.grid_squares = np.reshape(np.repeat(np.nan, num_lats*num_lons), (num_lats, num_lons)).astype('float')
		
		
		self.latlist = (np.arange(0,num_lats)*self.grid_res) + self.map_bounds[0]
		self.lonlist = (np.arange(0,num_lons)*self.grid_res) + self.map_bounds[2]
			
		print('Grids Created')
		
	def oversample_pixel(self, pixel_lat, pixel_lon, pixel_snum, values, errors, pixel_theta):
		# Takes single pixel with dimensions
		# Takes resolution of constant grid and initial coordinates
		# Determines overlaps of single pixel
		# Adds to Sum [0] and Denom [1] of 2d lat/lon array
		
		# Pixel sizes for iasi and cris
		if self.instrt.lower() == 'iasi':
			semi_major = 12 + (29/60)*abs(pixel_snum-60)
			semi_minor = 12 + (8/60)*abs(60-pixel_snum)
		elif self.instrt.lower() == 'cris':
			semi_major = 14 + (29/60)*abs(pixel_snum-60)
			semi_minor = 14 + (8/60)*abs(60-pixel_snum)
		
		# Future pixel size calculations for different instruments
		# Use pixel_snum to represent across-track number (0-120 with nadir at 60 for IASI/CRIS)
		# For multiple modes, pixel_snum can be interpreted as an array

		# Radius calculation
		ll_major = (semi_major*km2deg)/2
		ll_minor = (semi_minor*km2deg)/2 # Correction
        
        # Apply sizes to shape generation
		circle = Point(pixel_lon, pixel_lat).buffer(1)
		ellipse = shapely.affinity.scale(circle, xfact=ll_major, yfact=ll_minor)
		ellipse = shapely.affinity.rotate(ellipse, pixel_theta)
		
		# Constant used in all calculations for this pixel (area fraction requires full pixel area)
		uncertainty_factor = PI*semi_major*semi_minor
        
        # Calculate lat/lon for pixels and locations on the grid
		relative_lat = pixel_lat - self.map_bounds[0] - ll_major*1.5
		delta_lat    = self.map_bounds[1] - self.map_bounds[0]
		data_lat     = int(round( len(self.calcdatae[0][0]) * (relative_lat/delta_lat) ))
		
		relative_lon = pixel_lon - self.map_bounds[2]
		delta_lon    = self.map_bounds[3] - self.map_bounds[2]
		data_lon     = int(round( len(self.calcdatae[0][0][0]) * (relative_lon/delta_lon) ))
		
		direction_matrix = [[-1,0],[1,1]]
		# 0th indicates right/left longitude displacement
		# 1th indicates starting position
		# See l2_l2 documentation.pptx
		
		
		# Count squares either side of starting point until zero-area overlap reached
		is_complete = False
		row_count = 0
		# row count above/below pixel (lat)
		
		while not is_complete: # All squares considered for pixel
			# lat_check dependent on current square being looked at
			
			data_lat_check = int(round(data_lat + row_count))
			lat_check = self.map_bounds[0] + (data_lat_check/len(self.calcdatae[0][0]))*(self.map_bounds[1]-self.map_bounds[0])
			
			if data_lat_check < len(self.calcdatae[0][0]):
				for component in direction_matrix:
					## Count to either side of the square
					
					is_zero = False
					column_count = component[1]
					
					# delta lon (direction_matrix) parts [-1, 0] and [1,1]
					# where 0th part is starting point for column count (0-left, 1-right so 0 is not double counted)
					# and 1th part indicates (negative-left, positive-right)
					while not is_zero:
						# Test new square overlaps to left/right until overlap is zero
						
						data_lon_check = int(round(data_lon + column_count*component[0]))
						
						lon_check = self.map_bounds[2] + (data_lon_check/len(self.calcdatae[0][0][0]))*(self.map_bounds[3]-self.map_bounds[2])
						if data_lon_check < len(self.calcdatae[0][0][0]):
											
							if True: # For future filters
								
								# Assemble rectangle shape for intersection testing (must use lat/lon actual values, rather than data values)
								rectangle = Polygon([(lon_check    , lat_check    ),
													(lon_check+self.grid_res, lat_check    ),
													(lon_check+self.grid_res    , lat_check+self.grid_res),
													(lon_check, lat_check+self.grid_res) ])
								area_A = ellipse.intersection(rectangle).area
								if area_A > 0:
									progress = True
									
									# Future SD Filtering at Oversampling stage
									
#									nh3_mean = 0
									
									# Standard Deviation Filtering - Reference Grid (disabled)
#									if self.is_sd:
										# SD Filter
#										ref_lat, ref_lon = remap_to_new_grid(lat_check, lon_check, self.ref_bounds, self.ref_data)
#										nh3_mean = self.nh3_ref[ref_lat, ref_lon]
#										if values[0] < nh3_mean + 3*SD and values[0] > nh3_mean - 3*SD:
#											progress = True
#										else:
#											progress = False
#									else:
#										progress = True


									if progress:
										for iv, value in enumerate(values):
											error = errors[iv]
											if not np.isnan(value) and not np.isnan(error):
												# Save all values to grid
												if iv == 1:
													# Adjustment
													value = values[0] - (values[2]*SCALE + OFFSET)
													error = errors[0]
												# Sum numerator/denominator parts to the oversampling equation
												self.calcdatae[iv][0][data_lat_check][data_lon_check] += area_A/(error*uncertainty_factor)
												self.calcdatae[iv][1][data_lat_check][data_lon_check] += (value*area_A)/(error*uncertainty_factor)
												self.calcdatae[iv][2][data_lat_check][data_lon_check] += 1

									# Add 1 to the side count for next square
									column_count += 1
								else:
									# Zero overlap - edge reached
									is_zero = True
									# If edge reached at 0 left/right displacement, pixel is finished
									is_complete = (column_count == 0 and lat_check > pixel_lat)
									if is_complete:
										return None
						else:
							is_zero = True
			else:
				is_complete = True
						    
			row_count += 1
			# Add 1 to the row count for completed row
		return None
		
	def average_oversample(self, calcdata):
	    # Takes gridded data averaged
	    # Divide [0] by [1] for all lat/lon grid squares
	    # Returns basedata
	    
	   data = np.take(calcdata,1, axis=0)
	   
	   area = np.take(calcdata,0, axis=0)
	   basedata = data/area
	   basedata[np.isnan(basedata)]=0
	   
	   return basedata
		
	def add_frames_to_oversample(self,os_file):
	    # Takes single frame basedata
	    # Iterate through frame for each pixel
	    # Pass pixel to oversample pixel
	    # Receive gridded data
	    
		os_ncf = Dataset(os_file,'r',format="NETCDF4")
		
		# Extract necessary variables from base-filtered files
		
		lat       = np.array(os_ncf['lat'])
		lon       = np.array(os_ncf['lon'])
		dflag     = np.array(os_ncf['day_flag'])
		
		nh3       = np.array(os_ncf['nh3'])
		nh3_err   = np.array(os_ncf['nh3_err'])

		tpw       = np.array(os_ncf['tpw'])
		tpw_err   = np.array(os_ncf['tpw'])
		tpw_base  = np.array(os_ncf['tpw_base'])
		tpw_base_err = tpw_base * abs((tpw_err/tpw))
		
		dt1km     = np.array(os_ncf['dt1000'])
		error_corr = True # scanlinenum and theta are swapped in some files
		if error_corr:
			scl       = np.array(os_ncf['theta'])
			theta     = np.array(os_ncf['scanlinenum'])
		else:
			scl       = np.array(os_ncf['scanlinenum'])
			theta     = np.array(os_ncf['theta'])
		
		scl[scl is ma.masked] = 1
		
		# Testing dt1000 filter is performing as expected!
		
		#print('Comparison of values')
		#print('Initial:')
		#print('dt1000 mean: {}, tpwvza mean: {}, nh3 mean: {}'.format(np.nanmean(dt1km), np.nanmean(tpw), np.nanmean(nh3)))
		#dt1000_f = dt1km > 15
		#print('Filtered:')
		#print('dt1000 mean: {}, tpwvza mean: {}, nh3 mean: {}'.format(np.nanmean(dt1km[dt1000_f]), np.nanmean(tpw[dt1000_f]), np.nanmean(nh3[dt1000_f])))
		#x=input()
		
		## Application of filters
		
		# Mandatory Filters
		if self.daynight == 'd': # day
			zen_filter = dflag == 1
		elif self.daynight == 'n': # night
			zen_filter = dflag == 0
		else: # Day and night
			zen_filter = [True for f in range(len(dflag))]
		tpw_filter_h = tpw < 40
		tpw_filter_l = tpw > 5
		
		man_filters = zen_filter & tpw_filter_h & tpw_filter_l
		
		# Optional Filters
		scl_factor = 0 # Scan angle filtering
		dt1000_factor = 0 # Dt1000 filtering
		
		full_filter = man_filters	
		if scl_factor != 0:
			scl_filter = np.absolute(scl-60) < scl_factor
			full_filter = full_filter & scl_filter
		if dt1000_factor != 0:
			dt1000_filter = dt1km > dt1000_factor
			full_filter = full_filter & dt1000_filter
			
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
			self.oversample_pixel(lat[index], lon[index], scl[index], [nh3i, 0, tpwi, tpwbi],
			                                                          [nh3ei, 0, tpwei, tpwbei],
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
			base_file_format = '{}_oversampled_partial_l2_{}{}{}_{}_{}_{}.nc'.format(self.instrt.lower(), self.date[0], self.date[1], self.date[2], daynight, self.version, self.subversion)
			out_dir = self.outpath + '/'
		else:
			base_file_format = '{}_oversampled_l2_{}{}_{}_{}_{}.nc'.format(self.instrt.lower(), self.date[0], self.date[1], daynight, self.version, self.subversion)
			
			out_dir = self.outpath + '/'
		if not(os.access(out_dir, os.F_OK)):
			os.makedirs(out_dir)
		if not os.path.isfile(out_dir+base_file_format):
			os.system('touch {}{}'.format(out_dir, base_file_format))
		print('Written to ',out_dir+base_file_format)
		print(np.nanmax(self.calcdatae[0]), np.nanmin(self.calcdatae[0]))
		ncfile = Dataset(out_dir+base_file_format, 'w', format = 'NETCDF4')
		
		lat_dim = ncfile.createDimension('lat',len(self.latlist))
		lon_dim = ncfile.createDimension('lon',len(self.lonlist))
		
		lat_var = ncfile.createVariable('lat', np.float32, ('lat',))
		lat_var.long_name = 'latitude'
		lat_var[:] = self.latlist
		
		lon_var = ncfile.createVariable('lon', np.float32, ('lon',))
		lon_var.long_name = 'longitude'
		lon_var[:] = self.lonlist
		
		## Class system
		var_objs_arr = [None for var in range(len(self.variables))]
		if self.timeformat == 'Daily':
			for iv in range(len(self.variables)):
				var_objs_arr[iv] = var_save(ncfile, self.variables[iv],self.longs[iv],self.units[iv],self.calcdatae[iv])
				ncfile = var_objs_arr[iv].ncfile
			
		else:
			for iv in range(len(self.variables)):
				var_objs_arr[iv] = single_var(ncfile,self.variables[iv],self.longs[iv],self.units[iv],self.basedatae[iv])
				ncfile = var_objs_arr[iv].ncfile
		
		ncfile.close()
		 
	def add_frames_to_average(self,filein):
		
		ncfile = Dataset(filein, 'r', format = 'NETCDF4')
		for iv in range(len(self.variables)):
			denlabel = self.variables[iv] + '_den'
			numlabel = self.variables[iv] + '_num'
			self.calcdatae[iv]   = self.addto_calcdata(self.calcdatae[iv], np.array(ncfile[denlabel]), np.array(ncfile[numlabel]))
		
	def addto_calcdata(self,calcdata, newdata0, newdata1,label=False):
		
		calcdata[0] = calcdata[0] + newdata0
		calcdata[1] = calcdata[1] + newdata1
		return calcdata
		 
	def main(self):
		print('Main')
		print(self.timeformat)
		if self.timeformat == 'Daily':
			# Do single file
			filein = self.inpath + '/{}_valid_l2a_{}{}{}_{}.nc'.format(self.instrt.lower(), self.date[0], self.date[1], self.date[2], self.version)
			self.add_frames_to_oversample(filein)
			self.save_to_nc()
		elif self.timeformat == 'Monthly':
			if self.daynight == 'd':
				daynight = 'day'
			elif self.daynight == 'n':
				daynight = 'night'
			else:
				daynight = 'daynight'
			
			
			for day in range(0,31):
				print(day)
				base_file_format = '{}_oversampled_partial_l2_{}{}{}_{}_{}_{}.nc'.format(self.instrt.lower(), self.date[0], self.date[1], format(day+1,'02d'), daynight, self.version, self.subversion)
			
				out_dir = self.outpath + '/'
				try:
					self.add_frames_to_average(out_dir+base_file_format)
				except:
					pass
			for iv in range(len(self.variables)):
				self.basedatae[iv] = self.average_oversample(self.calcdatae[iv])
			self.save_to_nc()
		elif self.timeformat == 'STM':
			if self.daynight == 'd':
				daynight = 'day'
			elif self.daynight == 'n':
				daynight = 'night'
			else:
				daynight = 'daynight'
			for day in range(1,32):
				print('DAY: {}'.format(day))
				filein = self.inpath + '/{}_valid_l2a_{}{}{}_{}.nc'.format(self.instrt.lower(), self.date[0], self.date[1], format(day,'02d'), self.version)
				try:
					self.add_frames_to_oversample(filein)
				except FileNotFoundError:
					pass
			# Do SD filtering
			for iv in range(len(self.variables)):
				self.basedatae[iv] = self.average_oversample(self.calcdatae[iv])
			self.save_to_nc()

if __name__ == '__main__':
    variables = ['nh3','nh3tpwvza','tpw','tpwb']
    longs = ['Oversampled Ammonia',
    				'Oversampled Adjusted Ammonia',
    				'Oversampled TPWVZA',
    				'Oversampled TPW']
    units = ['ppbv','ppbv','mm','mm']
    
    options, operands = getopt(sys.argv[1:], "", ["date=","instrt=","version=","subversion=","time="])
    # Determine options and processing method
    
    date = operands[0]
    version = operands[2]
    subversion = operands[3]
    instrt = operands[1]
    time = operands[4]
    
    # CRIS L2 NH3 TPW Adjustment Values
    if instrt.lower() == 'cris':
    	SCALE = 0.0012
    	OFFSET = 0.0024
    
    # IASI L2 NH3 TPW Adjustment Values
    elif instrt.lower() == 'iasi':
    	SCALE = 0.003
    	OFFSET = -0.06
    SD = 0.09
    
    l3f = l3_file(date, instrt, version, subversion, time, variables, longs, units)
