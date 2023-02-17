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

def remap_to_new_grid(lat, lon, bounds, data):
	abs_lat = lat - bounds[0]
	abs_lon = lon - bounds[2]
	ratio_lat = (data[1]-data[0])/(bounds[1]-bounds[0])
	ratio_lon = (data[3]-data[2])/(bounds[3]-bounds[2])
	
	data_lat = int(abs_lat * ratio_lat)
	data_lon = int(abs_lon * ratio_lon)
	return data_lat, data_lon

class var_save:
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
	def __init__(self,ncfile,vshort,vlong, vunits, vdata):
		self.var = ncfile.createVariable(vshort,np.float32,('lat','lon'))
		self.var.long_name = vlong
		self.var.units = vunits
		self.var[:,:] = vdata
		
		self.ncfile = ncfile

class l3_file: # All files going into a single end file
	def __init__(self, date, instrt, version, subversion, time, variables, longs, units):
		
		self.grid_res = 0.1
		self.map_bounds = [-89.5, 89.5, -179.5, 179.5]
		self.timeformat = time
		if time == 'Daily':
			self.arr_len = 5
		else:
			self.arr_len = 50
			
		self.variables = variables
		self.longs = longs
		self.units = units	
		
		self.version = version
		self.subversion = subversion
		
		self.create_grids()
		
		self.instrt = instrt.upper()
		instrt_full = instrt
		
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
		print(self.outpath)
		
		print('Processing files for {}/{}/{}'.format(self.date[2],self.date[1],self.date[0]))
		self.main()
		
	def create_grids(self):
		self.lonlist = []
		self.latlist = []
		
		print('Creating Grids')
		
		num_lats = int((self.map_bounds[1] - self.map_bounds[0])/self.grid_res)
		num_lons = int((self.map_bounds[3] - self.map_bounds[2])/self.grid_res)
		lvs = len(self.variables)
		self.basedatae = np.reshape( np.repeat(0, num_lats*num_lons*lvs), (lvs, num_lats, num_lons)).astype('float')
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
		if self.instrt.lower() == 'iasi':
			semi_major = 12 + (29/60)*abs(pixel_snum-60)
			semi_minor = 12 + (8/60)*abs(60-pixel_snum)
		elif self.instrt.lower() == 'cris':
			semi_major = 14 + (34/135)*abs(pixel_snum-135)
			semi_minor = 14 + (9.3/135)*abs(135-pixel_snum)
		
		ll_major = (semi_major*km2deg)/2
		ll_minor = (semi_minor*km2deg)/2 # Correction
        
		circle = Point(pixel_lon, pixel_lat).buffer(1)
		ellipse = shapely.affinity.scale(circle, xfact=ll_major, yfact=ll_minor)
		ellipse = shapely.affinity.rotate(ellipse, pixel_theta)
		
		
		uncertainty_factor = PI*semi_major*semi_minor
        
		relative_lat = pixel_lat - self.map_bounds[0] - ll_major*1.5
		delta_lat    = self.map_bounds[1] - self.map_bounds[0]
		data_lat     = int(round( len(self.calcdatae[0][0]) * (relative_lat/delta_lat) ))
		
		relative_lon = pixel_lon - self.map_bounds[2]
		delta_lon    = self.map_bounds[3] - self.map_bounds[2]
		data_lon     = int(round( len(self.calcdatae[0][0][0]) * (relative_lon/delta_lon) ))
		
		is_complete = False
		square_count = 0
		while not is_complete:
			data_lat_check = int(round(data_lat + square_count))
			lat_check = self.map_bounds[0] + (data_lat_check/len(self.calcdatae[0][0]))*(self.map_bounds[1]-self.map_bounds[0])
			if data_lat_check < len(self.calcdatae[0][0]):
				for dl in [[-1,0],[1,1]]:
			
					is_zero = False
					side_count = dl[1]
					while not is_zero:
						data_lon_check = int(round(data_lon + side_count*dl[0]))
						
						lon_check = self.map_bounds[2] + (data_lon_check/len(self.calcdatae[0][0][0]))*(self.map_bounds[3]-self.map_bounds[2])
						if data_lon_check < len(self.calcdatae[0][0][0]):
											
							if True:
							# Temporary latlon check	
								rectangle = Polygon([(lon_check    , lat_check    ),
													(lon_check+self.grid_res, lat_check    ),
													(lon_check+self.grid_res    , lat_check+self.grid_res),
													(lon_check, lat_check+self.grid_res) ])
								area_A = ellipse.intersection(rectangle).area
								if area_A > 0:
									nh3_mean = 0
									if self.is_sd:
										# SD Filter
										ref_lat, ref_lon = remap_to_new_grid(lat_check, lon_check, self.ref_bounds, self.ref_data)
										nh3_mean = self.nh3_ref[ref_lat, ref_lon]
										if values[0] < nh3_mean + 3*SD and values[0] > nh3_mean - 3*SD:
											progress = True
										else:
											progress = False
									else:
										progress = True
									if progress:
										for iv, value in enumerate(values):
											error = errors[iv]
											if not np.isnan(value) and not np.isnan(error):
												if iv == 1:
													# Adjustment
													value = values[0] - (values[2]*SCALE + OFFSET)
													error = errors[0]
												self.calcdatae[iv][0][data_lat_check][data_lon_check] += area_A/(error*uncertainty_factor)
												self.calcdatae[iv][1][data_lat_check][data_lon_check] += (value*area_A)/(error*uncertainty_factor)
												self.calcdatae[iv][2][data_lat_check][data_lon_check] += 1

										
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
	   basedata[np.isnan(basedata)]=0
	   
	   return basedata
		
	def add_frames_to_oversample(self,os_file):
	    # Takes single frame basedata
	    # Iterate through frame for each pixel
	    # Pass pixel to oversample pixel
	    # Receive gridded data
	    
		os_ncf = Dataset(os_file,'r',format="NETCDF4")
		
		lat       = np.array(os_ncf['lat'])
		lon       = np.array(os_ncf['lon'])
		dflag     = np.array(os_ncf['day_flag'])
		
		nh3       = np.array(os_ncf['nh3'])
		nh3_err   = np.array(os_ncf['nh3_err'])

		tpw       = np.array(os_ncf['tpw'])
		tpw_err   = np.array(os_ncf['tpw'])
		tpw_base  = np.array(os_ncf['tpw_base'])
		tpw_base_err = tpw_base * abs((tpw_err/tpw))
		
		scl       = np.array(os_ncf['theta'])
		dt1km     = np.array(os_ncf['dt1000'])
		theta     = np.array(os_ncf['scanlinenum'])
		
		scl[scl is ma.masked] = 1
		
		
		dt1000_fval = 0
		# Application of filters
		
		if self.daynight == 'd':
			zen_filter = dflag == 1
		elif self.daynight == 'n':
			zen_filter = dflag == 0
		else:
			zen_filter = [True for f in range(len(dflag))]
			
		tpw_filter_h = tpw < 40
		tpw_filter_l = tpw > 5
		
		if dt1000_fval != 0:
			dt1000_filter = dt1km > dt1000_fval
			full_filter = zen_filter & tpw_filter_h & tpw_filter_l & dt1000_filter # Ignore for no filtering
		else:
			full_filter = zen_filter & tpw_filter_h & tpw_filter_l
		# All filters applied
		
		index_array = np.array([i for i in range(len(lat))])
		print('Total Number of pixels: {}'.format(len(index_array)))
		index_array = index_array[full_filter]
			
		print('{} pixels after filtering'.format(len(index_array)))
		idx = 0
		while idx < len(index_array):
			index = index_array[idx]
			if idx %1000 == 0:
				print(idx, index)
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
			
			#self.grid_squares       = np.concatenate(self.grid_squares, np.array(ncfile['nh3-values']),axis=2)
			
			# Combine new grid squares
			
	#def do_sd_filtering(self):
		# Grid Squares has dimensions X,Y,5/50
		# Calculate mean of each square
		#grid_means = np.nanmean(self.grid_squares,axis=2)
		#grid_maxes = grid_means + SD*3
		#grid_mins = grid_means - SD*3
		
		#max_filter = self.grid_squares > grid_maxes
		#min_filter = self.grid_squares < grid_mins
		
		#self.calcdatae[:][:][max_filter & min_filter] = np.nan # Set above/below SD range to nan
		#self.grid_squares[max_filter & min_filter] = np.nan
		
		#self.grid_sds = np.nanstd(self.grid_squares,axis=2) # Recompute standard deviations
		
		
	def addto_calcdata(self,calcdata, newdata0, newdata1,label=False):
		#nd0_filter = np.isnan(newdata0)
		#nd1_filter = np.isnan(newdata1)
		#nd0 = newdata0[np.logical_not(nd0_filter)]
		#nd1 = newdata1[np.logical_not(nd1_filter)]
		#newdata0[nd0_filter | nd1_filter] = 0
		#newdata1[nd0_filter | nd1_filter] = 0
		#ndfilter = nd0_filter | nd1_filter
		
		#print(np.nanmax(newdata0), np.nanmin(newdata0))
		#print(np.nanmax(newdata1), np.nanmin(newdata1))
		
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
			#self.do_sd_filtering()	
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
			#self.do_sd_filtering()
			for iv in range(len(self.variables)):
				self.basedatae[iv] = self.average_oversample(self.calcdatae[iv])
			print('CALC: {}'.format(self.basedatae[0][1470][1670]))
			self.save_to_nc()
			
variables = ['nh3','nh3tpwvza','tpw','tpwb']
longs = ['Oversampled Ammonia',
				'Oversampled Adjusted Ammonia',
				'Oversampled TPWVZA',
				'Oversampled TPW']
units = ['ppbv','ppbv','mm','mm']
				
def main(date):

	l3f = l3_file(date)	

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
