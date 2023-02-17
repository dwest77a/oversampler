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

# CRIS L2 NH3 TPW Adjustment Values
#SCALE = 0.0012
#OFFSET = 0.0024

# IASI L2 NH3 TPW Adjustment Values
SCALE = 0.003
OFFSET = -0.06

def calculate_theta_angles(lats, lons):
	print('calculating angles')
	delta_lat = np.array(lats[1:len(lats)]) - np.array(lats[0:len(lats)-1])
	dlat = ( np.insert(delta_lat, 0, 0) + np.append(delta_lat, 0) ) /2
	
	delta_lon = np.array(lons[1:len(lons)]) - np.array(lons[0:len(lons)-1])
	dlon = ( np.insert(delta_lon, 0, 0) + np.append(delta_lon, 0) ) /2
	
	theta_angles = (180/PI)* np.arctan(dlat/dlon)
	theta_angles = np.where(np.isnan(theta_angles), 0, theta_angles)
	
	return theta_angles

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
	def __init__(self, date, instrt, version, subversion, time):
		
		self.grid_res = 0.1
		self.map_bounds = [-89.5, 89.5, -179.5, 179.5]
		
		self.version = version
		self.subversion = subversion
		
		self.variables = ['surfT','em741','em1140','em1290']
		self.longs = ['Surface Temperature',
						'Emissivity at 741',
						'Emissivity at 1140',
						'Emissivity at 1290']
		self.units = ['K','','','']
		self.create_grids()
		
		self.instrt = instrt.upper()
		instrt_full = instrt
		
		self.date = [date[0:4],date[4:6],format(int(date[6:len(date)-1]),'02d')]
		self.daynight = date[len(date)-1:len(date)] # Last character
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
		
		
		self.latlist = (np.arange(0,num_lats)*self.grid_res) + self.map_bounds[0]
		self.lonlist = (np.arange(0,num_lons)*self.grid_res) + self.map_bounds[2]
		lvs = len(self.variables)
		
		self.basedatae = np.reshape( np.repeat(0, num_lats*num_lons*lvs), (lvs, num_lats, num_lons)).astype('float')
		self.calcdatae = np.reshape( np.repeat(0, 3*num_lats*num_lons*lvs), (lvs, 3, num_lats, num_lons)).astype('float')
	
		print('Grids Created')
		
	def oversample_pixel(self, pixel_lat, pixel_lon, pixel_snum, values, errors, pixel_theta): # nh3 tpw tpwb
		# Takes single pixel with dimensions
		# Takes resolution of constant grid and initial coordinates
		# Determines overlaps of single pixel
		# Adds to Sum [0] and Denom [1] of 2d lat/lon array
		semi_major = 12 + (29/60)*abs(pixel_snum-60)
		semi_minor = 12 + (8/60)*abs(60-pixel_snum)
		
		ll_major = (semi_major*km2deg)/2
		ll_minor = (semi_minor*km2deg)/2
        
		circle = Point(pixel_lon, pixel_lat).buffer(1)
		ellipse = shapely.affinity.scale(circle, ll_major, ll_minor)
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
							
							rectangle = Polygon([(lon_check    , lat_check    ),
												(lon_check+self.grid_res, lat_check    ),
												(lon_check+self.grid_res    , lat_check+self.grid_res),
												(lon_check, lat_check+self.grid_res) ])
							area_A = ellipse.intersection(rectangle).area
							if area_A > 0:
								for iv, value in enumerate(values):
									error = errors[iv]
									if not np.isnan(value) and not np.isnan(error):
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
	    
		dt1000_f1 = True
		dt1000_f2 = False
		dt1000_f3 = False
	    
		os_ncf = Dataset(os_file,'r',format="NETCDF4")
		
		lat       = np.array(os_ncf['latitude'])
		lon       = np.array(os_ncf['longitude'])
		surft     = np.array(os_ncf['tsk'])
		surft_err     = np.array(os_ncf['tsk_err'])
		
		em1       = np.array(os_ncf['em'])[:,40]
		em2       = np.array(os_ncf['em'])[:,121]
		em3       = np.array(os_ncf['em'])[:,159]
		em1_err       = np.array(os_ncf['em_err'])[:,4]
		em2_err       = np.array(os_ncf['em_err'])[:,6]
		em3_err       = np.array(os_ncf['em_err'])[:,8]
		
		cfr       = np.array(os_ncf['cfr'])
		cost      = np.array(os_ncf['jy'])
		
		scl       = np.array(os_ncf['ixt'])
		theta     = calculate_theta_angles(lat, lon)
		
		scl[scl is ma.masked] = 1
		
		cfr_filter  = cfr < 0.2
		cost_filter = cost < 1000
		# No day filter

		full_filter = cfr_filter & cost_filter
		# All filters applied
		
		index_array = np.array([i for i in range(len(lat))])
		print('Total Number of pixels: {}'.format(len(index_array)))
		index_array = index_array[full_filter]
			
		print('{} pixels after filtering'.format(len(index_array)))
		idx = 0
		while idx < len(index_array):
			index = index_array[idx]
			if idx %1000 == 0:
				print(idx,index)
			values = [surft[index],em1[index],em2[index],em3[index]]
			errors = [surft_err[index],em1_err[index],em2_err[index],em3_err[index]]
			self.oversample_pixel(lat[index], lon[index], scl[index], values, errors,theta[index])
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
		# nh3 = var_class(ncfile, self.nh3calcdata)
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
		try:
			ncfile = Dataset(filein, 'r', format = 'NETCDF4')
		
			for iv in range(len(self.variables)):
				denlabel = self.variables[iv] + '_den'
				numlabel = self.variables[iv] + '_num'
				self.calcdatae[iv]   = self.addto_calcdata(self.calcdatae[iv], np.array(ncfile[denlabel]),np.array(ncfile[numlabel]))
				print(iv, np.count_nonzero(self.calcdatae[iv][0] - self.calcdatae[iv][1]), np.size(self.calcdatae[iv][0]))
		except IndexError:
			print('missing {}'.format(filein))
		
	def addto_calcdata(self,calcdata, newdata0, newdata1,label=False):
		nd0_filter = np.isnan(newdata0)
		nd1_filter = np.isnan(newdata1)
		newdata0[nd0_filter | nd1_filter] = 0
		newdata1[nd0_filter | nd1_filter] = 0
		
		calcdata[0] = calcdata[0] + newdata0
		calcdata[1] = calcdata[1] + newdata1
		return calcdata
		 
	def main(self):
		print('Main')
		print(self.timeformat)
		if self.timeformat == 'Daily':
			# Do single file
			
			# Replace this with base input
			day_base_in = '/gws/pw/j05/rsg_share/public/projects/ims/data/lv2/output_ghg_cv9_v123_lam_nat_nfo33_fgsnwp_nobc_newbc_rbc_ram4_rem20-ncam_rcl5_raer6st_rmg_oap_cnv2_ixam7-8_o3f110_swir2/metopa/'
			time = '{}/{}/{}'.format(self.date[0],self.date[1],self.date[2])
			filein = day_base_in + time
			print(filein)
			files_array     = ff.list_files(filein, starts='ral',ends='.nc')
			print(len(files_array), 'files')
			for fileN in files_array: # Add all base file data
				self.add_frames_to_oversample(fileN)
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
				if True:
					self.add_frames_to_average(out_dir+base_file_format)
				#except:
					#pass
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
			self.save_to_nc()
				
def main(date):

	l3f = l3_file(date)	

options, operands = getopt(sys.argv[1:], "", ["date=","instrt=","version=","subversion=","time="])
# Determine options and processing method

date = operands[0]
version = operands[2]
subversion = operands[3]
instrt = operands[1]
time = operands[4]
l3f = l3_file(date, instrt, version, subversion, time)
