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

import math

sys.path.append('/home/users/dwest77/Documents/std_py')
import file_import as fm
import find_files as ff

PI = 3.1415926
km2deg = 9e-3

TITLE_DICT    = { 'fontsize':'large',
                  'fontweight':'normal',
                  'color':'black',
                  'verticalalignment':'baseline',
                  'horizontalalignment':'center'}

class Oversample:
	def __init__(self, grid_res, map_bounds, variable, vshort, instrt, timeformat, year, month):
		self.grid_res = grid_res
		self.map_bounds = map_bounds
		
		self.variable = variable
		self.vshort = vshort

		
		self.instrt = instrt
		self.timeformat = timeformat
		self.year = year
		self.month = month
		
		self.oversample_main()
	
	def restart_grids(self):
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
		
	def plot_current(self):
		
		out_fig = plt.figure(figsize=(8,8))
		out_fig.set_facecolor('white')
		map_ax = out_fig.add_axes()
		
		title = 'L3 NH3 Oversampled IASI - 02/02/2020'
		
		plt.title(title,fontdict=TITLE_DICT,loc='center')
		map_bounds = self.map_bounds
		ecv_map = Basemap(projection='cyl',llcrnrlat=map_bounds[0], urcrnrlat=map_bounds[1], llcrnrlon=map_bounds[2], urcrnrlon=map_bounds[3], lat_ts=20, resolution='l')
		ecv_map.drawcoastlines()
		ecv_map.drawcountries()

		ecv_map.drawparallels(np.arange(-90.,91.,30.))
		ecv_map.drawmeridians(np.arange(-180.,181.,60.))

		ecv_map.drawmapboundary(fill_color='grey')
		
		print(len(self.lonlist), len(self.latlist))
        
		lon2d,lat2d = np.meshgrid( self.lonlist, self.latlist)
		grid_x,grid_y = ecv_map(lon2d,lat2d)
		ecv_map_color = ecv_map.pcolor(grid_x,grid_y,self.basedata, cmap = 'seismic', alpha=0.7)
		
		colorbar_ax = out_fig.add_axes()
		colorbar = plt.colorbar( ecv_map_color, cax = colorbar_ax, orientation = 'horizontal', extend='both')
            

		plt.clim(self.min,self.max)
		colorbar.ax.set_xlabel('Units: ppbv', verticalalignment='center', horizontalalignment='center')
        
		plt.show()

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
		
		npixel_lat = os_ncf['latitude'][0]
		npixel_lon = os_ncf['longitude'][0]
		
		pixel_lat = os_ncf['latitude'][1]
		pixel_lon = os_ncf['longitude'][1]
		frac = (npixel_lat - pixel_lat) / (npixel_lon - pixel_lon)
		pixel_theta = math.atan(frac)*(180/PI)
		
		pixels = []
		start_time = time.time()
		mcount = 0
		for index in range(0, len(os_ncf['latitude'])):
			if True:
				lat = os_ncf['latitude'][index]
				lon = os_ncf['longitude'][index]
				if lon < self.map_bounds[2] or lon > self.map_bounds[3]:
					pass
				elif lat < self.map_bounds[0] or lat > self.map_bounds[1]:
					pass
				else:
					val = os_ncf[self.variable][index][1]
					val_err = os_ncf[self.variable+'_err'][index][1]
					if not val is ma.masked and not val_err is ma.masked and not np.isnan(val) and not np.isnan(val_err):
						pixels.append(index)
			if not os_ncf[self.variable][index][1] is ma.masked:
				mcount += 1
		print('Running {} pixels withing range'.format(len(pixels)))
		for pdex, index in enumerate(pixels):
			if not os_ncf[self.variable][index][1] is ma.masked:
				snum = os_ncf['ixt'][index]
				if snum is ma.masked:
					snum = 1
				self.oversample_pixel(os_ncf['latitude'][index],
								os_ncf['longitude'][index],
								snum,
								(os_ncf[self.variable][index][1]),
								abs(os_ncf[self.variable+'_err'][index][1]),
								pixel_theta)
			else:
				print('ma')

		#end_time = time.time()
		#ltime = end_time - start_time
		if len(pixels) == 0:
			print('File Contains No relevant data')
		return None
		
	def assemble_file_name(self,day):
		# IASI_NH3_Daily_cv9_v12_2020_02_02.nc
		base_file_format = '{}_{}_{}_cv9_v12_{}_{}_{}.nc'.format(self.instrt, self.vshort, self.timeformat, self.year, self.month, day)
		file_dir = '/home/users/dwest77/Documents/IASI_os_l3/'
		return file_dir+base_file_format
		
	def write_to_nc(self,day):
		out_dir = self.assemble_file_name(day)
		ncfile = Dataset(out_dir, 'w', format = 'NETCDF4')
		
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
		print(len(self.basedata), len(self.basedata[0]))
		try:
			print(len(self.basedata[0][0]))
		except:
			x=1
		print(len(self.latlist), len(self.lonlist))
		data_var[:,:] = self.basedata
		
		ncfile.close()
		print('write successful')
		
	
	def oversample_main(self):
	    # Takes base location of frame files
	    # Takes year, month
	    # Define constant grid
	    # Use list_files (find_files) or other to get all files in a list
	    # Pass each file to add_frames_to_oversample
	    # Accumulate constant grid
	    # Perform average oversample
		
		base_input = '/gws/pw/j05/rsg_share/public/projects/ims/data/lv2/output_ghg_cv9_v12_lam_fgsnwp_nobc_newbc_rbc_ram2_rem20-ncam_rcl5_raer6st_rmg_oap_cnv2_o3f113_csyn_swir2/npp/{}/{}'
		#temp = '/gws/pw/j05/rsg_share/public/projects/ims/data/lv2/output_ghg_cv9_v123_lam_nat_mwnrt_nfo33_fgsnwp_nobc_newbc_rbc_ram4_rem20-ncam_rcl5_raer6st_rmg_oap_cnv2_ixam7-8_o3f110_swir2/metopa/'
		
		base_input = base_input.format(self.year,self.month)
		if self.timeformat == 'Daily':
			for day in range(1,32):
				print('    DAY: ',day)
				if True:
					full_input = base_input + f'/{day:02d}/'
					files_array = ff.list_files(full_input,starts='ral',ends='.nc')
					self.restart_grids()
					print('found {} files'.format(len(files_array)))
					if len(files_array) > 0:
			
						for ix in range(0,len(files_array)):
							print('File {} of {} '.format(ix+1, len(files_array)))
							self.add_frames_to_oversample(files_array[ix])
						self.average_oversample()
						self.write_to_nc(f'{day:02d}')
					else:
						print('{} skipped'.format(day))
				#except:
					#print('{} skipped'.format(day))
		return None
	    
test = Oversample(0.1,[40, 60, -20 ,5],'mgf','nh3','IASI','Daily', '2020', '03')
