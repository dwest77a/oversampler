## L2 Data Base Filter
## 	- Version 1.0 - 10/05/2021 09:15
## 	- Formatted - 28/06/2021
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
from getopt import getopt

import math

# Standard dw package tools
from pyanalysis import fileIO as ff
from pyanalysis import datasetMath as dm
#from pyanalysis import visual as vis

mmair = 28.964001
mmh2o = 18
g = 9.80665
dw = 1000
PI = 3.1415926


options, operands = getopt(sys.argv[1:], "", ["date=", "version=", "instrt=", "extension="])
# Determine options and processing method

date = operands[0]
version = operands[1]
instrt = operands[2]
try:
	extension = operands[3]
except:
	extension = ''

INPUTS = {'iasi':['/gws/pw/j05/rsg_share/public/transfer/barry/6ad76b04-b3df-11eb-a4d3-024ad8e814ad/iasi_l2/valid_l2a/v00.00/{}/{}/{}',
                  '/gws/pw/j05/rsg_share/public/projects/ims/data/lv2/output_ghg_cv9_v123_lam_nat_nfo33_fgsnwp_nobc_newbc_rbc_ram4_rem20-ncam_rcl5_raer6st_rmg_oap_cnv2_ixam7-8_o3f110_swir2/metopa/{}/{}/{}'],
          'cris':['/gws/pw/j05/rsg_share/public/transfer/barry/6ad76b04-b3df-11eb-a4d3-024ad8e814ad/cris_l2/valid_l2a/v00.00c/{}/{}/{}','']}

def calculate_dt1k(t, tsk, press):
	print('calculating dt1000')
	def z(p):
		return 16 * (3 - math.log10(p))
	
	z1 = z(852.78)
	z2 = z(878.62)
	z3 = z(865.96)
	
	t1k = t[:,91] + (t[:,92] - t[:,91]) * ( (z3 - z1)/(z2-z1) )
	dt1km     = np.array(tsk - t1k)
	
	return dt1km
	
def calculate_theta_angles(lats, lons):
	print('calculating angles')
	delta_lat = np.array(lats[1:len(lats)]) - np.array(lats[0:len(lats)-1])
	dlat = ( np.insert(delta_lat, 0, 0) + np.append(delta_lat, 0) ) /2
	
	delta_lon = np.array(lons[1:len(lons)]) - np.array(lons[0:len(lons)-1])
	dlon = ( np.insert(delta_lon, 0, 0) + np.append(delta_lon, 0) ) /2
	
	theta_angles = (180/PI)* np.arctan(dlat/dlon)
	theta_angles = np.where(np.isnan(theta_angles), 0, theta_angles)
	
	return theta_angles
	
def calculate_tpw_base(press, surf, wat):
	print('calculating tpw values')
	# Tesselate pressure and surface pressure to 2d grids
	dp_layers     = np.array(press[1:101]) - np.array(press[0:100])
	
	press_tess = np.reshape(np.tile(press, len(surf)), (len(surf),len(press)) )
	surf_tess = np.reshape(np.repeat(surf, len(press)), (len(surf), len(press)) )
	
	# Dp level tesselation into 2d grid
	dp_levels = ( np.insert(dp_layers, 0, 0) + np.append(dp_layers, 0) ) /2
	dp_levels_tess = np.reshape( np.tile(dp_levels, len(surf)), (len(surf), len(press)) )
	
	# Boolean integration filter for p vs sp
	integration_filter = press_tess < surf_tess
	
	# Water column value calculation
	wat_mmr = (mmh2o/mmair)*np.exp(wat)
	wat_cm = wat_mmr * dp_levels_tess
	
	# Apply integration filter and sum
	wat_cm_filter = np.where(integration_filter, wat_cm, 0)
	tpw_base = np.sum(wat_cm_filter, axis=1) / (g*dw*10)
	
	
	return tpw_base
	
def calculate_land_flag(lats, lons):
	lg = Dataset('/home/users/dwest77/Documents/l2_l3_edit/base_validation/land_grid.nc','r',format='NETCDF4')
	land_grid = np.array(lg['land_flag'])
	lg.close()
	
	grid_lats = ((lats + 89.5)*10).astype('int')
	grid_lats = np.where(grid_lats < 1790, grid_lats, 1789)
	grid_lats = np.where(grid_lats > -1, grid_lats, 0)
	
	grid_lons = ((lons + 179.5)*10).astype('int')
	grid_lons = np.where(grid_lons < 3580, grid_lons, 3579)
	grid_lons = np.where(grid_lons > -1, grid_lons, 0)
	
	return land_grid[grid_lats, grid_lons]
	
def write_errs(arr, arr_err):
	print('Values:',np.nanmax(arr), np.nanmin(arr))
	print('Errs:',np.nanmax(arr_err), np.nanmin(arr_err))

def get_valid_data(variables, ifile, graph_data, date):
	
	ncf = Dataset(ifile, 'r',format="NETCDF4")
	ncf_vars = {}
	print('Getting ncf info')
	write_errs(np.array(ncf['w']), np.array(ncf['w_err']))
	write_errs(np.array(ncf['t']), np.array(ncf['t_err']))
	write_errs(np.array(ncf['mgf']), np.array(ncf['mgf_err']))
	t1 = np.array(ncf['w']) > 5*np.array(ncf['w_err'])
	print(np.count_nonzero(t1), t1.size)
	
	## Extract all necessary variables for calculations ##
	press     = np.array(ncf['p'])
	vza       = np.array(ncf['satzen'])
	lats      = np.array(ncf['latitude'])
	lons      = np.array(ncf['longitude'])
	zen       = np.array(ncf['solzen'])
	cfrs      = np.array(ncf['cfr'])
	cost      = np.array(ncf['jy'])
	btd_flg   = np.array(ncf['btd_flag'])
	wat       = np.array(ncf['w'])
	surf      = np.array(ncf['sp'])
	t         = np.array(ncf['t'])
	tsk       = np.array(ncf['tsk'])
	scl       = np.array(ncf['ixt']) 
	elev      = np.array(ncf['iasi_altitude']) # Elevation L2 data
	
	for variable in variables:
		if 'mgf' in variable:
			ncf_vars[variable] = [np.array(ncf[variable])[index][1] for index in range(len(ncf[variable]))]
			ncf_vars[variable + '_err'] = [np.array(ncf[variable + '_err'])[index][1] for index in range(len(ncf[variable]))]
		else:
			ncf_vars[variable] = np.array(ncf[variable][:])
			ncf_vars[variable + '_err'] = np.array(ncf[variable + '_err'][:])
	
	ncf.close()
	
	## Filtering ##
	part_data     = [[] for var in range(len(variables)*2 + 2 + 9)]
	
	var_max_min   = [[-999,999] for var in variables]
	tpw_max_min   = [-999,999]
	
	dates_arr = np.repeat(float(str(date[0]) + str(date[1]) + str(date[2])), len(surf))
	# Set and create full filter
	print('assembling filters')
	
	#cfr_filter         = cfrs < 0.2
	cfr_filter         = cfrs < 0.05
	cost_filter        = cost < 1000
	elev_filter        = elev < 400 # Initial filter
	
	#btd_flg_filter     = btd_flg == 0 # Reconfigured for elevation 08/07/2021 15:20
	
	full_filter        = cfr_filter & cost_filter & elev_filter #btd_flg_filter
	
	## End Filtering ##
	
	## Calculations ##
	
	day_flag      = np.array(zen < 90).astype('int')
	land_flag     = calculate_land_flag(lats, lons)
	
	dt1km         = calculate_dt1k(t, tsk, press)
	theta_angles  = calculate_theta_angles(lats, lons)

	sec_arr       = np.cos ( np.array(vza) * (math.pi/180) )
	tpw           = calculate_tpw_base(press, surf, wat)
	tpwsec        = tpw / sec_arr

	## Error Calculations ##
	print('calculating errors')
	wat_col_means = np.nanmean(wat, axis=1) # Mean of all water columns
	wat_mean_tess = np.reshape( np.repeat(wat_col_means, len(wat[0])), (len(wat_col_means), len(wat[0])) )
	wat_std_errs = np.sqrt( np.sum( np.square(wat - wat_mean_tess) ) ) / len(wat[0]) # 6000 values
	
	press_mean = np.nanmean(press)
	press_std_err = np.sqrt( np.sum( np.square(press - press_mean)) ) / len(press) # 1 value
	
	tpw_fraction_errs = np.sqrt( np.square(wat_std_errs/wat_col_means) + (press_std_err/press_mean)**2 )
	tpw_std_errs = tpw_fraction_errs * tpwsec
	
	## End All Calculations ##
	
	## Save data to parts array ##
	
	# Input list of data
	# Dates, latitudes, longitudes, land_flag, day_flag, 
	# nh3, tpw, tpw_err etc. (See write_to_nc)
	# angles, scanlinenum, dt1000, tpw
	
	part_data[0] = np.append(part_data[0], dates_arr[full_filter])
	part_data[1] = np.append(part_data[1], lats[full_filter])
	part_data[2] = np.append(part_data[2], lons[full_filter])
	part_data[3] = np.append(part_data[3], land_flag[full_filter])
	part_data[4] = np.append(part_data[4], day_flag[full_filter])
	
	for vindex in range(len(variables)):
		part_data[2*vindex + 5] = np.append(part_data[2*vindex + 5], np.array(ncf_vars[variables[vindex]])[full_filter])
		part_data[2*vindex + 6] = np.append(part_data[2*vindex + 6], np.array(ncf_vars[variables[vindex] + '_err'])[full_filter])
		
		var_max_min[vindex][0] = np.nanmax(np.array(ncf_vars[variables[vindex]])[full_filter])
		var_max_min[vindex][1] = np.nanmin(np.array(ncf_vars[variables[vindex]])[full_filter])
			
	lp = len(part_data)
	part_data[lp-6] = np.append(part_data[lp-6], theta_angles[full_filter])
	part_data[lp-5] = np.append(part_data[lp-5], scl[full_filter])
	part_data[lp-4] = np.append(part_data[lp-4], dt1km[full_filter])		
	part_data[lp-3] = np.append(part_data[lp-3], tpwsec[full_filter])		
	part_data[lp-2] = np.append(part_data[lp-2], tpw_std_errs[full_filter])		
	part_data[lp-1] = np.append(part_data[lp-1], tpw[full_filter])	
	#part_data[lp-1] = np.append(part_data[lp-7], elev[full_filter])

	tpw_max_min[0] = np.nanmax(tpwsec[full_filter])
	tpw_max_min[1] = np.nanmin(tpwsec[full_filter])
	## End data saving ##
	
	## Add data to existing array ##
	
	if len(part_data[0]) == 0:
		print('No valid pixels')
		return graph_data, None, None
	else:
		print('Valid pixels: ',len(part_data[0]))
		for qindex in range(len(part_data)):
			graph_data[qindex] = np.append(graph_data[qindex], part_data[qindex])
		print('Added to dataset')
		return graph_data, tpw_max_min, var_max_min
		
	## End function get_valid_data
	
def write_to_nc(graph_data, filename, var_max_min):
	ncf_new = Dataset(filename, 'w', format='NETCDF4')
	
	track_dim = ncf_new.createDimension('track',len(graph_data[0]))
	meta_dim = ncf_new.createDimension('meta',2)
	
	date = ncf_new.createVariable('date', np.float32, ('track',))
	date.long_name = 'Date index yyyy/mm/dd'
	date.units = ''
	date[:] = graph_data[0]
	
	lat = ncf_new.createVariable('lat', np.float32, ('track',))
	lat.long_name = 'latitude'
	lat.units = 'deg'
	lat[:] = graph_data[1]
	
	lon = ncf_new.createVariable('lon', np.float32, ('track',))
	lon.long_name = 'longitude'
	lon.units = 'deg'
	lon[:] = graph_data[2]
	
	land_flag = ncf_new.createVariable('land_flag', np.float32, ('track',))
	land_flag.long_name = 'Land Flag (1 for land, 0 for sea)'
	land_flag.units = ''
	land_flag[:] = graph_data[3]
	
	day_flag = ncf_new.createVariable('day_flag', np.float32, ('track',))
	day_flag.long_name = 'Day Flag (1 for day, 0 for night)'
	day_flag.units = ''
	day_flag[:] = graph_data[4]
	
	nh3 = ncf_new.createVariable('nh3', np.float32, ('track',))
	nh3.long_name = 'Ammonia Concentration'
	nh3.units = 'ppbv'
	nh3[:] = graph_data[5]
	
	nh3_meta = ncf_new.createVariable('nh3_meta', np.float32, ('meta',))
	nh3_meta[:] = var_max_min[0]
	
	nh3_err = ncf_new.createVariable('nh3_err', np.float32, ('track',))
	nh3_err.long_name = 'Estimated Error'
	nh3_err.units = 'ppbv'
	nh3_err[:] = graph_data[6]
	
	tpw = ncf_new.createVariable('tpw', np.float32, ('track',))
	tpw.long_name = 'Total Precipitable Water vapour'
	tpw.units = 'mm'
	tpw[:] = graph_data[10]
	
	tpw_meta = ncf_new.createVariable('tpw_meta', np.float32, ('meta',))
	tpw_meta[:] = var_max_min[1]
	
	tpw_err = ncf_new.createVariable('tpw_err', np.float32, ('track',))
	tpw_err.long_name = 'Estimated Error'
	tpw_err.units = 'mm'
	tpw_err[:] = graph_data[11]
	
	tpw_base = ncf_new.createVariable('tpw_base', np.float32, ('track',))
	tpw_base.long_name = 'Nadir Total Precipitable Water vapour'
	tpw_base.units = 'mm'
	tpw_base[:] = graph_data[12]
	
	dt1000 = ncf_new.createVariable('dt1000', np.float32, ('track',))
	dt1000.long_name = 'Surface to 1000 temperature difference'
	dt1000.units = 'K'
	dt1000[:] = graph_data[9]
	
	scl = ncf_new.createVariable('scanlinenum',np.float32, ('track',))
	scl.long_name = 'Scan Line Number'
	scl[:] = graph_data[8]
	print(np.nanmax(graph_data[8]), np.nanmin(graph_data[8]))
	
	theta = ncf_new.createVariable('theta',np.float32, ('track',))
	theta.long_name = 'Pixel alignment angle'
	theta[:] = graph_data[7]

	ncf_new.close()

def main(date, version):
	
	instrt_long  = instrt
	
	base_input   = INPUTS[instrt][0]
	second_input = INPUTS[instrt][1]
	
	date_arr = [date[0:4],date[4:6],format(int(date[6:len(date)]),'02d')]
	
	units           = ['ppbv']
	variables       = ['mgf']
	override_names  = ['nh3']
	
	graph_data  = [[] for var in range(len(variables)*2 + 2 + 9)] # Elevation
	tpw_max_min = [-999,999]
	var_max_min = [[-999,999] for var in variables]
		
	if True: # Condition for later
		inpath          = base_input.format(date_arr[0], date_arr[1], date_arr[2])
		inpath_alpha    = second_input.format(date_arr[0], date_arr[1], date_arr[2])
		
		files_array     = ff.list_files(inpath, starts='ral',ends='.nc')
		files_arr_alpha = ff.list_files(inpath_alpha, starts='ral',ends='.nc')
		
		filled_arr = None
		if len(files_array) > 0:
			filled_arr = files_array
		elif len(files_arr_alpha) > 0:
			filled_arr = files_arr_alpha
			
		if filled_arr != None:
			files_array = filled_arr
			
			for ix in range(0,len(files_array)): ###
				print('File {} of {} '.format(ix+1, len(files_array)))
				graph_data, tpw_mm, var_mm = get_valid_data(variables, files_array[ix], graph_data, date_arr)
				
				# Assign Max Mins
				if tpw_mm != None:
					if tpw_mm[0] > tpw_max_min[0]:
						tpw_max_min[0] = tpw_mm[0]
					if tpw_mm[1] < tpw_max_min[1]:
						tpw_max_min[1] = tpw_mm[1]
				if var_mm != None:
					for vindex in range(len(variables)):
						if var_mm[vindex][0] > var_max_min[vindex][0]:
							var_max_min[vindex][0] = var_mm[vindex][0]
						if var_mm[vindex][1] < var_max_min[vindex][1]:
							var_max_min[vindex][1] = var_mm[vindex][1]
							
			# Kept variable tpw separate for calculations	
			units.append('mm')
			variables.append('tpw')
			override_names.append('tpw')
			var_max_min.append(tpw_max_min)
			
			filename = '{}_valid_l2a_{}_{}.nc'.format(instrt, date, version)
			filepath = '/gws/pw/j05/rsg_share/public/transfer/barry/6ad76b04-b3df-11eb-a4d3-024ad8e814ad/{}_l2/valid_l2a/{}/{}/{}/'.format(instrt_long, version, date_arr[0],date_arr[1])
			if not os.path.isdir(filepath):
				os.makedirs(filepath)
			write_to_nc(graph_data, filepath+filename, var_max_min)
		else:
			print('No files for specified date for either case - {}'.format(date))
	# Send to text file
	
date_arr = [date[0:4],date[4:6],format(int(date[6:len(date)]),'02d')]
filename = '{}_valid_l2a_{}_{}{}.nc'.format(instrt, date, version, extension)
filepath = '/gws/pw/j05/rsg_share/public/transfer/barry/6ad76b04-b3df-11eb-a4d3-024ad8e814ad/{}_l2/valid_l2a/{}/{}/{}/'.format(instrt, version, date_arr[0],date_arr[1])

if not os.path.isfile(filename+filepath):
	main(date,version)
else:
	print('skipping existing file')
