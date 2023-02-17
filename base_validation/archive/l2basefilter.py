## Graph Creator for L2 files
## Version 1.0 - 10/05/2021 09:15
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

#from mpl_toolkits import mplot3d
#from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap
#from matplotlib.patches import Ellipse, Rectangle

from getopt import getopt

import math

# Standard dw package tools
STD_PY_LOC = '/home/users/dwest77/Documents/std_py'
try:
    sys.path.append(STD_PY_LOC)
    import file_import as fm
    import find_files as ff
    import pmath as pm
    import output_data as od
except: 
    print('ImportError: Std_py library missing: requires file_import.py, find_files.py and pmath.py as minimum')
    sys.exit()
g = 9.81
map_bounds = [0,0,0,0]
bound = False
isfilter1=True
isfilter2=False
isfilter3=True
isfilter4=True
istpwsec = True

mmair = 28.964001
mmh2o = 18
g = 9.80665
dw = 1000
PI = 3.1415926


options, operands = getopt(sys.argv[1:], "", ["date=", "version=", "instrt="])
# Determine options and processing method

date = operands[0]
version = operands[1]
instrt = operands[2]

INPUTS = {'iasi':['/gws/pw/j05/rsg_share/public/transfer/barry/6ad76b04-b3df-11eb-a4d3-024ad8e814ad/iasi_l2/valid_l2a/v00.00/{}/{}/{}',
                  '/gws/pw/j05/rsg_share/public/projects/ims/data/lv2/output_ghg_cv9_v123_lam_nat_nfo33_fgsnwp_nobc_newbc_rbc_ram4_rem20-ncam_rcl5_raer6st_rmg_oap_cnv2_ixam7-8_o3f110_swir2/metopa/{}/{}/{}'],
          'cris':['/gws/pw/j05/rsg_share/public/transfer/barry/6ad76b04-b3df-11eb-a4d3-024ad8e814ad/cris_l2/valid_l2a/v00.00c/{}/{}/{}','']}

map_example = Basemap(projection='cyl',llcrnrlat=-89.5, urcrnrlat=89.5, 
									   llcrnrlon=-179.5, urcrnrlon=179.5, lat_ts=20, resolution='l')

def determine_land(lat, lon):
	xpt, ypt = map_example(lon, lat)
	return int(map_example.is_land(xpt, ypt))

def get_valid_data(variables, ifile, graph_data, date):
	# Open file and retrieve contents
	ncf = Dataset(ifile, 'r',format="NETCDF4")
	
	ncf_vars = {}
	### Extraction of Necessary Variables ###
	tpw = []
	tpw_errs = []
	
	# ts - surface temperature
	# find location of minimum of abs(p-sp)
	
	t         = np.array(ncf['t'])
	tsk       = np.array(ncf['tsk'])
	
	def z(p):
		return 16 * (3 - math.log10(p))
	
	z1 = z(852.78)
	z2 = z(878.62)
	z3 = z(865.96)
	
	t1k = t[:,91] + (t[:,92] - t[:,91]) * ( (z3 - z1)/(z2-z1) )
	
	
	dt1km     = np.array(tsk - t1k)
	
	#fig = plt.figure()
	#fig = od.output_two_var_graph(fig, [dt1km, dt1km], graph_title='Histogram', ylabel='Number', xlabel='DT1000', do_reg=False, do_errs=False, do_hist=True, do_bins=True, nbins=30, graph_lims = [[-20,60],[-10,4000]])
	#plt.show()

	press     = np.array(ncf['p'])
	vza       = np.array(ncf['satzen'])
	lats      = np.array(ncf['latitude'])
	lons      = np.array(ncf['longitude'])
	zens      = np.array(ncf['solzen'])
	cfrs      = np.array(ncf['cfr'])
	cost      = np.array(ncf['jy'])
	btd_flg   = np.array(ncf['btd_flag'])
	wat       = np.array(ncf['w'])
	surf      = np.array(ncf['sp'])
	
	temp = np.nanmean(t, axis=0)
	
	# For oversampling specifically
	scl       = np.array(ncf['ixt'])
	# scanlinenum
	# pixel_angle
	
	# Surface at 97 index
	# 1km at 95
	
	dp_layers = np.array(press[1:101]) - np.array(press[0:100])
	
	press_mean    = pm.simple_mean(press)
	press_std_err = pm.mean_error(press_mean, press)
	
	sec_arr       = np.cos ( np.array(vza) * (math.pi/180) )
	
	# Could add other variables later
	for variable in variables:
		if 'mgf' in variable:
			ncf_vars[variable] = [ncf[variable][index][1] for index in range(len(ncf[variable]))]
			ncf_vars[variable + '_err'] = [ncf[variable + '_err'][index][1] for index in range(len(ncf[variable]))]
		else:
			ncf_vars[variable] = np.array(ncf[variable][:])
			ncf_vars[variable + '_err'] = np.array(ncf[variable + '_err'][:])
	
	### End Variable Extraction ###
	
	# Number of parts of data to be saved
	part_data  = [[] for var in range(len(variables)*2 + 2 + 9)]
	
	
	is_day = False
	day_night = True
	
	var_max_min = [[-999,999] for var in variables]
	tpw_max_min = [-999,999]
	
	print('Running calculations')
	index = 0
	while index < len(lats)-1 and day_night:
		
		# Python requires dp_levels to be overwritten with each loop
		dp_levels = ( np.insert(dp_layers, 0, 0) + np.append(dp_layers, 0) ) /2
		# Annoyingly
		
		lat = lats[index]
		lon = lons[index]
		zen = zens[index]
		dt1000 = dt1km[index]
		dayfilter = False
		# Sun filter
		
		# Cloud Fraction filter (>20% is ignored)
		if not np.isnan(cfrs[index]):
			filter1 = isfilter1 and (cfrs[index] > 0.2) # Ignore > 20%
		else:
			filter1 = False
		# Alternate Cloud Fraction filter (> 5% ignored)
		if not np.isnan(cfrs[index]):
			filter2 = isfilter2 and (cfrs[index] > 0.05) # Ignore > 5%
		else:
			filter2 = False
		# Cost filter (>1000 is ignored)
		if not np.isnan(cost[index]):
			filter3 = isfilter3 and (cost[index] > 1000) # Ignore > 1000
		else:
			filter3 = False
		# btd_flag filter (==1 ignored)
		if not np.isnan(btd_flg[index]):
			filter4 = isfilter4 and (btd_flg[index] == 1) # Ignore = 0
		else:
			filter4 = False
		# Lat/Lon filters
		if bound and (lon < map_bounds[2] or lon > map_bounds[3]):
			pass
		elif bound and (lat < map_bounds[0] or lat > map_bounds[1]):
			pass
		elif filter1 or filter2 or filter3 or filter4 or dayfilter:
			pass
		# Passed all filters
		else:

			## Calculate TPW
			
			wat_mmr = (mmh2o/mmair)*np.exp(wat[index]) # mmr at index (latlon)
			wat_cm = 0
			for idx in range(len(wat_mmr)):
				if press[idx] < surf[index]:
					wat_cm += wat_mmr[idx]*dp_levels[idx]
			tpw = wat_cm / (g*dw*10) # conversion factors from m to mm and hpa to pa incur an extra factor of 10
			print(wat[index])
			wat_mean    = pm.simple_mean(wat[index])
			wat_std_err = pm.mean_error(wat_mean, wat[index])
			
			means       = [wat_mean, press_mean]
			errs        = [wat_std_err, press_std_err]
	
			if istpwsec:
				tpwsec = tpw / sec_arr[index]
			else:
				tpwsec = tpw
				
			print(index, wat_mmr)
			x=input()
			tpw_std_err = pm.simple_errors(means, errs)* tpwsec	
			## End TPW Calculations
			
			## Pixel Theta Calculations
			idx = index
			valid_angle=False
			while not valid_angle:
				npixel_lat = lats[idx]
				npixel_lon = lons[idx]
		
				pixel_lat = lats[idx+1]
				pixel_lon = lons[idx+1]
				frac = (npixel_lat - pixel_lat) / (npixel_lon - pixel_lon)
				pixel_theta = math.atan(frac)*(180/PI)
				if not pixel_theta is ma.masked and not np.isnan(pixel_theta):
					valid_angle = True
				idx += 1
			## End Pixel Theta Calculations
			
			# Second set of filters for different nan values
			if (not tpwsec is ma.masked and not np.isnan(tpwsec)):
				isAppend = False
				for variable in variables:
					val = ncf_vars[variable][index]
					if not val is ma.masked and not np.isnan(val):
						isAppend = True
					else:
						print(val)
				if isAppend:
					# Made it through all filters
					land_flag = determine_land(lat, lon)
					
					part_data[0].append(float(str(date[0]) + str(date[1]) + str(date[2])) ) # Date
					part_data[1].append(lat) # Lat
					part_data[2].append(lon) # Lon
					part_data[3].append(land_flag) # land_flag
					part_data[4].append(int(zen < 90)) # Day flag
					
					for vindex in range(len(variables)):
						part_data[2*vindex + 5].append( ncf_vars[variables[vindex]][index])
						part_data[2*vindex + 6].append( ncf_vars[variables[vindex] + '_err'][index])
						
						if ncf_vars[variables[vindex]][index] > var_max_min[vindex][0]:
							var_max_min[vindex][0] = ncf_vars[variables[vindex]][index]
						if ncf_vars[variables[vindex]][index] < var_max_min[vindex][1]:
							var_max_min[vindex][1] = ncf_vars[variables[vindex]][index]
						
					part_data[len(part_data)-6].append(pixel_theta)
					part_data[len(part_data)-5].append(scl[index])
					part_data[len(part_data)-4].append(dt1000)
					part_data[len(part_data)-3].append(tpwsec)
					part_data[len(part_data)-2].append(tpw_std_err)
					part_data[len(part_data)-1].append(tpw)
					
					if tpwsec > tpw_max_min[0]:
						tpw_max_min[0] = tpwsec
					if tpwsec < tpw_max_min[1]:
						tpw_max_min[1] = tpwsec
				
		index += 1
	ncf.close()
	# Temp wat avg calculation
	x=input('File complete')
	if len(part_data[0]) == 0:
		print('No valid pixels')
		return graph_data, None, None
	else:
		print('Valid pixels: ',len(part_data[0]))
		for qindex in range(len(part_data)):
			if part_data[qindex] != None:
				graph_data[qindex] = graph_data[qindex] + part_data[qindex] 
		print('Added to dataset')
		return graph_data, tpw_max_min, var_max_min
	
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
	scl[:] = graph_data[7]
	
	theta = ncf_new.createVariable('theta',np.float32, ('track',))
	theta.long_name = 'Pixel alignment angle'
	theta[:] = graph_data[8]
	
	ncf_new.close()

def main(date, version):
	
	instrt_long = instrt
	
	base_input = INPUTS[instrt][0]
	second_input = INPUTS[instrt][1]
	
	date_arr = [date[0:4],date[4:6],format(int(date[6:len(date)]),'02d')]
	ymd = date
	
	units = ['ppbv']
	variables = ['mgf']
	override_names = ['nh3']
	
	graph_data = [[] for var in range(len(variables)*2 + 2 + 9)]
	tpw_max_min = [-999,999]
	var_max_min = [[-999,999] for var in variables]
		
	if True: # Condition later
		inpath = base_input.format(date_arr[0], date_arr[1], date_arr[2])
		
		inpath_alpha = second_input.format(date_arr[0], date_arr[1], date_arr[2])
		files_array = ff.list_files(inpath, starts='ral',ends='.nc')
		files_arr_alpha = ff.list_files(inpath_alpha, starts='ral',ends='.nc')
		farr = None
		if len(files_array) > 0:
			farr = files_array
		elif len(files_arr_alpha) > 0:
			farr = files_arr_alpha
		if farr != None:
			files_array = farr
			for ix in range(0,len(files_array)):
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

main(date,version)
