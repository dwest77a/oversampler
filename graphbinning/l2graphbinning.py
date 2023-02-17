import os
import sys
from datetime import datetime
sys.path.append('/home/users/dwest77/Documents/ecvatool/config')
from defaultpci import default_plotconfig
from getopt import getopt

import numpy as np
import math
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

time = '40:00'
min_time = '10:00'
BATCH = True

instrt_long='iasi'
instrt = 'iasi'


VREC = '/gws/pw/j05/rsg_share/public/transfer/barry/6ad76b04-b3df-11eb-a4d3-024ad8e814ad/{}_l2/graphbinned_l2/{}/'

bash_string = '''#!/bin/bash
#SBATCH --partition=short-serial,short-serial-4hr
#SBATCH --job-name=l2graphbin
#SBATCH -o outs/%j.out
#SBATCH -e errs/%j.err
#SBATCH --time={}
#SBATCH --time-min={}
#SBATCH --mem=1G
'''.format(time, min_time)

#options, operands = getopt(sys.argv[1:], "", ["year=","month=","fdex="])

#year = operands[0]
#month = operands[1]
#fdex = operands[2]

#version = record_version()
filename_base = instrt + '_graph_l2b_{}_{}_'

filters = [ [True, 0, True, 0, filename_base.format('day+night', 'land+sea')],
						[False, 0, True, 0, filename_base.format('day+night','sea')],
						[False, 1, True, 0, filename_base.format('day+night','land')],
						[True, 0, False, 1, filename_base.format('day','land+sea')],
						[False, 0, False, 1, filename_base.format('day','sea')],
						[False, 1, False, 1, filename_base.format('day','land')],
						[True, 0, False, 0, filename_base.format('night','land+sea')],
						[False, 0, False, 0, filename_base.format('night','sea')],
						[False, 1, False, 0, filename_base.format('night','land')] ]

def write_to_nc(bin_values, bin_errs, file_base, filename):
	if not os.path.isdir(file_base):
		os.makedirs(file_base)
	if not os.path.isfile(file_base+filename):
		os.system('touch {}'.format(file_base+filename))
		
	ncf_new = Dataset(file_base+filename,'w',format='NETCDF4')
	
	track_dim = ncf_new.createDimension('track',len(bin_values[0]))
	
	
	nh3 = ncf_new.createVariable('nh3', np.float32, ('track',))
	nh3.long_name = 'Ammonia Concentration'
	nh3.units = 'ppbv'
	nh3[:] = bin_values[1]
	
	nh3_err = ncf_new.createVariable('nh3_err', np.float32, ('track',))
	nh3_err.long_name = 'Ammonia Concentration binned error'
	nh3_err.units = 'ppbv'
	nh3_err[:] = bin_errs[1]
	
	tpw = ncf_new.createVariable('tpw', np.float32, ('track',))
	tpw.long_name = 'Total Precipitable Water vapour'
	tpw.units = 'mm'
	tpw[:] = bin_values[0]
	
	tpw_err = ncf_new.createVariable('tpw_err', np.float32, ('track',))
	tpw_err.long_name = 'Total Precipitable Water vapour binned error'
	tpw_err.units = 'mm'
	tpw_err[:] = bin_errs[0]
	
	bin_pop = ncf_new.createVariable('bin_pop', np.float32, ('track',))
	bin_pop.long_name = 'Number of values per bin'
	bin_pop.units = ''
	bin_pop[:] = bin_values[2]
	
	ncf_new.close()

year_init = 2015
year_fin = 2020
month_init = 1
month_fin = 12

def record_subversion(version): # Subversion
	versions_rec = VREC.format(instrt_long, version)
	if not os.path.isdir(versions_rec):
		os.makedirs(versions_rec)
	versions_rec += 'subversions.txt'
	if not os.path.isfile(versions_rec):
		os.system('touch {}'.format(versions_rec))
		
	versions = []
	f = open(versions_rec,'r')
	c = f.readlines()
	rec_string = ''
	if len(c) == 0:
		rec_string = 'Version ID	description\n'
	else:
		for line in c:
			rec_string += line
			versions.append(line.split("	")[0])
	f.close()
	isvalid = False
	while not isvalid:
		if versions != []:
			print('Latest Version: {}'.format(versions[len(versions)-1]))
		new_version = input('Version ID: ')
		if new_version not in versions:
			isvalid = True
		elif new_version == 'replace':
			return versions[len(versions)-1]
		else:
			print('Warning: Version number already in use')
			clob = input('Overwrite possible existing data? (y/n): ')
			if clob == 'y':
				isvalid = True
	if new_version != 'replace':
		desc = input('Description: ')
	
		rec_string += '{}	{}\n'.format(new_version, desc)
		f = open(versions_rec,'w')
		f.write(rec_string)
		f.close()
	else:
		new_version = versions[len(versions)-1]
	return new_version

def main(var_range, binned, indep, sv, v):
	
	nbins = 50
	x_range = var_range[0] - var_range[1]
	binwidth = x_range/nbins
	
	binarr = [ [var_range[1] + binwidth/2 + binwidth*ibin for ibin in range(nbins)] for i in range(len(filters)) ]
	
	indep_errs_arr = [ [ [] for ibin in range(nbins)] for i in range(len(filters)) ]
	
	bin_errs_arr = [ [ [] for ibin in range(nbins)] for i in range(len(filters)) ]
	
	indepsumarr = [ [ 0 for ibin in range(nbins)] for i in range(len(filters))]
	# Bin data within this range of tpw values
	# Open each file within this set of time
	# Perform filters
	# Assemble set of data to be binned
	# Bin data
	
	# Save binned data
	for year in range(year_init, year_fin+1):
		print('year: {}'.format(year))
		for mth in range(month_init, month_fin+1):
			month = format(mth,'02d')
			uk_bounds = [49.5, 60.5, -10, 2]
			
			filtered_points = []
			
			l2_file_dir = '/gws/pw/j05/rsg_share/public/transfer/barry/6ad76b04-b3df-11eb-a4d3-024ad8e814ad/{}_l2/valid_l2a/{}/{}/{}'.format(instrt_long, v,year, month)
			files_array = ff.list_files(l2_file_dir, starts='', ends='')
			print('found {} for {}/{}'.format(len(files_array), month, year))
			for idx, fileN in enumerate(files_array):
				# Open file to view contents
				ncf = Dataset(fileN,'r',format='NETCDF4')
				
				lat = np.array(ncf['lat'])
				lon = np.array(ncf['lon'])
				tpwb = np.array(ncf['tpw_base'])
				dayf = np.array(ncf['day_flag'])
				landf = np.array(ncf['land_flag'])
				binned_arr = np.array(ncf[binned])
				indep_arr = np.array(ncf[indep])
				
				num_arr = np.array( [ i for i in range(len(binned_arr))])
				
				spatial1 = np.array(lat > uk_bounds[0])
				spatial2 = np.array(lat < uk_bounds[1])
				spatial3 = np.array(lon > uk_bounds[2])
				spatial4 = np.array(lon < uk_bounds[3])
				
				tpw1 = np.array(tpwb < 70)
				for fdex, filterN in enumerate(filters):
					isday = np.array((filterN[3] == dayf) | (filterN[2]))
					island = np.array((filterN[1] == landf) | (filterN[0]))
				
					all_cons = spatial1 & spatial2 & spatial3 & spatial4 & tpw1 & isday & island
				
					bindex_arr = np.array( np.round((binned_arr - var_range[1])/binwidth) )
					indep_reduced = indep_arr[all_cons]
					binned_reduced = binned_arr[all_cons]
					bda = bindex_arr[all_cons]
					nums = num_arr[all_cons]
				
					for bindex, bd in enumerate(bda):
						if bd < nbins and bd >= 0:
							indep_errs_arr[fdex][int(bd)].append(indep_reduced[bindex])
							bin_errs_arr[fdex][int(bd)].append(binned_reduced[bindex])
							indepsumarr[fdex][int(bd)] += indep_reduced[bindex]
	for fdex in range(len(filters)):
		filterN = filters[fdex]
		bin_sizes = []
		x_err = []
		y_err = []
		avg_arr = []
		for ibin in range(nbins):
			bin_sizes.append(len(indep_errs_arr[fdex][ibin]))
			if len(indep_errs_arr[fdex][ibin]) > 0:
				indep_mean = indepsumarr[fdex][ibin] / len(indep_errs_arr[fdex][ibin])
				
				avg_arr.append(indep_mean) # Compute bin average
				
				bin_mean = np.sum(bin_errs_arr[fdex][ibin])/len(bin_errs_arr[fdex][ibin])
				
				x_err.append(pm.mean_error(bin_mean, bin_errs_arr[fdex][ibin]))
				
				y_err.append(pm.mean_error(indep_mean, indep_errs_arr[fdex][ibin]))
				

			else:
				avg_arr.append(np.nan)
				x_err.append(np.nan)
				y_err.append(np.nan)
			# Calculate binned errors
		file_base = '/gws/pw/j05/rsg_share/public/transfer/barry/6ad76b04-b3df-11eb-a4d3-024ad8e814ad/{}_l2/graphbinned_l2/{}/'.format(instrt_long, v)
		filename = filterN[4] + '{}{}_{}{}_{}_{}_{}_{}.ncf'.format(year_init,year_fin, month_init,month_fin, fdex, binned, v, sv)
		write_to_nc([binarr[fdex], avg_arr, bin_sizes], 
		            [x_err, y_err],file_base, filename)
						
						
				
				# For filter in filters
				# If filter is met, add whole point data to filtered points
version = 'v00.03'
subversion = record_subversion(version)

main([5, 40], 'tpw','nh3',subversion, version)

