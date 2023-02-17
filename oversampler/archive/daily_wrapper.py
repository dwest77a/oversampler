import os
import sys
from getopt import getopt
import numpy as np

if not os.path.isdir('outs'):
    os.makedirs('outs')
if not os.path.isdir('errs'):
	os.makedirs('errs')
if not os.path.isdir('jbs_sbatch'):
	os.makedirs('jbs_sbatch')


sys.path.append('/home/users/dwest77/Documents/std_py')
import find_files as ff
from netCDF4 import Dataset

mystring = '#!/bin/bash\n'
mystring += '#SBATCH --partition=short-serial,short-serial-4hr\n'
mystring += '#SBATCH --job-name=ammonial3osfile\n'
mystring += '#SBATCH -o outs/%j.out\n' 
mystring += '#SBATCH -e errs/%j.err\n' 
mystring += '#SBATCH --time=30:00\n'
mystring += '#SBATCH --time-min=30:00\n'
mystring += '#SBATCH --mem=4G\n'

instrt_long = 'cris'
instrt = 'cris'
year_init = 2015
year_fin = 2020

month_init = 1
month_fin = 12


VREC = '/gws/pw/j05/rsg_share/public/transfer/barry/6ad76b04-b3df-11eb-a4d3-024ad8e814ad/{}_l2/oversampled_l2/{}/'
options, operands = getopt(sys.argv[1:], "", ["batch="])

dlc = []

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

version = 'v00.05'
subversion = record_subversion(version)

def setup_first(fileN):
	ncf = Dataset(fileN, 'r', format='NETCDF4')
	nh3 = np.array(ncf['nh3'])
	tpw = np.array(ncf['tpw'])
	tpwb = np.array(ncf['tpw_base'])
	
	basedata = [[nh3], [tpw], [tpwb]]
	
	ncf.close()
	return basedata
	
def addto(fileN, basedata):
	ncf = Dataset(fileN, 'r', format='NETCDF4')
	nh3_new = np.array(ncf['nh3'])
	tpw_new = np.array(ncf['tpw'])
	tpwb_new = np.array(ncf['tpw_base'])
	
	basedata[0].append(nh3_new)
	basedata[1].append(tpw_new)
	basedata[2].append(tpwb_new)
	ncf.close()
	return basedata
	
def average(basedata):
	# np.nanmean function
	avgdata = np.nanmean(basedata, axis=0)
	return avgdata
	
def write_to_nc_daily_wrapper(basedata, filename, ex_file):
	# Specific name to avoid function clashes between scripts
	# Write monthly averaged oversampled data to nc file
	
	# Extract standard lat/lon lists
	ncf = Dataset(ex_file,'r',format='NETCDF4')
	latlist = np.array(ncf['lat'])
	lonlist = np.array(ncf['lon'])
	ncf.close()
	
	# Create new file if not exists
	if not os.path.isfile(filename):
		os.system('touch {}'.format(filename))
	
	ncf = Dataset(filename, 'w', format='NETCDF4')
	
	# Average data for the three variables
	nh3base = average(basedata[0])
	tpwbase = average(basedata[1])
	tpwbbase = average(basedata[2])
	
	# Write new dimensions and variables
	lat_dim = ncf.createDimension('lat',len(latlist))
	lon_dim = ncf.createDimension('lon',len(lonlist))
		
	lat_var = ncf.createVariable('lat', np.float32, ('lat',))
	lat_var.long_name = 'latitude'
	lat_var[:] = latlist
		
	lon_var = ncf.createVariable('lon', np.float32, ('lon',))
	lon_var.long_name = 'longitude'
	lon_var[:] = lonlist
		
	nh3 = ncf.createVariable('nh3', np.float32, ('lat','lon',))
	nh3.long_name = 'Oversampled Ammonia'
	nh3.units = 'ppbv'
	nh3[:,:] = nh3base
	
	tpw = ncf.createVariable('tpw', np.float32, ('lat','lon',))
	tpw.long_name = 'Oversampled Total Precipitable Water with angle adjustment'
	tpw.units = 'mm'
	tpw[:,:] = tpwbase
	
	tpwb = ncf.createVariable('tpw_base', np.float32, ('lat','lon',))
	tpwb.long_name = 'Oversampled TPW'
	tpwb.units = 'mm'
	tpwb[:,:] = tpwbbase
	
	ncf.close()
	print('finished writing to file')
	print(filename)

days_in_month = [31,28,31,30,31,30,31,31,30,31,30,31] # 31, 28
for year in range(year_init, year_fin+1):
	for month in range(month_init, month_fin+1):
		if operands[0] == 'y':
			# Do batch processing of daily files
			
			mth = format(month, '02d')
			if True:
				day = 1
				# Assemble sbatch file contents
				ymd = str(year) + mth + format(day,'02d') + 'd'
				monthstring = mystring
				monthstring += 'module add jaspy\n'
				monthstring += 'python ostool_daily.py '+ymd+' {} {} {}\n'.format(instrt, version, subversion)
				os.system('touch jbs_sbatch/jbs'+ymd+'.sbatch')
				f = open('jbs_sbatch/jbs'+ymd+'.sbatch','w')
				f.write(monthstring)
				f.close()
				os.system('sbatch jbs_sbatch/jbs'+ymd+'.sbatch')
				x=input()
		else:
			# Do processing to combine daily files to monthly files
			mth = format(month, '02d')
			oversampled_files = '/gws/pw/j05/rsg_share/public/transfer/barry/6ad76b04-b3df-11eb-a4d3-024ad8e814ad/iasi_mhs_l2/oversampled_l2/{}/{}/{}'.format(version, year, mth)
			files_array = ff.list_files(oversampled_files,starts='iasi',ends='.nc')
			if len(files_array) > 0:
				print('Monthly average for {}/{}'.format(mth, year))
				basedata = setup_first(files_array[0])
				for index in range(1, len(files_array)):
					print(index, '/', len(files_array))
					fileN = files_array[index]
					basedata = addto(fileN, basedata)
				print('averaging')
			
				exfile = files_array[0]
				newfile = 'iasi_oversampled_l2_months_{}{}_day_{}_{}.nc'.format(year, mth, version, subversion)
				write_to_nc_daily_wrapper(basedata, oversampled_files+'/'+newfile, exfile)
			else:
				print('Skipped {}/{}'.format(mth, year))
			
print('jobs submitted')
