import os
import sys
from datetime import datetime
from getopt import getopt

import numpy as np
from netCDF4 import Dataset

# Standard dw package tools
STD_PY_LOC = '/home/users/dwest77/Documents/std_py'
try:
    sys.path.append(STD_PY_LOC)
    import find_files as ff
    import pmath as pm
except: 
    print('ImportError: Std_py library missing: requires file_import.py, find_files.py and pmath.py as minimum')
    sys.exit()

time = '40:00'
min_time = '10:00'
BATCH = True
instrt = 'iasi'
versions_rec = '/gws/pw/j05/rsg_share/public/transfer/barry/6ad76b04-b3df-11eb-a4d3-024ad8e814ad/{}_l2/valid_l2a/{}_l2_versions.txt'.format(instrt, instrt)

bash_string = '''#!/bin/bash
#SBATCH --partition=short-serial,short-serial-4hr
#SBATCH --job-name=l2basefilter
#SBATCH -o filter_outs/%j.out
#SBATCH -e filter_errs/%j.err
#SBATCH --time={}
#SBATCH --time-min={}
#SBATCH --mem=1G
'''.format(time, min_time)

options, operands = getopt(sys.argv[1:], "", ["batch="])
year_init = 2012
year_fin = 2015
month_init = 1
month_fin = 12

def record_version(): # Main Base Version (L2 valid)
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
	
def write_to_nc(filtered_points, max_mins,  filename):
	os.system('touch {}'.format(filename))
	ncf_new = Dataset(filename,'w',format='NETCDF4')
	
	track_dim = ncf_new.createDimension('track',len(filtered_points[0]))
	meta_dim = ncf_new.createDimension('meta',2)
	
	date = ncf_new.createVariable('date', np.float32, ('track',))
	date.long_name = 'Date index yyyy/mm/dd'
	date.units = ''
	date[:] = filtered_points[0]
	
	nh3 = ncf_new.createVariable('nh3', np.float32, ('track',))
	nh3.long_name = 'Ammonia Concentration'
	nh3.units = 'ppbv'
	nh3[:] = filtered_points[1]
	
	nh3_meta = ncf_new.createVariable('nh3_meta', np.float32, ('meta',))
	nh3_meta[:] = max_mins[0]
	
	nh3_err = ncf_new.createVariable('nh3_err', np.float32, ('track',))
	nh3_err.long_name = 'Estimated Error'
	nh3_err.units = 'ppbv'
	nh3_err[:] = filtered_points[2]
	
	tpw = ncf_new.createVariable('tpw', np.float32, ('track',))
	tpw.long_name = 'Total Precipitable Water vapour'
	tpw.units = 'mm'
	tpw[:] = filtered_points[3]
	
	tpw_meta = ncf_new.createVariable('tpw_meta', np.float32, ('meta',))
	tpw_meta[:] = max_mins[1]
	
	tpw_err = ncf_new.createVariable('tpw_err', np.float32, ('track',))
	tpw_err.long_name = 'Estimated Error'
	tpw_err.units = 'mm'
	tpw_err[:] = filtered_points[4]
	
	ncf_new.close()

version = record_version()

jobs_folder = '/home/users/dwest77/Documents/l2_l3_edit/base_validation/jbs_sbatch/'

for year in range(year_init, year_fin+1):
	for month in range(month_init, month_fin+1):
		if operands[0] == 'y':
			for day in range(1,32):
				ymd = str(year)+format(month,'02d')+format(day,'02d')
				monthstring = bash_string
				monthstring += 'module add jaspy\n'
				monthstring += 'python l2basefilter_2.py {} {} {}'.format(ymd, version, instrt)
				os.system('touch '+jobs_folder+'jbs'+ymd+version+'.sbatch')
				f = open(jobs_folder+'jbs'+ymd+version+'.sbatch','w')
				f.write(monthstring)
				f.close()
				os.system('sbatch '+jobs_folder+'jbs'+ymd+version+'.sbatch')
		
