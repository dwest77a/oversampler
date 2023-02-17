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

instrt = 'cris'
instrt_long = instrt

time = '40:00'
min_time = '10:00'
BATCH = True
VREC = '/gws/pw/j05/rsg_share/public/transfer/barry/6ad76b04-b3df-11eb-a4d3-024ad8e814ad/{}_l2/gridded_l2/{}/'

bash_string = '''#!/bin/bash
#SBATCH --partition=ukcds
#SBATCH --account=rsgcds
#SBATCH --job-name=l2gridding
#SBATCH -o grid_outs/%j.out
#SBATCH -e grid_errs/%j.err
#SBATCH --time={}
#SBATCH --time-min={}
#SBATCH --mem=2G
'''.format(time, min_time)

options, operands = getopt(sys.argv[1:], "", ["batch="])
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

jobs_folder = '/home/users/dwest77/Documents/l2_l3_edit/gridding/jbs_sbatch/'
for filter_val in [-9999]:
	subversion = record_subversion(version)
	for year in range(year_init, year_fin+1):
		for month in range(month_init, month_fin+1):
			if operands[0] == 'y':
				for dn in ['day','night']:
					ymd = str(year)+format(month,'02d')
					monthstring = bash_string
					monthstring += 'module add jaspy\n'
					monthstring += 'python l2gridding_rep.py {} {} {} {} {} {} {}'.format(year, format(month,'02d'), dn, instrt, version, subversion, filter_val)
					os.system('touch '+jobs_folder+'jbs'+ymd+version+'.sbatch')
					f = open(jobs_folder+'jbs'+ymd+version+'.sbatch','w')
					f.write(monthstring)
					f.close()
					os.system('sbatch '+jobs_folder+'jbs'+ymd+version+'.sbatch')
			
		
