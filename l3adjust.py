## Level 3 File Variable Editor
## Version 1.0 - 10/05/2021 09:15
## Daniel Westwood - daniel.westwood@stfc.ac.uk

## NetCDF4
from netCDF4 import Dataset

## Matplotlib Packages
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Path, PathPatch
from matplotlib import colors
import matplotlib.widgets as wg
import matplotlib as m
m.use('TkAgg') ## Faster rendering

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
#import matplotlib.path as Path

## System
import os
import sys
from datetime import datetime
from getopt import getopt

## Shape
import shapely.affinity
from shapely.geometry import Point, Polygon

## Numpy
import numpy as np
import numpy.ma as ma
import warnings
import math

# Standard dw package tools
STD_PY_LOC = '/home/users/dwest77/Documents/std_py'
try:
    sys.path.append(STD_PY_LOC)
    import file_import as fm
    import find_files as ff
    import pmath as pm
    import NCFList as ncl
except: 
    print('ImportError: Std_py library missing: requires file_import.py, find_files.py and pmath.py as minimum')
    sys.exit()
    
# IASI Scale properties
scale = 0.007
offset = 0.11
factor = 1000

# CRIS Scale Properties
#scale = 0.0011
#offset = 0.2
    
def open_file(filename, newpath, newfile):
	if not os.path.isdir(newpath):
		os.makedirs(newpath)
	os.system('touch {}'.format(newpath+newfile))
	
	newname_and_path = newpath + newfile
	ncf_original = Dataset(filename, 'r', format='NETCDF4')
	# Get only required variables
	
	nh3 = np.array(ncf_original['nh3'])
	tpw = np.array(ncf_original['tpw']) # Added vza
	sza = np.cos ( np.array(ncf_original['satzen']) * (math.pi/180) )
	tpwsec = np.array(ncf_original['tpw']) / sza
	lat = ncf_original['latitude']
	lon = ncf_original['longitude']
	
	nh3_adj = nh3 - (tpwsec*scale + offset)/factor
	
	ncf_new = Dataset(newpath+newfile, 'w', format='NETCDF4')
	
	lat_dim = ncf_new.createDimension('lat',len(lat))
	lon_dim = ncf_new.createDimension('lon',len(lon))
		
	lat_var = ncf_new.createVariable('lat', np.float32, ('lat',))
	lat_var.long_name = 'latitude'
	lat_var[:] = lat[:]
		
	lon_var = ncf_new.createVariable('lon', np.float32, ('lon',))
	lon_var.long_name = 'longitude'
	lon_var[:] = lon[:]
		
	data_var = ncf_new.createVariable('nh3_a', np.float32, ('lat','lon',))
	data_var.long_name = 'Mean NH3 column average mixing ratio adjusted for TPWSECVZA * {} + {}'.format(scale, offset)
	data_var.units = 'ppmv'
	data_var[:,:] = nh3_adj
	ncf_original.close()
		
	ncf_new.close()
	
# List all files (iasi path)
# Convert files
# Save to 

iasi_new_dir = '/gws/pw/j05/rsg_share/public/transfer/barry/6ad76b04-b3df-11eb-a4d3-024ad8e814ad/iasi_l3/nh3_v00.03/'
cris_new_dir = '/gws/pw/j05/rsg_share/public/transfer/barry/6ad76b04-b3df-11eb-a4d3-024ad8e814ad/cris_l3/nh3_v00.03/'

iasi_pre = '/gws/pw/j05/rsg_share/public/transfer/barry/6ad76b04-b3df-11eb-a4d3-024ad8e814ad/iasi_l3/nh3_v00.01.1/'
cris_pre = '/gws/pw/j05/rsg_share/public/transfer/barry/6ad76b04-b3df-11eb-a4d3-024ad8e814ad/cris_l3/v00.01/'

#filename = iasi_pre_16 + 'ral-l3h-tqoe-iasi_mhs_amsu_metopa-tir_mw-201811-day_reg9-v1000.nc'
#filename_adj = iasi_new_dir + 'ral-l3h-tqoe-iasi_mhs_amsu_metopa-tir_mw-201811-day_reg9-v1000_v00.03_nh3_a.nc'

files_array = ff.list_files(iasi_pre, starts='ral',ends='.nc')
for index in range(0,len(files_array)):
	print(index,'/', len(files_array))
	filei = files_array[index]
	newname = filei.replace(iasi_pre,'')
	newname = newname.replace('.nc','_v00.03_nh3_a.nc')
	newpath = iasi_new_dir
	try:
		open_file(filei, newpath, newname)
	except:
		print('Error with file',index)
