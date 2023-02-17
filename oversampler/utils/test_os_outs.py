## NetCDF4
from netCDF4 import Dataset

## Matplotlib Packages
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Path, PathPatch
from matplotlib import colors
import matplotlib.widgets as wg
import matplotlib as m
import matplotlib.ticker as ticker
m.use('TkAgg') ## Faster rendering

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
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
	import find_files as ff
	import pmath as pm
	import output_data as od
except: 
	print('ImportError: Std_py library missing: requires find_files.py, pmath.py and output_data.py as minimum')
	sys.exit()
	
	
base = '/gws/pw/j05/rsg_share/public/transfer/barry/6ad76b04-b3df-11eb-a4d3-024ad8e814ad/cris_l2/oversampled_l2/v00.05/{}/{}/cris_oversampled_l2_{}{}_day_v00.05_{}.nc'
	
for year in range(2015,2021):
	for mth in range(1,13):
		month = format(mth,'02d')
		try:
			
			basef = base.format(year, month, year, month,'sv34a')
			#baseg = base.format(year, month, year, month,'sv34')
			ncf = Dataset(basef, 'r',format='NETCDF4')
			
			
			#ncf_2 = Dataset(baseg, 'r',format='NETCDF4')
			
			nh3 = np.array(ncf['nh3tpwvza'])
			#nh3_2 = np.array(ncf_2['nh3tpwvza_den'])
			#print(nh3.size, nh3_2.size)
			#nh3_diff = nh3 - nh3_2
			#print(nh3_diff.size, np.count_nonzero(nh3_diff))
			#x=input()
			nh3[np.isnan(nh3)] = 0
			print(np.count_nonzero(nh3), year, month)
			if np.count_nonzero(nh3) == 0:
				dele = input('Delete '+basef+'? Y/N')
				if dele == 'Y':
					os.system('rm '+basef)
		except:
			print('failed for ',year, month)
