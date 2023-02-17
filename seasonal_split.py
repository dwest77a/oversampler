## Graph dataEditor program
import math
import numpy as np
import os
import sys

sys.path.append('/home/users/dwest77/Documents/std_py')
import file_import as fm
import find_files as ff

date_sets = [['3','4','5'],['6','7','8'],['9','10','11'],['12','1','2']]
date_labels = ['35','68','911','122']

def split_data_seasonal(input_file, output_file, accepted_dates):  
	f = open(input_file, 'r')
	instring = f.readlines()
	f.close()
	outstring = instring[0]
	for index in range(1,len(instring)):
		line = instring[index].replace('\n','')
		line = line.replace(' ','').split(',')
		mth = line[0].split('/')[0]
		if mth in accepted_dates:
			outstring += instring[index]
	os.system('touch {}'.format(output_file))
	f = open(output_file, 'w')
	f.write(outstring)
	f.close()
		
def create_seasonal_data():
	base_path = '/home/users/dwest77/Documents/CrisIasi/Ofiles'
	## Find files
	in_files = ff.list_files(base_path, starts='iasi',ends='112g.txt')
	# /iasi_ecv_daynight_landsea_yr1yr2_mth1mth2g.txt
	for input_file in in_files:
		for index in range(0,4):
			output_file = input_file.replace('112',date_labels[index])
			split_data_seasonal(input_file, output_file, date_sets[index])
	
create_seasonal_data()
