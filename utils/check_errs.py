import os
import sys

sys.path.append('/home/users/dwest77/Documents/std_py')
import file_import as fm
import find_files as ff

errs = 'filter_errs/'
err_files = ff.list_files(errs)
num_errs = 0
for efile in err_files:
	f = open(efile,'r')
	arr = f.readlines()
	elen = len(arr)
	if elen > 1:
		num_errs += 1
	f.close()
		
print('Errors: ',num_errs)
