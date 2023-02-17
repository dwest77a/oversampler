import os
import sys

if not os.path.isdir('outs'):
    os.makedirs('outs')
if not os.path.isdir('errs'):
	os.makedirs('errs')
if not os.path.isdir('jbs_sbatch'):
	os.makedirs('jbs_sbatch')


mystring = '#!/bin/bash\n'
mystring += '#SBATCH --partition=short-serial,short-serial-4hr\n'
mystring += '#SBATCH --job-name=ammonial3osfilemonthly\n'
mystring += '#SBATCH -o outs/%j.out\n' 
mystring += '#SBATCH -e errs/%j.err\n' 
mystring += '#SBATCH --time=30:00\n'
mystring += '#SBATCH --time-min=30:00\n'
mystring += '#SBATCH --mem=1G\n'

year = '2015'
mystring += '#SBATCH --array=1-12\n'
mystring += 'module add jaspy\n'
mystring += 'python ostool_monthly.py '+year+'${SLURM_ARRAY_TASK_ID}n\n'
os.system('touch jbs_sbatch/jbs'+year+'.sbatch')
f = open('jbs_sbatch/jbs'+year+'.sbatch','w')
f.write(mystring)
f.close()
os.system('sbatch /home/users/dwest77/Documents/oversampler/jbs_sbatch/jbs'+year+'.sbatch')
print('jobs submitted')
