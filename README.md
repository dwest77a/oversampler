# Oversampler Routine
Repository Collated 16/02/2023

# Instructions for first use
 - clone the repository
 - create a virtual environment `python -m venv /path/to/venv`
 - activate the virtual environment `source /path/to/venv/bin/activate`
 - install dependencies `pip install -r requirements.txt`

# Known issues
With outdated files, new file changes are as follows:
```
sys.path.append('/home/users/dwest77/Documents/std_py')
import find_files as ff
import pmath as pm
```
This syntax block in the imports has been replaced with pyanalysis library:
```
from pyanalysis import fileIO as ff
from pyanalysis import datasetMath as dm
```

Please check https://github.com/dwest77a/pyanalysis for the latest scripts under pyanalysis to reference.

Note: Current visualisation functions use matplotlib basemaps which require backdating. Visualisation scripts thus no longer function properly, so please comment out any mentions to mpl_toolkits in this repository.

For current development, focus on base_validation (l2basefilter.py) and oversampler (ostool.py) scripts

# Contents:

## Refs
.nc files for reference (elevation, land/sea values etc.)

## Utils
Useful scripts for analysing results etc.
 - check_errs.py for checking error contents
 - elev_gridder.py creates the reference grid by averaging values from multiple files to create a single .nc file with elevation data saved as `Estimated Surface Elevation`
 - land_gridder.py creates the reference file for land/sea values

# Oversampling Pipeline

## 1. Base Validation
Scripts for filtering L2 data
 - l2basefilter.py is the main script for applying filters to given data and saving. Requires high computation time to run through all files.
 - l2basewrapper.py can be used to push filtering to the batch job system.

## 2. Oversampling
Scripts for performing main oversampling method to L2 data to produce daily gridded files.
 - ostool.py is the main script
 - ostool_wrapper.py can be used to push oversampling to the batch job system

## 3. Gridding
Scripts for performing averaging to concatenate daily gridded files to monthly files
 - l2gridding.py is the main script (with variations for different purposes)
 - l2osgridding.py will create monthly averaged files from the oversampled files.

## 4. Graphbinning
Contains a script for binning l2 graph data - needs tidying so ignore for now.

Note: You can use scdiff.py to test differences between files.
Usage: Create a new script, `from pyanalysis.scdiff import difference` then run difference('/path/to/dir/', file1, file2, write=True)
