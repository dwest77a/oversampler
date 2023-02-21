# Oversampler Routine
Repository Collated 16/02/2023

# Instructions for first use
 - clone the repository
 - create a virtual environment `python -m venv /path/to/venv`
 - activate the virtual environment `source /path/to/venv/bin/activate`
 - install dependencies `pip install -r requirements.txt`

## Contents:

### Refs
.nc files for reference (elevation, land/sea values etc.)

### Utils
Useful scripts for analysing results etc.
 - check_errs.py for checking error contents
 - elev_gridder.py creates the reference grid by averaging values from multiple files to create a single .nc file with elevation data saved as `Estimated Surface Elevation`
 - land_gridder.py creates the reference file for land/sea values

### Base Validation
Scripts for filtering L2 data
 - l2basefilter.py is the main script for applying filters to given data and saving. Requires high computation time to run through all files.
 - l2basewrapper.py can be used to push filtering to the batch job system.

### Graphbinning
Contains a script for binning l2 graph data - needs tidying so ignore for now.

Note: You can use scdiff.py to test differences between files.
Usage: Create a new script, `from pyanalysis.scdiff import difference` then run difference('/path/to/dir/', file1, file2, write=True)

 