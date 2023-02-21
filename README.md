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

### Base Validation
Scripts for filtering L2 data
 - elev_gridder.py creates the reference grid by averaging values from multiple files to create a single .nc file with elevation data saved as `Estimated Surface Elevation`
 - land_gridder.py creates the reference file for land/sea values

 