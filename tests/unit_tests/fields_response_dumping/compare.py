#!/usr/bin/env python3
import argparse
from netCDF4 import Dataset
import numpy as np
import sys

parser = argparse.ArgumentParser(description="Compare variables from two netcdf files.")

parser.add_argument('file1', help='Path to first GS2 output .nc file (with .out.nc extension).')
parser.add_argument('file2', help='Path to second GS2 output .nc file (with .out.nc extension).')

args = parser.parse_args()

d1 = Dataset(args.file1)
d2 = Dataset(args.file2)

var_names = ['phi2', 'apar2', 'bpar2']

results_match = True
for v in var_names:
        results_match = results_match and np.allclose(d1.variables[v], d2.variables[v])
        
if not results_match:
	sys.exit(1)
else:
	sys.exit(0)

