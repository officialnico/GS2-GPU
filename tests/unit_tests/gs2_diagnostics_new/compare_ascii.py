#!/usr/bin/env python3
import argparse
import numpy as np

parser = argparse.ArgumentParser(description="Compare variables from two netcdf files.")

parser.add_argument('file1', help='Path to first GS2 output .nc file (with .out.nc extension).')
parser.add_argument('file2', help='Path to second GS2 output .nc file (with .out.nc extension).')

args = parser.parse_args()

var1 = np.loadtxt(args.file1)
var2 = np.loadtxt(args.file2)

text_file = open("tmpdata.dat", "w")
#Only compare the last N lines in the second file
#to account for differences between old and new diagnostics
if not np.allclose(var1, var2[-var1.shape[0]:,:]):
	text_file.write("F")
else:
	text_file.write("T")
text_file.close()

