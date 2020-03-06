#!/usr/bin/env python3

"""
Script to compare variables from two netcdf files.
"""

import argparse
from netCDF4 import Dataset
import numpy as np

def main():
    parser = argparse.ArgumentParser(description="Compare variables from two netcdf files. If the variables agree to within tolerence, this prints T in the file tmpdata.dat, otherwise it prints F.")
    parser.add_argument('file1', help='Path to first GS2 output .nc file (with .out.nc extension).')
    parser.add_argument('var1', help='Name of variable in file1')
    parser.add_argument('file2', help='Path to second GS2 output .nc file (with .out.nc extension).')
    parser.add_argument('var2', help='Name of variable in file2')

    args = parser.parse_args()


    d1 = Dataset(args.file1)
    d2 = Dataset(args.file2)
    var1 = np.array(d1.variables[args.var1])
    var2 = np.array(d2.variables[args.var2])

    if not (var1 == var2).all():
            exit(1)
    else:
            exit(0)

if __name__ == '__main__':
    main()
