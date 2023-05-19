#!/usr/bin/env python3
import argparse
from netCDF4 import Dataset
import numpy as np

parser = argparse.ArgumentParser(description="Compare variables from two netcdf files.")

parser.add_argument(
    "file1", help="Path to first GS2 output .nc file (with .out.nc extension)."
)
parser.add_argument("var1", help="Name of variable in file1")
parser.add_argument(
    "file2", help="Path to second GS2 output .nc file (with .out.nc extension)."
)
parser.add_argument("var2", help="Name of variable in file2")

args = parser.parse_args()

are_same = False
try:
    d1 = Dataset(args.file1)
    d2 = Dataset(args.file2)

    var1 = np.array(d1.variables[args.var1])
    var2 = np.array(d2.variables[args.var2])

    are_same = np.allclose(var1, var2)
finally:
    with open("tmpdata.dat", "w") as text_file:
        text_file.write("T" if are_same else "F")
