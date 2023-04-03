#!/usr/bin/env python3
import argparse
import numpy as np

parser = argparse.ArgumentParser(description="Compare variables from two netcdf files.")

parser.add_argument("file1", help="Path to first GS2 output .nc file (with .out.nc extension).")
parser.add_argument("file2", help="Path to second GS2 output .nc file (with .out.nc extension).")

args = parser.parse_args()
print(f"Reading {args.file1}, {args.file2}", flush=True)

# If we can't parse the file as an array of numbers, fall back to simple text comparison.
# Bear in mind that this might be fragile!
try:
    var1 = np.loadtxt(args.file1)
    var2 = np.loadtxt(args.file2)
    are_same = np.allclose(var1, var2[-var1.shape[0] :, :])
except ValueError:
    with open(args.file1) as f:
        text1 = f.read()
    with open(args.file2) as f:
        text2 = f.read()
    are_same = text1 == text2


with open("tmpdata.dat", "w") as text_file:
    if not are_same:
        text_file.write("F")
    else:
        text_file.write("T")
