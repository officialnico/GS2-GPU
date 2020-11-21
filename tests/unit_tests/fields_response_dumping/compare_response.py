#!/usr/bin/env python3
import argparse
from netCDF4 import Dataset
import numpy as np
from glob import glob
import sys

parser = argparse.ArgumentParser(description="Compare variables from two netcdf files.")

parser.add_argument('dir1', help='Path to first GS2 response directory.')
parser.add_argument('dir2', help='Path to second GS2 response directory.')

args = parser.parse_args()

def get_sorted_data(dir_in):
    files = glob(dir_in+"/*")
    ik_vals = [int(x.split("_ik_")[1].split("_is_")[0]) for x in files]
    is_vals = [int(x.split("_is_")[1].split("_dt_")[0]) for x in files]
    data = [Dataset(f).variables['response'] for f in files]
    return data, ik_vals, is_vals, files

data1,ik1,is1,f1 = get_sorted_data(args.dir1)
v1 = [(data1[i], ik1[i], is1[i], f1[i]) for i in range(len(data1))]
data2,ik2,is2,f2 = get_sorted_data(args.dir2)
v2 = [(data2[i], ik2[i], is2[i], f2[i]) for i in range(len(data2))]

results_match = True
for d1, d2 in zip(v1, v2):
    this_case = np.allclose(d1[0], d2[0])
    if not this_case:
        print("Fail : {a}, {b}".format(a=d1[-1], b=d2[-1]))
    results_match = results_match and this_case

if not results_match:
    sys.exit(1)
else:
    sys.exit(0)
