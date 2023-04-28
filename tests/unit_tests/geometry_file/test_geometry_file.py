#!/usr/bin/env python3
import sys
sys.path.append('../../../scripts/')
from theta_grid_file_parser import write_theta_grid_data, read_gs2_netcdf_into_dictionary
from netCDF4 import Dataset
import numpy as np
import os
import subprocess
import pathlib

def get_test_directory() -> pathlib.Path:
    """Get the directory of this test file"""
    return pathlib.Path(__file__).parent

def get_gs2_path() -> pathlib.Path:
    """Returns the absolute path to the gs2 executable

    Can be controlled by setting the GS2_EXE_PATH environment variable
    """
    default_path = get_test_directory() / "../../../bin/gs2"
    gs2_path = pathlib.Path(os.environ.get("GS2_EXE_PATH", default_path))
    return gs2_path.absolute()

def run_gs2(gs2_path: str, input_filename: str):
    """Run gs2 with a given input file"""
    subprocess.run([gs2_path, input_filename]).check_returncode()

#Run GS2 on eik_grid.in
run_gs2(get_gs2_path(), 'eik_grid.in')

#Read in data
eik_data = read_gs2_netcdf_into_dictionary('eik_grid.out.nc')

#Write out to a grid file
write_theta_grid_data('grid.file', eik_data)

#Run GS2 on grid_file.in
run_gs2(get_gs2_path(), 'grid_file.in')

#Read in data
grid_data = read_gs2_netcdf_into_dictionary('grid_file.out.nc')

#Check everything is close
passing = True
for k in grid_data.keys():
    passing = passing and np.allclose(grid_data[k], eik_data[k])

if not passing:
    print("Failed")
    sys.exit(1)
else:
    print("Passed")
    sys.exit(0)
