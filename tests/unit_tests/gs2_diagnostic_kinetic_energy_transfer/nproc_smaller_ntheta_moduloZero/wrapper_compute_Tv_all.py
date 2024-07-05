#! /usr/bin/env python3

import subprocess
import os

input_filename = os.path.join(os.getcwd(), "test_kinetic_energy_transfer.out.nc")
it_start = -2
it_end = -1
nproc = 1

subprocess.run(['python', 'compute_Tv_theta.py', f"{input_filename}", f"{it_start}", f"{it_end}", f"{nproc}"])
subprocess.run(['python', 'compute_Tv_theta_kxsource.py', f"{input_filename}", f"{it_start}", f"{it_end}", f"{nproc}"])
subprocess.run(['python', 'compute_Tv_theta_kxtarget.py', f"{input_filename}", f"{it_start}", f"{it_end}", f"{nproc}"])
subprocess.run(['python', 'compute_Tv_theta_kysource.py', f"{input_filename}", f"{it_start}", f"{it_end}", f"{nproc}"])