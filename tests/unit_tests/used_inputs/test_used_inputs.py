#!/usr/bin/env python3
# This test should be run with `pytest`

import os
import pathlib
import re
import shutil
import subprocess

import netCDF4
import numpy as np


HEADER_CREATED_RE = re.compile(
    r"! Created by GS2 at (\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}.\d+(Z\+\d{2}:\d{2}))"
)
HEADER_UUID_RE = re.compile(
    r"! Run UUID: ([0-9a-fA-F]{8}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{12})"
)
HEADER_FILENAME_RE = re.compile(r"! Input parameters as used by [a-zA-Z0-9_.]+")


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


def run_gs2(gs2_path: str, input_filename: str) -> int:
    """Run gs2 with a given input file"""
    subprocess.run([gs2_path, input_filename]).check_returncode()


def copy_input_file(input_filename: str, destination):
    """Copy input_filename to destination directory"""
    shutil.copyfile(get_test_directory() / input_filename, destination / input_filename)


def compare_two_output_files(file1, file2, tolerance=1.0e-6):
    """Compares all the variables in file1 and file2 using `numpy.allclose`"""

    file1_data = netCDF4.Dataset(file1)
    file2_data = netCDF4.Dataset(file2)

    assert file1_data["t"][-1] == file2_data["t"][-1]

    # Remove variables that we know can change from run to run. These
    # are stored as attributes in new diagnostics
    variables = list(file1_data.variables.keys())
    try:
        variables.remove("code_info")
    except ValueError:
        pass
    try:
        variables.remove("input_file")
    except ValueError:
        pass

    for variable in variables:
        print(f"Checking {variable}... ", end="")
        assert np.all(file1_data[variable][:] == file2_data[variable][:])
        print("Pass")


def compare_two_input_files(file1, file2):
    """Compares two input files, ignoring comments"""
    with open(file1, "r") as f:
        file1_contents = f.read()
    with open(file2, "r") as f:
        file2_contents = f.read()

    def clean_file(text):
        text = HEADER_CREATED_RE.sub("Timestamp", text)
        text = HEADER_UUID_RE.sub("UUID", text)
        text = HEADER_FILENAME_RE.sub("Filename", text)
        return text

    assert clean_file(file1_contents) == clean_file(file2_contents)


def test_used_inputs(tmp_path):
    """Run a short test case, and then re-run it using the "used_inputs.in"
    file. Results should be identical.
    """

    input_filename = "test_input.in"
    copy_input_file(input_filename, tmp_path)
    os.chdir(tmp_path)
    run_gs2(get_gs2_path(), input_filename)

    used_inputs_filename = "test_input.used_inputs.in"
    run_gs2(get_gs2_path(), used_inputs_filename)

    compare_two_output_files(
        tmp_path / "test_input.out.nc", tmp_path / "test_input.used_inputs.out.nc"
    )
    compare_two_input_files(
        tmp_path / used_inputs_filename,
        tmp_path / "test_input.used_inputs.used_inputs.in",
    )
