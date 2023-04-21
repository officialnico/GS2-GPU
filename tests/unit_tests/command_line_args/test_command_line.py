#!/usr/bin/env python3
# This test should be run with `pytest`

import os
import pathlib
import re
import shutil
import subprocess


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


def test_version(tmp_path):
    """Check the version flag prints a version string we expect"""
    os.chdir(tmp_path)

    out = subprocess.run(
        [get_gs2_path(), "--version"], capture_output=True, text=True, check=True
    )
    version_re = re.compile(
        r"GS2 version \d+\.\d+\.\d+(-\d+)?(-g[a-f0-9]+)?(-dirty)?\n"
    )

    assert version_re.match(out.stdout) is not None


def test_help(tmp_path):
    """Check the help flag lists some flags we expect"""
    os.chdir(tmp_path)

    out = subprocess.run(
        [get_gs2_path(), "--help"], capture_output=True, text=True, check=True
    )

    assert "--help" in out.stdout
    assert "--version" in out.stdout
