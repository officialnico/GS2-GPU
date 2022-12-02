#!/usr/bin/env python3

import os
import pathlib
import subprocess

from config_generation import constructFromFile
from get_all_matching_files import get_all_matching_files

import pytest


def get_gs2_top_level() -> pathlib.Path:
    """Returns the absolute path to the top-level GS2 directory

    Can be controlled by setting the GK_HEAD_DIR environment variable
    """
    default_path = pathlib.Path(__file__).parent / "../../.."
    gs2_path = pathlib.Path(os.environ.get("GK_HEAD_DIR", default_path))
    return gs2_path.absolute()


def get_gs2_bin_path() -> pathlib.Path:
    """Returns the absolute path to the GS2 binary directory

    Can be controlled by setting the GS2_BIN_PATH environment variable
    """
    default_path = get_gs2_top_level() / "bin"
    gs2_path = pathlib.Path(os.environ.get("GS2_BIN_PATH", default_path))
    return gs2_path.absolute()


@pytest.fixture(scope="module")
def default_input_file(tmp_path_factory) -> str:
    tmp_dir = tmp_path_factory.mktemp("namelist_collection")
    os.chdir(tmp_dir)
    write_default_input_file = (
        get_gs2_bin_path() / "write_default_input_file"
    ).absolute()
    subprocess.run([write_default_input_file], check=True)
    with open(tmp_dir / "defaults.in") as f:
        return f.read()


def get_build_config():
    gs2_bin = (get_gs2_bin_path() / "gs2").absolute()
    output = subprocess.run(
        [gs2_bin, "--build-config"], check=True, text=True, capture_output=True
    ).stdout.splitlines()

    return {"eigensolver": "Eigensolver: T" in output}


def test_all_namelists_present(default_input_file):

    top_dir = get_gs2_top_level()

    files = get_all_matching_files(
        paths=[top_dir / "src/", top_dir / "externals/utils/"],
        patterns=["*.[fF]??"],
        exclude=[top_dir / "src/templates"],
    )

    files = list(set(files))

    namelists = set()

    for file in files:
        if (tmp := constructFromFile(file, quiet=True)) is None:
            continue
        for v in tmp.values():
            namelists.add(v.namelistName)

    build_config = get_build_config()
    if not build_config["eigensolver"]:
        namelists.discard("eigval_knobs")

    missing = [
        namelist for namelist in namelists if (namelist not in default_input_file)
    ]

    assert missing == [], f"Namelists missing from default input file: {missing}"
