"""@file psdi-data-conversion/tests/cli_test.py

Created 2025-01-15 by Bryan Gillis.

Tests of the command-line interface
"""

import os
import pytest
import shlex
import sys
from unittest.mock import patch

from psdi_data_conversion.main import (DEFAULT_COORD_GEN, DEFAULT_COORD_GEN_QUAL, DEFAULT_LISTING_LOG_FILE,
                                       FileConverterInputException, LOG_EXT, parse_args)


def get_parsed_args(s):
    """Performs argument parsing on a string which represents what the arguments would be after the function call
    """
    l_args = shlex.split("test " + s)
    with patch.object(sys, 'argv', l_args):
        return parse_args()


def test_input_validity():
    """Unit tests to ensure that the CLI properly checks for valid input
    """

    # Test that we get what we put in for a standard execution
    cwd = os.getcwd()
    args = get_parsed_args(f"file1 file2 -f mmcif -i {cwd} -t pdb -a {cwd}/.. -w 'Atomsk' "
                           r"--from-flags '\-ab \-c \--example' --to-flags '\-d' " +
                           "--coord-gen Gen3D best")
    assert args.l_args[0] == "file1"
    assert args.l_args[1] == "file2"
    assert args.input_dir == cwd
    assert args.to_format == "pdb"
    assert args.output_dir == f"{cwd}/.."
    assert args.converter == "Atomsk"
    assert args.from_flags == "-ab -c --example"
    assert args.to_flags == "-d"
    assert args.coord_gen == "Gen3D"
    assert args.coord_gen_qual == "best"

    # It should fail with no arguments
    with pytest.raises(FileConverterInputException):
        get_parsed_args("")

    # It should fail if the output format isn't specified
    with pytest.raises(FileConverterInputException):
        get_parsed_args("file1.mmcif")

    # It should fail if the input directory doesn't exist
    with pytest.raises(FileConverterInputException):
        get_parsed_args("file1.mmcif -i /no/where -t pdb")

    # It should fail with bad or too many arguments to --coord-gen
    with pytest.raises(FileConverterInputException):
        get_parsed_args("file1.mmcif -t pdb --coord-gen Gen1D")
    with pytest.raises(FileConverterInputException):
        get_parsed_args("file1.mmcif -t pdb --coord-gen Gen3D worst")
    with pytest.raises(FileConverterInputException):
        get_parsed_args("file1.mmcif -t pdb --coord-gen Gen3D best quality")

    # It should work if we just ask for a list
    args = get_parsed_args("--list")
    assert args.list


def test_input_processing():
    """Unit tests to ensure that the CLI properly processes input arguments to determine values that are needed but
    weren't provided
    """

    # Check that input format is determined from the first file in a list
    args = get_parsed_args("file1.mmcif file2.pdb -t pdb")
    assert args.from_format == "mmcif"

    # Check that input dir defaults to the current directory
    cwd = os.getcwd()
    assert args.input_dir == cwd

    # Check that output dir defaults to match input dir
    output_check_args = get_parsed_args(f"file1.mmcif -i {cwd}/.. -t pdb")
    assert output_check_args.output_dir == f"{cwd}/.."

    # Check that we get the default coordinate generation options
    assert args.coord_gen == DEFAULT_COORD_GEN
    assert args.coord_gen_qual == DEFAULT_COORD_GEN_QUAL
    assert get_parsed_args("file1.mmcif -t pdb --coord-gen Gen3D").coord_gen_qual == DEFAULT_COORD_GEN_QUAL

    # Check that log file is based off of the first file name in normal mode
    assert args.log_file == "file1" + LOG_EXT

    # Check that the log file uses the expected default value in list mode
    list_check_args = get_parsed_args("--list")
    assert list_check_args.log_file == DEFAULT_LISTING_LOG_FILE
