"""@file psdi-data-conversion/tests/cli_test.py

Created 2025-01-15 by Bryan Gillis.

Tests of the command-line interface
"""

import os
import pytest
import shlex
import sys
from unittest.mock import patch

from psdi_data_conversion.main import FileConverterException, parse_args


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
                           r"--flags '\-ab \-c \--example'")
    assert args.l_args[0] == "file1"
    assert args.l_args[1] == "file2"
    assert args.input_dir == cwd
    assert args.to_format == "pdb"
    assert args.output_dir == f"{cwd}/.."
    assert args.converter == "Atomsk"
    assert args.flags == "-ab -c --example"

    # It should fail with no arguments
    with pytest.raises(FileConverterException):
        get_parsed_args("")

    # It should fail if the output format isn't specified
    with pytest.raises(FileConverterException):
        get_parsed_args("file1.mmcif")

    # It should fail if the input directory doesn't exist
    with pytest.raises(FileConverterException):
        get_parsed_args("file1.mmcif -i /no/where -t pdb")

    # Expect that it should work if we just ask for a list
    args = get_parsed_args("--list")
    assert args.list


def test_input_processing():
    """Unit tests to ensure that the CLI properly processes input arguments to determine values that are needed but
    weren't provided
    """
    # TODO
    pass
