"""@file psdi-data-conversion/tests/cli_test.py

Created 2025-01-15 by Bryan Gillis.

Tests of the command-line interface
"""

import os
import pytest
import shlex
import sys
from unittest.mock import patch

from psdi_data_conversion.converter import L_ALLOWED_CONVERTERS, LOGGING_NONE
from psdi_data_conversion.main import (DEFAULT_COORD_GEN, DEFAULT_COORD_GEN_QUAL, DEFAULT_LISTING_LOG_FILE,
                                       FileConverterInputException, main, LOG_EXT, parse_args)


def get_parsed_args(s):
    """Performs argument parsing on a string which represents what the arguments would be after the function call
    """
    l_args = shlex.split("test " + s)
    with patch.object(sys, 'argv', l_args):
        return parse_args()


def run_with_arg_string(s):
    """Runs the convert script with the provided argument string
    """
    l_args = shlex.split("test " + s)
    with patch.object(sys, 'argv', l_args):
        main()


def test_input_validity():
    """Unit tests to ensure that the CLI properly checks for valid input
    """

    # Test that we get what we put in for a standard execution
    cwd = os.getcwd()
    args = get_parsed_args(f"file1 file2 -f mmcif -i {cwd} -t pdb -a {cwd}/.. -w 'Atomsk' -d " +
                           r"--from-flags '\-ab \-c \--example' --to-flags '\-d' " +
                           "--coord-gen Gen3D best -q --log-file text.log")
    assert args.l_args[0] == "file1"
    assert args.l_args[1] == "file2"
    assert args.input_dir == cwd
    assert args.to_format == "pdb"
    assert args.output_dir == f"{cwd}/.."
    assert args.converter == "Atomsk"
    assert args.delete_input is True
    assert args.from_flags == "-ab -c --example"
    assert args.to_flags == "-d"
    assert args.coord_gen == "Gen3D"
    assert args.coord_gen_qual == "best"
    assert args.quiet is True
    assert args.log_file == "text.log"
    assert args.logging_mode == LOGGING_NONE

    # It should fail with no arguments
    with pytest.raises(FileConverterInputException):
        get_parsed_args("")

    # It should fail if the output format isn't specified
    with pytest.raises(FileConverterInputException):
        get_parsed_args("file1.mmcif")

    # It should fail if the input directory doesn't exist
    with pytest.raises(FileConverterInputException):
        get_parsed_args("file1.mmcif -i /no/where -t pdb")

    # It should fail if the converter isn't recognized
    with pytest.raises(FileConverterInputException):
        get_parsed_args("file1.mmcif -t pdb -w Ato")

    # It should fail with bad or too many arguments to --coord-gen
    with pytest.raises(FileConverterInputException):
        get_parsed_args("file1.mmcif -t pdb --coord-gen Gen1D")
    with pytest.raises(FileConverterInputException):
        get_parsed_args("file1.mmcif -t pdb --coord-gen Gen3D worst")
    with pytest.raises(FileConverterInputException):
        get_parsed_args("file1.mmcif -t pdb --coord-gen Gen3D best quality")

    # It should fail if it doesn't recognise the logging mode
    with pytest.raises(FileConverterInputException):
        get_parsed_args("file1.mmcif -t pdb --logging-mode max")

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

    # Check that trying to get the log file raises an exception due to the test file not existing
    with pytest.raises(FileConverterInputException):
        assert args.log_file == "file1" + LOG_EXT

    # Check that the log file uses the expected default value in list mode
    list_check_args = get_parsed_args("--list")
    assert list_check_args.log_file == DEFAULT_LISTING_LOG_FILE


def test_list_converters(capsys):
    """Test the option to list available converters
    """
    run_with_arg_string("--list")
    captured = capsys.readouterr()
    assert "Available converters are:" in captured.out
    for converter_name in L_ALLOWED_CONVERTERS:
        assert converter_name in captured.out


def test_detail_converter(capsys):
    """Test the option to provide detail on a converter
    """

    # Test all converters are recognised and don't raise an error
    for converter_name in L_ALLOWED_CONVERTERS:
        run_with_arg_string(f"--list {converter_name}")
        captured = capsys.readouterr()
        assert "not recognized" not in captured.err
        assert "Converter use detailing is still TBD" in captured.out

    # Test we do get an error for a bad converter name
    run_with_arg_string("--list bad_converter")
    captured = capsys.readouterr()
    assert "not recognized" in captured.err


def test_convert(tmp_path_factory, capsys, test_data_loc):
    """Test running file conversions
    """
    input_dir = tmp_path_factory.mktemp("input")
    output_dir = tmp_path_factory.mktemp("output")

    test_filename_base = "1NE6"

    from_format = "mmcif"
    to_format = "pdb"

    input_filename = f"{test_filename_base}.{from_format}"
    output_filename = f"{test_filename_base}.{to_format}"

    # Symlink the input file from the test_data directory to the input directory
    os.symlink(os.path.join(test_data_loc, input_filename),
               os.path.join(input_dir, input_filename))

    # Run a basic conversion
    basic_arg_string = f"{input_filename} -t {to_format} -i {input_dir} -a {output_dir}"
    run_with_arg_string(basic_arg_string)

    # Check that the expected output file has been created
    ex_output_file = os.path.join(output_dir, f"{output_filename}")
    assert os.path.isfile(ex_output_file), f"Expected output file {ex_output_file} does not exist"

    # Check logs and output
    captured = capsys.readouterr()
    assert "Success!" in captured.out
    assert "ERROR" not in captured.err

    # Check that running in quiet mode suppresses output
    run_with_arg_string(basic_arg_string + " -q")
    captured = capsys.readouterr()
    assert "Success!" not in captured.out
    assert "ERROR" not in captured.err

    # Test a call we expect to fail due to invalid input type being provided
    run_with_arg_string(basic_arg_string + " -f pdb")
    captured = capsys.readouterr()
    assert "Success!" not in captured.out
    assert "ERROR" in captured.err

    # Check that we can specify a file with its format instead of extension
    run_with_arg_string(f"{test_filename_base} -f {from_format} -t {to_format} -i {input_dir} -a {output_dir}")
    captured = capsys.readouterr()
    assert "Success!" in captured.out
    assert "ERROR" not in captured.err
