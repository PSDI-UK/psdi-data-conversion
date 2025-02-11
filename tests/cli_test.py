"""@file psdi-data-conversion/tests/cli_test.py

Created 2025-01-15 by Bryan Gillis.

Tests of the command-line interface
"""

import os
import pytest
import shlex
import sys
from unittest.mock import patch

from psdi_data_conversion import constants as const
from psdi_data_conversion.converter import L_REGISTERED_CONVERTERS
from psdi_data_conversion.converters.atomsk import CONVERTER_ATO
from psdi_data_conversion.converters.openbabel import CONVERTER_OB
from psdi_data_conversion.database import (get_conversion_quality, get_converter_info, get_in_format_args,
                                           get_out_format_args, get_possible_converters, get_possible_formats)
from psdi_data_conversion.main import FileConverterInputException, main, parse_args


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
    args = get_parsed_args(f"file1 file2 -f mmcif -i {cwd} -t pdb -o {cwd}/.. -w 'Atomsk' " +
                           r"--delete-input --from-flags '\-ab \-c \--example' --to-flags '\-d' " +
                           "--coord-gen Gen3D best -q --log-file text.log")
    assert args.l_args[0] == "file1"
    assert args.l_args[1] == "file2"
    assert args.input_dir == cwd
    assert args.to_format == "pdb"
    assert args.output_dir == f"{cwd}/.."
    assert args.name == "Atomsk"
    assert args.delete_input is True
    assert args.from_flags == "-ab -c --example"
    assert args.to_flags == "-d"
    assert args.coord_gen == "Gen3D"
    assert args.coord_gen_qual == "best"
    assert args.quiet is True
    assert args.log_file == "text.log"
    assert args.log_mode == const.LOG_NONE

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
        get_parsed_args("file1.mmcif -t pdb --log-mode max")

    # It should work if we just ask for a list, and set log mode to stdout
    args = get_parsed_args("--list")
    assert args.list
    assert args.log_mode == const.LOG_STDOUT

    # We should also be able to ask for info on a specific converter
    args = get_parsed_args("-l Open Babel")
    assert args.name == "Open Babel"
    args = get_parsed_args("--list 'Open Babel'")
    assert args.name == "Open Babel"
    args = get_parsed_args("-l Atomsk")
    assert args.name == "Atomsk"


def test_input_processing():
    """Unit tests to ensure that the CLI properly processes input arguments to determine values that are needed but
    weren't provided
    """

    # Check that different ways of specifying converter are all processed correctly
    converter_name = "Open Babel"
    args = get_parsed_args(f"file1.mmcif -t pdb -w {converter_name}")
    assert args.name == converter_name
    args = get_parsed_args(f"file1.mmcif -t pdb -w '{converter_name}'")
    assert args.name == converter_name

    # Check that input dir defaults to the current directory
    cwd = os.getcwd()
    assert args.input_dir == cwd

    # Check that output dir defaults to match input dir
    output_check_args = get_parsed_args(f"file1.mmcif -i {cwd}/.. -t pdb")
    assert output_check_args.output_dir == f"{cwd}/.."

    # Check that we get the default coordinate generation options
    assert args.coord_gen == const.DEFAULT_COORD_GEN
    assert args.coord_gen_qual == const.DEFAULT_COORD_GEN_QUAL
    assert get_parsed_args("file1.mmcif -t pdb --coord-gen Gen3D").coord_gen_qual == const.DEFAULT_COORD_GEN_QUAL

    # Check that trying to get the log file raises an exception due to the test file not existing
    with pytest.raises(FileConverterInputException):
        assert args.log_file == "file1" + const.LOG_EXT

    # Check that the log file uses the expected default value in list mode
    list_check_args = get_parsed_args("--list")
    assert list_check_args.log_file == const.DEFAULT_LISTING_LOG_FILE


def test_list_converters(capsys):
    """Test the option to list available converters
    """
    run_with_arg_string("--list")
    captured = capsys.readouterr()
    assert "Available converters:" in captured.out
    for converter_name in L_REGISTERED_CONVERTERS:
        assert converter_name in captured.out, converter_name

    # Check that no errors were produced
    assert not captured.err


def test_detail_converter(capsys):
    """Test the option to provide detail on a converter
    """

    # Test all converters are recognised, don't raise an error, and we get info on them
    for converter_name in L_REGISTERED_CONVERTERS:

        converter_info = get_converter_info(converter_name)

        run_with_arg_string(f"--list {converter_name}")
        captured = capsys.readouterr()
        compressed_out: str = captured.out.replace("\n", "").replace(" ", "")

        def string_is_present_in_out(s: str) -> bool:
            return s.replace("\n", " ").replace(" ", "") in compressed_out

        assert string_is_present_in_out(converter_name)

        if not converter_info.description:
            assert "available for this converter" in captured.out
        else:
            assert string_is_present_in_out(converter_info.description)

        # Check for URL
        assert converter_info.url in captured.out

        # Check for list of allowed input/output formats
        assert "   INPUT  OUTPUT" in captured.out

        l_allowed_in_formats, l_allowed_out_formats = get_possible_formats(converter_name)
        for in_format in l_allowed_in_formats:
            output_allowed = "yes" if in_format in l_allowed_out_formats else "no"
            assert string_is_present_in_out(f"{in_format}yes{output_allowed}")
        for out_format in l_allowed_out_formats:
            input_allowed = "yes" if out_format in l_allowed_in_formats else "no"
            assert string_is_present_in_out(f"{out_format}{input_allowed}yes")

        # Check that no errors were produced
        assert not captured.err

    # Test we do get an error for a bad converter name
    with pytest.raises(SystemExit):
        run_with_arg_string("--list bad_converter")
    captured = capsys.readouterr()
    assert "not recognized" in captured.err

    # Test that we can also provide the converter name with -w/--with
    run_with_arg_string(f"-l -w {CONVERTER_ATO}")
    captured = capsys.readouterr()
    assert not captured.err
    assert CONVERTER_ATO in captured.out
    assert const.CONVERTER_DEFAULT not in captured.out


def test_get_converters(capsys):
    """Test the option to get information on converters which can perform a desired conversion
    """
    in_format = "xyz"
    out_format = "inchi"
    l_converters = get_possible_converters(in_format, out_format)

    run_with_arg_string(f"-l -f {in_format} -t {out_format}")
    captured = capsys.readouterr()
    compressed_out: str = captured.out.replace("\n", "").replace(" ", "")

    def string_is_present_in_out(s: str) -> bool:
        return s.replace("\n", " ").replace(" ", "") in compressed_out

    assert not captured.err

    assert bool(l_converters) == string_is_present_in_out("The following converters can convert from "
                                                          f"{in_format} to {out_format}:")

    for converter_name in l_converters:
        if converter_name in L_REGISTERED_CONVERTERS:
            assert string_is_present_in_out(converter_name)
    for converter_name in L_REGISTERED_CONVERTERS:
        if converter_name not in l_converters:
            assert not string_is_present_in_out(converter_name)


def test_conversion_info(capsys):
    """Test the option to provide detail on degree of success and arguments a converter allows for a given conversion
    """

    converter_name = CONVERTER_OB
    in_format = "xyz"
    out_format = "inchi"
    qual_str, details, d_prop_conv_info = get_conversion_quality(converter_name, in_format, out_format)

    # Test a basic listing of arguments
    run_with_arg_string(f"-l {converter_name} -f {in_format} -t {out_format}")
    captured = capsys.readouterr()
    compressed_out: str = captured.out.replace("\n", "").replace(" ", "")

    def string_is_present_in_out(s: str) -> bool:
        return s.replace("\n", " ").replace(" ", "") in compressed_out

    assert not captured.err

    # Check that conversion quality details are in the output as expected
    assert string_is_present_in_out(f"Conversion from '{in_format}' to '{out_format}' with {converter_name} is "
                                    f"possible with {qual_str} conversion quality")
    assert string_is_present_in_out("Notes on this conversion:")
    assert string_is_present_in_out(const.QUAL_NOTE_OUT_MISSING % const.QUAL_2D_LABEL)
    assert string_is_present_in_out(const.QUAL_NOTE_OUT_MISSING % const.QUAL_3D_LABEL)
    assert string_is_present_in_out(const.QUAL_NOTE_IN_MISSING % const.QUAL_CONN_LABEL)

    l_in_flags, l_in_options = get_in_format_args(converter_name, in_format)
    l_out_flags, l_out_options = get_out_format_args(converter_name, out_format)

    # Check headings for input/output flags/options are present if and only if some of those flags/options exist
    assert bool(l_in_flags) == string_is_present_in_out(f"Allowed input flags for format '{in_format}':")
    assert bool(l_out_flags) == string_is_present_in_out(f"Allowed output flags for format '{out_format}':")
    assert bool(l_in_options) == string_is_present_in_out(f"Allowed input options for format '{in_format}':")
    assert bool(l_out_options) == string_is_present_in_out(f"Allowed output options for format '{out_format}':")

    # Check that info for each flag and option is printed as expected
    for flag_info in l_in_flags + l_out_flags:
        info = flag_info.info if flag_info.info and flag_info.info != "N/A" else ""
        assert string_is_present_in_out(f"{flag_info.flag}{flag_info.description}{info}")
    for option_info in l_in_options + l_out_options:
        info = option_info.info if option_info.info and option_info.info != "N/A" else ""
        assert string_is_present_in_out(f"{option_info.flag}<{option_info.brief}>{option_info.description}{info}")


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
    basic_arg_string = f"{input_filename} -t {to_format} -i {input_dir} -o {output_dir}"
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
    run_with_arg_string(f"{test_filename_base} -f {from_format} -t {to_format} -i {input_dir} -o {output_dir}")
    captured = capsys.readouterr()
    assert "Success!" in captured.out
    assert "ERROR" not in captured.err


def test_archive_convert(tmp_path_factory, capsys, test_data_loc):
    """Test running conversion on archives of files
    """

    test_filename_base = "caffeine-smi"
    l_archive_formats = ["zip", "tar", "tar.gz"]

    l_ex_filename_bases = ["caffeine-no-flags",
                           "caffeine-ia",
                           "caffeine-ia-ox",
                           "caffeine-ia-okx",
                           "caffeine-ia-okx-oof4",
                           "caffeine-ia-okx-oof4l5",]
    to_format = "inchi"

    for archive_format in l_archive_formats:
        input_filename = f"{test_filename_base}.{archive_format}"

        input_dir = tmp_path_factory.mktemp("input")
        output_dir = tmp_path_factory.mktemp("output")

        # Symlink the input file from the test_data directory to the input directory
        os.symlink(os.path.join(test_data_loc, input_filename),
                   os.path.join(input_dir, input_filename))

        # Run the conversion
        basic_arg_string = f"{input_filename} -t {to_format} -i {input_dir} -o {output_dir}"
        run_with_arg_string(basic_arg_string)

        # Check that all expected output files exist
        for ex_filename_base in l_ex_filename_bases:
            ex_output_filename = os.path.join(output_dir, f"{ex_filename_base}.{to_format}")
            assert os.path.isfile(ex_output_filename)
