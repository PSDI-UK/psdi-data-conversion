"""@file psdi-data-conversion/tests/cli_test.py

Created 2025-01-15 by Bryan Gillis.

Tests of the command-line interface
"""

import logging
import os
import pytest
import shlex
import sys
from unittest.mock import patch

from psdi_data_conversion import constants as const
from psdi_data_conversion.converter import D_CONVERTER_ARGS, L_REGISTERED_CONVERTERS
from psdi_data_conversion.converters.c2x import CONVERTER_C2X
from psdi_data_conversion.converters.openbabel import (CONVERTER_OB, COORD_GEN_KEY, COORD_GEN_QUAL_KEY,
                                                       DEFAULT_COORD_GEN, DEFAULT_COORD_GEN_QUAL)
from psdi_data_conversion.database import (get_conversion_quality, get_converter_info, get_in_format_args,
                                           get_out_format_args, get_possible_converters, get_possible_formats)
from psdi_data_conversion.main import FileConverterInputException, parse_args
from psdi_data_conversion.testing.utils import run_test_conversion_with_cla, run_with_arg_string
from psdi_data_conversion.testing import conversion_test_specs as specs


def test_unique_args():
    """Check that all converter-specific arguments have unique names
    """
    s_arg_names = set()
    for name in L_REGISTERED_CONVERTERS:
        for arg_name, _, _ in D_CONVERTER_ARGS[name]:
            assert arg_name not in s_arg_names, ("Name clash between converters, with multiple using the argument "
                                                 f"'{arg_name}'")
            s_arg_names.add(arg_name)


def get_parsed_args(s):
    """Performs argument parsing on a string which represents what the arguments would be after the function call
    """
    l_args = shlex.split("test " + s)
    with patch.object(sys, 'argv', l_args):
        return parse_args()


@pytest.fixture(autouse=True)
def setup_test() -> None:
    """Reset global aspects before a test, so that different tests won't interfere with each other"""

    # Remove the global log file if one exists
    try:
        os.remove(const.GLOBAL_LOG_FILENAME)
    except FileNotFoundError:
        pass

    # Clear any existing loggers so new ones will be created fresh
    logging.Logger.manager.loggerDict.clear()


def test_default():
    """Test that the default converter is registered.
    """
    assert const.CONVERTER_DEFAULT in L_REGISTERED_CONVERTERS


def test_basic_conversions():
    """Run a basic set of conversions with various converters and file formats which we expect to succeed without
    issue
    """
    run_test_conversion_with_cla(specs.basic_tests)


def test_archive_convert():
    """Run a test of converting an archive of files
    """
    run_test_conversion_with_cla(specs.archive_tests)


def test_archive_wrong_format():
    """Run a test that converting an archive but specifying the wrong input format will produce a warning
    """
    run_test_conversion_with_cla(specs.archive_wrong_format_test)


def test_log_mode():
    """Test that the various log modes result in the expected log files being created
    """
    run_test_conversion_with_cla(specs.log_mode_tests)


def test_stdout():
    """Test that the output is sent to stdout when requested
    """
    run_test_conversion_with_cla(specs.stdout_test)


def test_quiet():
    """Test that quiet mode suppresses all output
    """
    run_test_conversion_with_cla(specs.quiet_test)


def test_open_babel_warning():
    """Run a test that expected warnings from Open Babel are captured in the log
    """
    run_test_conversion_with_cla(specs.open_babel_warning_test)


def test_invalid_converter():
    """Run a test of the converter to ensure it reports an error properly if an invalid converter is requested
    """
    run_test_conversion_with_cla(specs.invalid_converter_test)


def test_quality_note():
    """Run a test of the converter on an `.xyz` to `.inchi` conversion which we expect to have warnings about data
    loss and extrapolation
    """
    run_test_conversion_with_cla(specs.quality_note_test)


def test_cleanup():
    """Test that input files are deleted if requested
    """
    run_test_conversion_with_cla(specs.cleanup_input_test)


def test_failed_conversion():
    """Run a test of the converter on a conversion we expect to fail
    """
    run_test_conversion_with_cla(specs.failed_conversion_test)


def test_format_args():
    """Run a test that format args are processed correctly
    """
    run_test_conversion_with_cla(specs.format_args_test)


def test_coord_gen():
    """Run a test that coordinate generation args are processed correctly
    """
    run_test_conversion_with_cla(specs.coord_gen_test)


def test_input_validity():
    """Unit tests to ensure that the CLI properly checks for valid input
    """

    # Test that we get what we put in for a standard execution
    cwd = os.getcwd()
    args = get_parsed_args(f"file1 file2 -f mmcif -i {cwd} -t pdb -o {cwd}/.. -w '{CONVERTER_C2X}' " +
                           r"--delete-input --from-flags '\-ab \-c \--example' --to-flags '\-d' " +
                           r"--from-options '-x xval --xopt xoptval' --to-options '-y yval --yopt yoptval' "
                           "--strict --nc --coord-gen Gen3D best -q --log-file text.log")
    assert args.l_args[0] == "file1"
    assert args.l_args[1] == "file2"
    assert args.input_dir == cwd
    assert args.to_format == "pdb"
    assert args.output_dir == f"{cwd}/.."
    assert args.name == CONVERTER_C2X
    assert args.no_check is True
    assert args.strict is True
    assert args.delete_input is True
    assert args.from_flags == "-ab -c --example"
    assert args.to_flags == "-d"
    assert args.from_options == "-x xval --xopt xoptval"
    assert args.to_options == "-y yval --yopt yoptval"
    assert args.quiet is True
    assert args.log_file == "text.log"
    assert args.log_mode == const.LOG_NONE

    # Test Open-Babel-specific arguments
    args = get_parsed_args(f"file1 -t pdb -w '{CONVERTER_OB}' --coord-gen Gen3D best")
    assert args.d_converter_args[COORD_GEN_KEY] == "Gen3D"
    assert args.d_converter_args[COORD_GEN_QUAL_KEY] == "best"

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
    assert args.d_converter_args[COORD_GEN_KEY] == DEFAULT_COORD_GEN
    assert args.d_converter_args[COORD_GEN_QUAL_KEY] == DEFAULT_COORD_GEN_QUAL
    assert (get_parsed_args("file1.mmcif -t pdb --coord-gen Gen3D").d_converter_args[COORD_GEN_QUAL_KEY] ==
            DEFAULT_COORD_GEN_QUAL)

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

    # Test we do get a simple error for a bad converter name
    with pytest.raises(SystemExit):
        run_with_arg_string("--list bad_converter")
    captured = capsys.readouterr()
    assert "not recognized" in captured.err
    assert "Traceback" not in captured.out
    assert "Traceback" not in captured.err

    # Test that we can also provide the converter name with -w/--with
    run_with_arg_string(f"-l -w {CONVERTER_C2X}")
    captured = capsys.readouterr()
    assert not captured.err
    assert CONVERTER_C2X in captured.out
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

    assert bool(l_converters) == string_is_present_in_out("The following registered converters can convert from "
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
    qual = get_conversion_quality(converter_name, in_format, out_format)

    # Test a basic listing of arguments
    run_with_arg_string(f"-l {converter_name} -f {in_format} -t {out_format}")
    captured = capsys.readouterr()
    compressed_out: str = captured.out.replace("\n", "").replace(" ", "")

    def string_is_present_in_out(s: str) -> bool:
        return s.replace("\n", " ").replace(" ", "") in compressed_out

    assert not captured.err

    # Check that conversion quality details are in the output as expected
    assert string_is_present_in_out(f"Conversion from '{in_format}' to '{out_format}' with {converter_name} is "
                                    f"possible with {qual.qual_str} conversion quality")
    assert string_is_present_in_out("WARNING: Potential data loss or extrapolation issues with this conversion:")
    assert string_is_present_in_out(const.QUAL_NOTE_OUT_MISSING.format(const.QUAL_2D_LABEL))
    assert string_is_present_in_out(const.QUAL_NOTE_OUT_MISSING.format(const.QUAL_3D_LABEL))
    assert string_is_present_in_out(const.QUAL_NOTE_IN_MISSING.format(const.QUAL_CONN_LABEL))

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
