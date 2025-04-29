"""@file psdi-data-conversion/tests/cli_test.py

Created 2025-01-15 by Bryan Gillis.

Tests of the command-line interface
"""

import logging
import os
import shlex
import sys
from unittest.mock import patch

import pytest

from psdi_data_conversion import constants as const
from psdi_data_conversion.converter import D_CONVERTER_ARGS, L_REGISTERED_CONVERTERS, get_registered_converter_class
from psdi_data_conversion.converters.atomsk import CONVERTER_ATO
from psdi_data_conversion.converters.c2x import CONVERTER_C2X
from psdi_data_conversion.converters.openbabel import (CONVERTER_OB, COORD_GEN_KEY, COORD_GEN_QUAL_KEY,
                                                       DEFAULT_COORD_GEN, DEFAULT_COORD_GEN_QUAL)
from psdi_data_conversion.database import (get_conversion_quality, get_converter_info, get_format_info,
                                           get_in_format_args, get_out_format_args, get_possible_conversions,
                                           get_possible_formats)
from psdi_data_conversion.main import FileConverterInputException, parse_args
from psdi_data_conversion.testing.conversion_test_specs import l_cla_test_specs
from psdi_data_conversion.testing.utils import run_test_conversion_with_cla, run_with_arg_string
from psdi_data_conversion.utils import regularize_name


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


@pytest.mark.parametrize("test_spec", l_cla_test_specs,
                         ids=lambda x: x.name)
def test_conversions(test_spec):
    """Run all conversion tests in the defined list of test specifications
    """
    run_test_conversion_with_cla(test_spec)


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
    assert args.name == regularize_name("Open Babel")
    args = get_parsed_args("--list 'Open Babel'")
    assert args.name == regularize_name("Open Babel")
    args = get_parsed_args("-l Atomsk")
    assert args.name == regularize_name("Atomsk")


def test_input_processing():
    """Unit tests to ensure that the CLI properly processes input arguments to determine values that are needed but
    weren't provided
    """

    # Check that different ways of specifying converter are all processed correctly
    converter_name = "Open Babel"
    args = get_parsed_args(f"file1.mmcif -t pdb -w {converter_name}")
    assert args.name == regularize_name(converter_name)
    args = get_parsed_args(f"file1.mmcif -t pdb -w '{converter_name}'")
    assert args.name == regularize_name(converter_name)

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
    for converter_rname in L_REGISTERED_CONVERTERS:
        converter_name = get_registered_converter_class(converter_rname).name
        assert converter_name in captured.out, converter_name

    # Check that no errors were produced
    assert not captured.err


def test_detail_converter(capsys):
    """Test the option to provide detail on a converter
    """

    # Test all converters are recognised, don't raise an error, and we get info on them
    for name in L_REGISTERED_CONVERTERS:

        converter_info = get_converter_info(name)
        converter_name = get_registered_converter_class(name).name

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

        l_allowed_in_formats, l_allowed_out_formats = get_possible_formats(name)
        for in_format in l_allowed_in_formats:
            output_allowed = "yes" if in_format in l_allowed_out_formats else "no"
            assert string_is_present_in_out(f"{in_format.disambiguated_name}yes{output_allowed}")
        for out_format in l_allowed_out_formats:
            input_allowed = "yes" if out_format in l_allowed_in_formats else "no"
            assert string_is_present_in_out(f"{out_format.disambiguated_name}{input_allowed}yes")

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


def test_get_conversions(capsys):
    """Test the option to get information on converters which can perform a desired conversion
    """
    in_format = "xyz-0"
    out_format = "inchi"
    l_conversions = get_possible_conversions(in_format, out_format)

    run_with_arg_string(f"-l -f {in_format} -t {out_format}")
    captured = capsys.readouterr()
    compressed_out: str = captured.out.replace("\n", "").replace(" ", "")

    def string_is_present_in_out(s: str) -> bool:
        return s.replace("\n", " ").replace(" ", "") in compressed_out

    assert not captured.err

    assert bool(l_conversions) == string_is_present_in_out("The following registered converters can convert from "
                                                           f"{in_format} to {out_format}:")

    for name, _, _ in l_conversions:
        if name in L_REGISTERED_CONVERTERS:
            assert string_is_present_in_out(get_registered_converter_class(name).name)
    for name in L_REGISTERED_CONVERTERS:
        if name not in [x[0] for x in l_conversions]:
            assert not string_is_present_in_out(get_registered_converter_class(name).name)


def test_conversion_info(capsys):
    """Test the option to provide detail on degree of success and arguments a converter allows for a given conversion
    """

    converter_name = CONVERTER_OB

    in_format = "xyz"
    out_format = "inchi"
    qual = get_conversion_quality(converter_name, in_format, out_format)

    # Test a basic listing of arguments, checking with the converter name in lowercase to be sure that works
    run_with_arg_string(f"-l {converter_name.lower()} -f {in_format} -t {out_format}")
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

    # Now try listing for converters which don't yet allow in/out args

    in_format = "pdb-0"
    out_format = "cif"
    for converter_name in [CONVERTER_C2X, CONVERTER_ATO]:
        qual = get_conversion_quality(converter_name, in_format, out_format)

        run_with_arg_string(f"-l {converter_name} -f {in_format} -t {out_format}")
        captured = capsys.readouterr()
        compressed_out: str = captured.out.replace("\n", "").replace(" ", "")

        assert not captured.err


def test_format_info(capsys):
    """Test that we can get information on formats
    """

    # Try to get info on an unambiguous format

    in_format = "cif"
    in_format_info = get_format_info(in_format)
    run_with_arg_string(f"-l -f {in_format}")

    captured = capsys.readouterr()
    compressed_out: str = captured.out.replace("\n", "").replace(" ", "")

    def string_is_present_in_out(s: str) -> bool:
        return s.replace("\n", " ").replace(" ", "") in compressed_out

    assert not captured.err

    assert string_is_present_in_out(f"{in_format_info.id}: {in_format_info.name} "
                                    f"({in_format_info.note})")

    # Try to get info on an ambiguous format

    out_format = "pdb"
    l_out_format_info = get_format_info(out_format, which="all")
    run_with_arg_string(f"-l -t {out_format}")

    captured = capsys.readouterr()
    compressed_out: str = captured.out.replace("\n", "").replace(" ", "")

    assert not captured.err

    assert string_is_present_in_out(f"WARNING: Format '{out_format}' is ambiguous")

    for out_format_info in l_out_format_info:
        assert string_is_present_in_out(f"{out_format_info.id}: {out_format_info.disambiguated_name} "
                                        f"({out_format_info.note})")

    # Test we get expected errors for unrecognised formats

    in_format = 99999
    with pytest.raises(SystemExit):
        run_with_arg_string(f"-l -f {in_format}")

    captured = capsys.readouterr()
    compressed_err: str = captured.err.replace("\n", "").replace(" ", "")

    def string_is_present_in_err(s: str) -> bool:
        return s.replace("\n", " ").replace(" ", "") in compressed_err

    assert string_is_present_in_err(f"ERROR: Format '{in_format}' not recognised")

    out_format = "not_a_format"

    with pytest.raises(SystemExit):
        run_with_arg_string(f"-l -t {out_format}")

    captured = capsys.readouterr()
    compressed_err: str = captured.err.replace("\n", "").replace(" ", "")

    assert string_is_present_in_err(f"ERROR: Format '{out_format}' not recognised")
