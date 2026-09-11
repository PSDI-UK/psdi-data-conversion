"""@file psdi-data-conversion/tests/cli_test.py

Created 2025-01-15 by Bryan Gillis.

Tests of the command-line interface
"""

import logging
import os
import shlex
import sys
from tempfile import TemporaryDirectory
from unittest.mock import patch

import pytest

from psdi_data_conversion import constants as const
from psdi_data_conversion.converter import D_CONVERTER_ARGS, L_REGISTERED_CONVERTERS, get_registered_converter_class
from psdi_data_conversion.converters.openbabel.converter import (COORD_GEN_KEY, COORD_GEN_QUAL_KEY, DEFAULT_COORD_GEN,
                                                                 DEFAULT_COORD_GEN_QUAL)
from psdi_data_conversion.database import (D_FORMAT_PROPERTY_ATTRS, get_conversion_pathway, get_conversion_quality,
                                           get_converter_info, get_format_info, get_in_format_args,
                                           get_out_format_args, get_possible_conversions, get_possible_formats)
from psdi_data_conversion.main import FileConverterInputException, parse_args
from psdi_data_conversion.testing.constants import FORMAT_INCHI, FORMAT_MOLDY
from psdi_data_conversion.testing.conversion_test_specs import l_cla_test_specs
from psdi_data_conversion.testing.utils import run_test_conversion_with_cla, run_with_arg_string
from psdi_data_conversion.utils import regularize_name, strip_control_codes


def _compress_text(s: str):
    """Strips whitespace and control codes from output to ease comparisons without worrying about things like line-
    wrapping"""
    return strip_control_codes(s.replace("\n", "").replace(" ", ""))


def _compressed_match(s1, s2):
    """Assert that s1 is contained in s2, ignoring control codes and whitespace"""
    s1_compressed = _compress_text(str(s1))
    s2_compressed = _compress_text(str(s2))
    return s1_compressed in s2_compressed


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
def setup_test():
    """Reset global aspects before a test, so that different tests won't interfere with each other"""

    # Remove the global log file if one exists
    try:
        os.remove(const.GLOBAL_LOG_FILENAME)
    except FileNotFoundError:
        pass

    # Clear any existing loggers so new ones will be created fresh
    logging.Logger.manager.loggerDict.clear()

    # Change directory to a temporary directory, so we can be sure that the script can be run from anywhere and not
    # just the project directory and/or $HOME
    old_cwd = os.getcwd()
    with TemporaryDirectory(prefix="test_cwd") as tmp_cwd:
        os.chdir(tmp_cwd)
        yield
    os.chdir(old_cwd)


@pytest.mark.parametrize("test_spec", l_cla_test_specs,
                         ids=lambda x: x.name)
def test_conversions(test_spec):
    """Run all conversion tests in the defined list of test specifications
    """
    run_test_conversion_with_cla(test_spec)


def test_general_arg_parsing():
    """Test that a standard argument string is parsed properly
    """
    cwd = os.getcwd()
    args = get_parsed_args(f"file1 file2 -f mmcif -i {cwd} -t pdb -o {cwd}/.. -w '{const.CONVERTER_C2X}' " +
                           r"--delete-input --from-flags '\-ab \-c \--example' --to-flags '\-d' " +
                           r"--from-options '-x xval --xopt xoptval' --to-options '-y yval --yopt yoptval' "
                           "--strict --nc --coord-gen Gen3D best -q --log-file text.log")
    assert args.l_args[0] == "file1"
    assert args.l_args[1] == "file2"
    assert args.input_dir == cwd
    assert args.to_format == "pdb"
    assert args.output_dir == f"{cwd}/.."
    assert args.name == const.CONVERTER_C2X
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


def test_open_babel_args():
    """Test that Open-Babel-specific arguments are parsed correctly
    """
    args = get_parsed_args(f"file1.mmcif -t pdb -w '{const.CONVERTER_OB}' --coord-gen Gen3D best")
    assert args.d_converter_args[COORD_GEN_KEY] == "Gen3D"
    assert args.d_converter_args[COORD_GEN_QUAL_KEY] == "best"


def test_fail_no_args():
    """Test that the parsing fails if no arguments are provided"""
    with pytest.raises(FileConverterInputException) as e:
        get_parsed_args("")
    assert _compressed_match("One or more names of files to convert must be provided", e.value)


def test_fail_no_to_format():
    """Test that the parsing fails if output format isn't specified"""
    with pytest.raises(FileConverterInputException) as e:
        get_parsed_args("file1.mmcif")
    assert _compressed_match("Output format (`-t/--to`) must be provided", e.value)


def test_fail_no_input_dir():
    """Test that the parsing fails if the input directory doesn't exist"""
    with pytest.raises(FileConverterInputException) as e:
        get_parsed_args("file1.mmcif -i /no/where -t pdb")
    assert _compressed_match("The provided input directory '/no/where' does not exist as a directory", e.value)


def test_fail_invalid_converter():
    """Test that the parsing fails if the converter isn't recognized"""
    with pytest.raises(FileConverterInputException) as e:
        get_parsed_args("file1.mmcif -t pdb -w FakeConverter")
    assert _compressed_match("Converter 'fakeconverter' not recognised", e.value)


def test_fail_bad_coord_gen_type():
    """Test that the parsing fails with bad --coord-gen type"""
    with pytest.raises(FileConverterInputException) as e:
        get_parsed_args(f"file1.mmcif -t pdb -w '{const.CONVERTER_OB}' --coord-gen Gen1D")
    assert _compressed_match("Coordinate generation type 'Gen1D' not recognised.", e.value)


def test_fail_bad_coord_gen_quality():
    """Test that the parsing fails with bad --coord-gen quality"""
    with pytest.raises(FileConverterInputException) as e:
        get_parsed_args(f"file1.mmcif -t pdb -w '{const.CONVERTER_OB}' --coord-gen Gen3D worst")
    assert _compressed_match("Coordinate generation quality 'worst' not recognised.", e.value)


def test_fail_bad_coord_gen_len():
    """Test that the parsing fails too many args to --coord-gen"""
    with pytest.raises(FileConverterInputException) as e:
        get_parsed_args(f"file1.mmcif -t pdb -w '{const.CONVERTER_OB}' --coord-gen Gen3D best quality")
    assert _compressed_match("At most two arguments may be provided to `--coord-gen`", e.value)


def test_fail_bad_logging_mode():
    """Test that the parsing fails if it doesn't recognise the logging mode"""
    with pytest.raises(FileConverterInputException) as e:
        get_parsed_args(f"file1.mmcif -t pdb -w '{const.CONVERTER_OB}' --log-mode max")
    assert _compressed_match("Unrecognised logging mode: 'max'", e.value)


def test_list_args():
    """Test that the parsing works if we just ask for a list, and set log mode to stdout"""
    args = get_parsed_args("--list")
    assert args.list
    assert args.log_mode == const.LOG_STDOUT


def test_list_converter():
    """Test that the parsing works if we ask for info on a specific converter"""
    args = get_parsed_args("-l Open Babel")
    assert args.name == regularize_name("Open Babel")
    args = get_parsed_args("--list 'Open Babel'")
    assert args.name == regularize_name("Open Babel")
    args = get_parsed_args("-l Atomsk")
    assert args.name == regularize_name("Atomsk")


def test_converter_input():
    """Test that the converter specified with -w/--with is properly parsed
    """
    args = get_parsed_args(f"file1.mmcif -t pdb -w {const.CONVERTER_OB}")
    assert args.name == regularize_name(const.CONVERTER_OB)
    args = get_parsed_args(f"file1.mmcif -t pdb -w '{const.CONVERTER_OB}'")
    assert args.name == regularize_name(const.CONVERTER_OB)


def test_default_input_dir():
    """Test that the input directory is set to the current directory if not specified"""
    args = get_parsed_args(f"file1.mmcif -t pdb -w {const.CONVERTER_OB}")
    assert args.input_dir == os.getcwd()


def test_default_output_dir():
    """Test that the output dir defaults to match input dir"""
    cwd = os.getcwd()
    args = get_parsed_args(f"file1.mmcif -i {cwd}/.. -t pdb -w {const.CONVERTER_OB}")
    assert args.output_dir == f"{cwd}/.."


def test_default_coord_gen():
    """Test that we get the default coordinate generation options if they aren't explicitly specified"""
    args = get_parsed_args(f"file1.mmcif -t pdb -w {const.CONVERTER_OB}")
    assert args.d_converter_args[COORD_GEN_KEY] == DEFAULT_COORD_GEN
    assert args.d_converter_args[COORD_GEN_QUAL_KEY] == DEFAULT_COORD_GEN_QUAL
    assert get_parsed_args(f"file1.mmcif -t pdb -w {const.CONVERTER_OB} --coord-gen Gen3D"
                           ).d_converter_args[COORD_GEN_QUAL_KEY] == DEFAULT_COORD_GEN_QUAL


def test_fail_log_file_not_set():
    """Test that trying to get the log file raises an exception due to the test file not existing"""
    args = get_parsed_args(f"file1.mmcif -t pdb -w {const.CONVERTER_OB}")
    with pytest.raises(FileConverterInputException) as e:
        _ = args.log_file
    assert _compressed_match(f"Input file '{os.getcwd()}/file1.mmcif' cannot be found", e.value)


def test_log_file_list_mode():
    """Test that the log file uses the expected default value in list mode"""
    args = get_parsed_args("--list")
    assert args.log_file == const.DEFAULT_LISTING_LOG_FILE


def _check_no_errors(captured):
    """Check that no errors were produced in output"""
    assert not captured.err
    assert "Traceback" not in captured.out


@pytest.mark.parametrize("auto_str", ["", "-w auto", "-w Auto", "--with AUTO"])
def test_auto_converter(auto_str):
    """Ensure that a converter can be properly determined automatically
    """

    # Test that Open Babel is chosen when expected
    args = get_parsed_args(f"file1.pdb -f pdb-0 -t inchi {auto_str}")
    assert args.name == regularize_name(const.CONVERTER_OB)

    # Test that c2x is chosen when expected
    args = get_parsed_args(f"file1.pdb -f pdb-0 -t xyz-0 {auto_str}")
    assert args.name == regularize_name(const.CONVERTER_C2X)


def test_auto_ambiguous_from_format():
    """Ensure that the proper error is raised if the input format is ambiguous when using 'auto' converter
    """
    with pytest.raises(FileConverterInputException) as e:
        get_parsed_args("file1 -f pdb -t xyz-0 -w auto")
    assert _compressed_match("the input format determined from the extension of the input file or specified "
                             "with `-f/--from` must unambiguously", e.value)


def test_auto_ambiguous_ext():
    """Ensure that the proper error is raised if the input format is ambiguous when using 'auto' converter
    """
    with pytest.raises(FileConverterInputException) as e:
        get_parsed_args("file1.pdb -t xyz-0 -w auto")
    assert _compressed_match("the input format determined from the extension of the input file or specified "
                             "with `-f/--from` must unambiguously", e.value)


def test_auto_multi_ambiguous_ext():
    """Ensure that the proper error is raised if the input format is ambiguous for one or more in a list of input files
    when using 'auto' converter
    """
    with pytest.raises(FileConverterInputException) as e:
        get_parsed_args("file1.pdb file2.cif -t xyz-0 -w auto")
    assert _compressed_match("input format must be uniquely identifiable for all input files.", e.value)


def test_auto_invalid_to_format():
    """Ensure that the proper error is raised if the output format is invalid when using 'auto' converter
    """
    with pytest.raises(FileConverterInputException) as e:
        get_parsed_args("file1 -f pdb-0 -t invalid_format -w auto")
    assert _compressed_match("is not recognised as a valid output format. To see supported formats", e.value)


def test_auto_ambiguous_to_format():
    """Ensure that the proper error is raised if the to format is ambiguous when using 'auto' converter
    """
    with pytest.raises(FileConverterInputException) as e:
        get_parsed_args("file1 -f pdb-0 -t xyz -w auto")
    assert _compressed_match("is ambiguous and can correspond to multiple possible output formats", e.value)


def test_auto_no_common_converter():
    """Ensure that the proper error is raised if no one converter can perform all conversions when using 'auto'
    converter.
    """
    with pytest.raises(FileConverterInputException) as e:
        get_parsed_args("file1.abi file2.inchi -t pdb-0 -w auto")
    assert _compressed_match("No converter is available which can perform a conversion of all input files", e.value)


def test_list_converters(capsys):
    """Test the option to list available converters
    """
    run_with_arg_string("--list")
    captured = capsys.readouterr()
    assert "Available converters" in captured.out
    for converter_rname in L_REGISTERED_CONVERTERS:
        converter_name = get_registered_converter_class(converter_rname).meta.name
        assert converter_name in captured.out, converter_name

    _check_no_errors(captured)


@pytest.mark.parametrize("converter_name", L_REGISTERED_CONVERTERS)
def test_detail_converter(capsys, converter_name):
    """Test the option to provide detail on a converter
    """

    converter_info = get_converter_info(converter_name)

    run_with_arg_string(f"--list {converter_name}")
    captured = capsys.readouterr()

    assert _compressed_match(converter_info.pretty_name, captured.out)

    if not converter_info.description:
        assert _compressed_match("available for this converter", captured.out)
    else:
        assert _compressed_match(converter_info.description, captured.out)

    # Check for URL
    assert converter_info.url in captured.out

    # Check for list of allowed input/output formats
    assert "    INPUT    OUTPUT    DESCRIPTION" in captured.out

    l_allowed_in_formats, l_allowed_out_formats = get_possible_formats(converter_name)
    for in_format in l_allowed_in_formats:
        output_allowed = "yes" if in_format in l_allowed_out_formats else "no"
        assert _compressed_match(f"{in_format.disambiguated_name}yes{output_allowed}{in_format.description}",
                                 captured.out)
    for out_format in l_allowed_out_formats:
        input_allowed = "yes" if out_format in l_allowed_in_formats else "no"
        assert _compressed_match(f"{out_format.disambiguated_name}{input_allowed}yes{out_format.description}",
                                 captured.out)

    _check_no_errors(captured)


def test_detail_converter_bad_name(capsys):
    """Test we do get a simple error for a bad converter name
    """
    with pytest.raises(SystemExit):
        run_with_arg_string("--list bad_converter")
    captured = capsys.readouterr()
    assert "not recognized" in captured.err
    assert "Traceback" not in captured.out
    assert "Traceback" not in captured.err


def test_detail_converter_with(capsys):
    """Test that we can also provide the converter name with -w/--with
    """
    run_with_arg_string(f"-l -w {const.CONVERTER_C2X}")
    captured = capsys.readouterr()
    _check_no_errors(captured)
    assert const.CONVERTER_C2X in captured.out
    assert const.CONVERTER_OB not in captured.out


def test_get_conversions(capsys):
    """Test the option to get information on converters which can perform a desired conversion
    """
    in_format = "xyz-1"
    out_format = "inchi"
    l_conversions = get_possible_conversions(in_format, out_format)

    run_with_arg_string(f"-l -f {in_format} -t {out_format}")
    captured = capsys.readouterr()

    _check_no_errors(captured)

    assert bool(l_conversions) == _compressed_match("The following registered converters can convert from "
                                                    f"{in_format} to {out_format}:", captured.out)

    for converter_info, _, _ in l_conversions:
        if converter_info.name in L_REGISTERED_CONVERTERS:
            assert _compressed_match(converter_info.pretty_name, captured.out)
    for name in L_REGISTERED_CONVERTERS:
        converter_info = get_converter_info(name)
        if converter_info not in [x[0] for x in l_conversions]:
            assert not _compressed_match(converter_info.pretty_name, captured.out)


def test_list_chain(capsys):
    """Test the ability to get a pathway for a chained conversion
    """
    in_format = get_format_info(FORMAT_MOLDY)
    out_format = get_format_info(FORMAT_INCHI)
    pathway = get_conversion_pathway(in_format, out_format)
    assert len(pathway) > 1

    run_with_arg_string(f"-l -f {in_format.id} -t {out_format.id}")
    captured = capsys.readouterr()

    _check_no_errors(captured)

    assert _compressed_match(f"No direct conversions are possible from {in_format.format_word()} to "
                             f"{out_format.format_word()}", captured.out)

    assert _compressed_match(f"A chained conversion is possible from {in_format.format_word()} to "
                             f"{out_format.format_word()} using registered converters:", captured.out)

    for i, step in enumerate(pathway):
        assert _compressed_match(f"{i+1}) Convert from {step[1].format_word()} to {step[2].format_word()} with "
                                 f"{step[0].format_word()}", captured.out)


def test_list_chain_impossible(capsys):
    """Test that we get the expected output when a chained conversion is not possible
    """

    in_format = "cif"
    out_format = "abinit"

    run_with_arg_string(f"-l -f {in_format} -t {out_format}")
    captured = capsys.readouterr()

    assert _compressed_match(f"No chained conversions are possible from {in_format} to {out_format}.", captured.out)

    # Check that igraph's warning is suppressed
    assert not _compressed_match("Couldn't reach some vertices", captured.out)

    _check_no_errors(captured)


def test_conversion_info_open_babel(capsys):
    """Test that we get the expected information on the 'Open Babel' converter
    """

    converter_name = const.CONVERTER_OB

    in_format = "xyz-1"
    out_format = "inchi"
    qual = get_conversion_quality(converter_name, in_format, out_format)

    # Test a basic listing of arguments, checking with the converter name in lowercase to be sure that works
    run_with_arg_string(f"-l {converter_name.lower()} -f {in_format} -t {out_format}")
    captured = capsys.readouterr()

    _check_no_errors(captured)

    # Check that conversion quality details are in the output as expected
    assert _compressed_match(f"Conversion from {in_format} to {out_format} with {converter_name} is "
                             f"possible with {qual.qual_str} conversion quality", captured.out)
    assert _compressed_match("WARNING: Potential data loss or extrapolation issues with this conversion:",
                             captured.out)
    assert _compressed_match(const.QUAL_NOTE_OUT_MISSING.format(const.QUAL_2D_LABEL), captured.out)
    assert _compressed_match(const.QUAL_NOTE_OUT_MISSING.format(const.QUAL_3D_LABEL), captured.out)
    assert _compressed_match(const.QUAL_NOTE_IN_MISSING.format(const.QUAL_CONN_LABEL), captured.out)

    l_in_flags, l_in_options = get_in_format_args(converter_name, in_format)
    l_out_flags, l_out_options = get_out_format_args(converter_name, out_format)

    # Check headings for input/output flags/options are present if and only if some of those flags/options exist
    assert bool(l_in_flags) == _compressed_match(f"Allowed input flags for format '{in_format}'", captured.out)
    assert bool(l_out_flags) == _compressed_match(f"Allowed output flags for format '{out_format}'", captured.out)
    assert bool(l_in_options) == _compressed_match(f"Allowed input options for format '{in_format}'", captured.out)
    assert bool(l_out_options) == _compressed_match(f"Allowed output options for format '{out_format}'", captured.out)

    # Check that info for each flag and option is printed as expected
    for flag_info in l_in_flags + l_out_flags:
        info = flag_info.info if flag_info.info and flag_info.info != "N/A" else ""
        assert _compressed_match(f"{flag_info.name}{flag_info.description}{info}", captured.out)
    for option_info in l_in_options + l_out_options:
        info = option_info.info if option_info.info and option_info.info != "N/A" else ""
        assert _compressed_match(f"{option_info.name}<{option_info.brief}>{option_info.description}{info}",
                                 captured.out)


@pytest.mark.parametrize("converter_name", [const.CONVERTER_C2X, const.CONVERTER_ATO])
def test_conversion_info_others(capsys, converter_name):
    """Test that we get the expected information on other converters
    """

    in_format = "pdb-0"
    out_format = "cif"
    qual = get_conversion_quality(converter_name, in_format, out_format)

    run_with_arg_string(f"-l {converter_name} -f {in_format} -t {out_format}")

    captured = capsys.readouterr()
    _check_no_errors(captured)

    # Check that conversion quality details are in the output as expected
    assert _compressed_match(f"Conversion from {in_format} to {out_format} with {converter_name} is "
                             f"possible with {qual.qual_str} conversion quality", captured.out)
    assert _compressed_match("WARNING: Potential data loss or extrapolation issues with this conversion:",
                             captured.out)
    assert _compressed_match(const.QUAL_NOTE_OUT_MISSING.format(const.QUAL_CONN_LABEL), captured.out)


def test_format_info(capsys):
    """Test that we can successfully get information on a file format"""

    # Try to get info on an unambiguous format

    in_format = "inchi"
    in_format_info = get_format_info(in_format)
    run_with_arg_string(f"-l -f {in_format}")

    captured = capsys.readouterr()

    _check_no_errors(captured)

    # Check for basic format information
    assert _compressed_match(f"{in_format_info.disambiguated_name} (ID {in_format_info.id}): " +
                             in_format_info.description, captured.out)

    # Check for property information
    for attr, label in D_FORMAT_PROPERTY_ATTRS.items():
        support_status = getattr(in_format_info, attr)
        if support_status:
            assert _compressed_match(label + " supported", captured.out)
        elif support_status is False:
            assert _compressed_match(label + " not supported", captured.out)
        else:
            assert _compressed_match(label + " unknown whether or not to be supported", captured.out)


def test_format_info_ambiguous(capsys):
    """Test that we get expected information for an ambiguous format"""

    out_format = "pdb"
    l_out_format_info = get_format_info(out_format, which="all")
    run_with_arg_string(f"-l -t {out_format}")

    captured = capsys.readouterr()

    _check_no_errors(captured)

    assert _compressed_match(f"WARNING: Format '{out_format}' is ambiguous", captured.out)

    for out_format_info in l_out_format_info:
        assert _compressed_match(out_format_info.format_oneline(), captured.out)


def test_format_info_in_unrecognised(capsys):
    """Test we get expected errors for unrecognised input format"""

    in_format = 99999
    with pytest.raises(SystemExit):
        run_with_arg_string(f"-l -f {in_format}")

    assert _compressed_match(f"ERROR: Format '{in_format}' not recognised", capsys.readouterr().err)


def test_format_info_out_unrecognised(capsys):
    """Test we get expected errors for unrecognised output format"""

    out_format = "not_a_format"

    with pytest.raises(SystemExit):
        run_with_arg_string(f"-l -t {out_format}")

    assert _compressed_match(f"ERROR: Format '{out_format}' not recognised", capsys.readouterr().err)
