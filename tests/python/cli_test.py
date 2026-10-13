"""@file psdi-data-conversion/tests/cli_test.py

Created 2025-01-15 by Bryan Gillis.

Tests of the command-line interface
"""

import logging
import os
import shlex
import sys
from itertools import product
from tempfile import TemporaryDirectory
from unittest.mock import patch

import pytest

from psdi_data_conversion import constants as const
from psdi_data_conversion.converter import D_CONVERTER_ARGS, L_REGISTERED_CONVERTERS, get_registered_converter_class
from psdi_data_conversion.converters.openbabel.converter import (COORD_GEN_KEY, COORD_GEN_QUAL_KEY, DEFAULT_COORD_GEN,
                                                                 DEFAULT_COORD_GEN_QUAL)
from psdi_data_conversion.database import (D_FORMAT_PROPERTY_ATTRS, MSG_CONVERSION_ONELINE, Conversion, format_weight,
                                           get_conversion_pathway, get_conversion_quality, get_converter_info,
                                           get_format_info, get_in_format_args, get_out_format_args,
                                           get_possible_conversion_pathways, get_possible_conversions,
                                           get_possible_formats)
from psdi_data_conversion.main import FileConverterInputException, parse_args
from psdi_data_conversion.testing.constants import FORMAT_INCHI, FORMAT_MOL, FORMAT_MOLDY, FORMAT_PDB_0
from psdi_data_conversion.testing.conversion_test_specs import l_cli_test_specs
from psdi_data_conversion.testing.utils import compressed_match, run_test_conversion_with_cli, run_with_arg_string
from psdi_data_conversion.utils import regularize_name


def test_unique_args(subtests):
    """Check that all converter-specific arguments have unique names
    """
    s_arg_names = set()
    for name in L_REGISTERED_CONVERTERS:
        for arg_name, _, _ in D_CONVERTER_ARGS[name]:
            with subtests.test("arg is unique", conv=name, arg=arg_name):
                assert arg_name not in s_arg_names, ("Name clash between converters, with multiple using the argument "
                                                     f"'{arg_name}'")
            s_arg_names.add(arg_name)


def _get_parsed_args(s):
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


@pytest.mark.parametrize("test_spec", l_cli_test_specs,
                         ids=lambda x: x.name)
def test_conversions(test_spec, subtests):
    """Run all conversion tests in the defined list of test specifications
    """
    run_test_conversion_with_cli(test_spec, subtests)


def test_general_arg_parsing(subtests):
    """Test that a standard argument string is parsed properly
    """
    cwd = os.getcwd()
    args = _get_parsed_args(f"file1 file2 -f mmcif -i {cwd} -t pdb -o {cwd}/.. -w '{const.CONVERTER_C2X}' " +
                            r"--delete-input --from-flags '\-ab \-c \--example' --to-flags '\-d' " +
                            r"--from-options '-x xval --xopt xoptval' --to-options '-y yval --yopt yoptval' "
                            "--strict --nc --coord-gen Gen3D best -q --log-file text.log")

    with subtests.test("input files are correct"):
        assert args.l_args[0] == "file1"
        assert args.l_args[1] == "file2"

    with subtests.test("input dir is correct"):
        assert args.input_dir == cwd

    with subtests.test("output format is correct"):
        assert args.to_format == "pdb"

    with subtests.test("output dir is correct"):
        assert args.output_dir == f"{cwd}/.."

    with subtests.test("converter is correct"):
        assert args.converter == const.CONVERTER_C2X

    with subtests.test("chain mode is correct"):
        assert args.chain is False

    with subtests.test("chain path is correct"):
        assert args.path is None

    with subtests.test("`no_check` is correct"):
        assert args.no_check is True

    with subtests.test("strict mode is correct"):
        assert args.strict is True

    with subtests.test("`delete_input` is correct"):
        assert args.delete_input is True

    with subtests.test("`from_flags` are correct"):
        assert args.from_flags == "-ab -c --example"

    with subtests.test("`to_flags` are correct"):
        assert args.to_flags == "-d"

    with subtests.test("`from_options` are correct"):
        assert args.from_options == "-x xval --xopt xoptval"

    with subtests.test("`to_options` are correct"):
        assert args.to_options == "-y yval --yopt yoptval"

    with subtests.test("quiet mode is correct"):
        assert args.quiet is True

    with subtests.test("log file is correct"):
        assert args.log_file == "text.log"

    with subtests.test("log mode is correct"):
        assert args.log_mode == const.LOG_NONE


def test_open_babel_args(subtests):
    """Test that Open-Babel-specific arguments are parsed correctly
    """
    args = _get_parsed_args(f"file1.mmcif -t pdb -w '{const.CONVERTER_OB}' --coord-gen Gen3D best")

    with subtests.test("coord. gen mode is correct"):
        assert args.d_converter_args[COORD_GEN_KEY] == "Gen3D"

    with subtests.test("coord. gen quality is correct"):
        assert args.d_converter_args[COORD_GEN_QUAL_KEY] == "best"


def test_fail_no_args(subtests):
    """Test that the parsing fails if no arguments are provided"""
    with subtests.test("get input exception"), pytest.raises(FileConverterInputException) as e:
        _get_parsed_args("")
    with subtests.test("correct error message"):
        assert compressed_match("One or more names of files to convert must be provided", e.value)


def test_fail_no_to_format(subtests):
    """Test that the parsing fails if output format isn't specified"""
    with subtests.test("get input exception"), pytest.raises(FileConverterInputException) as e:
        _get_parsed_args("file1.mmcif")
    with subtests.test("correct error message"):
        assert compressed_match("Output format must be provided", e.value)


def test_fail_no_input_dir(subtests):
    """Test that the parsing fails if the input directory doesn't exist"""
    with subtests.test("get input exception"), pytest.raises(FileConverterInputException) as e:
        _get_parsed_args("file1.mmcif -i /no/where -t pdb")
    with subtests.test("correct error message"):
        assert compressed_match("The provided input directory '/no/where' does not exist as a directory", e.value)


def test_fail_invalid_converter(subtests):
    """Test that the parsing fails if the converter isn't recognized"""
    with subtests.test("get input exception"), pytest.raises(FileConverterInputException) as e:
        _get_parsed_args("file1.mmcif -t pdb -w FakeConverter")
    with subtests.test("correct error message"):
        assert compressed_match("Converter 'fakeconverter' not recognised", e.value)


def test_fail_bad_coord_gen_type(subtests):
    """Test that the parsing fails with bad --coord-gen type"""
    with subtests.test("get input exception"), pytest.raises(FileConverterInputException) as e:
        _get_parsed_args(f"file1.mmcif -t pdb -w '{const.CONVERTER_OB}' --coord-gen Gen1D")
    with subtests.test("correct error message"):
        assert compressed_match("Coordinate generation type 'Gen1D' not recognised.", e.value)


def test_fail_bad_coord_gen_quality(subtests):
    """Test that the parsing fails with bad --coord-gen quality"""
    with subtests.test("get input exception"), pytest.raises(FileConverterInputException) as e:
        _get_parsed_args(f"file1.mmcif -t pdb -w '{const.CONVERTER_OB}' --coord-gen Gen3D worst")
    with subtests.test("correct error message"):
        assert compressed_match("Coordinate generation quality 'worst' not recognised.", e.value)


def test_fail_bad_coord_gen_len(subtests):
    """Test that the parsing fails too many args to --coord-gen"""
    with subtests.test("get input exception"), pytest.raises(FileConverterInputException) as e:
        _get_parsed_args(f"file1.mmcif -t pdb -w '{const.CONVERTER_OB}' --coord-gen Gen3D best quality")
    with subtests.test("correct error message"):
        assert compressed_match("At most two arguments may be provided to `--coord-gen`", e.value)


def test_fail_bad_logging_mode(subtests):
    """Test that the parsing fails if it doesn't recognise the logging mode"""
    with subtests.test("get input exception"), pytest.raises(FileConverterInputException) as e:
        _get_parsed_args(f"file1.mmcif -t pdb -w '{const.CONVERTER_OB}' --log-mode max")
    with subtests.test("correct error message"):
        assert compressed_match("Unrecognised logging mode: 'max'", e.value)


def test_list_args(subtests):
    """Test that the parsing works if we just ask for a list, and set log mode to stdout"""
    args = _get_parsed_args("--list")
    assert args.list
    assert args.log_mode == const.LOG_STDOUT


def test_list_converter(subtests):
    """Test that the parsing works if we ask for info on a specific converter"""
    args = _get_parsed_args("-l Open Babel")
    assert args.converter == regularize_name("Open Babel")
    args = _get_parsed_args("--list 'Open Babel'")
    assert args.converter == regularize_name("Open Babel")
    args = _get_parsed_args("-l Atomsk")
    assert args.converter == regularize_name("Atomsk")


def test_converter_input(subtests):
    """Test that the converter specified with -w/--with is properly parsed
    """
    args = _get_parsed_args(f"file1.mmcif -t pdb -w {const.CONVERTER_OB}")
    assert args.converter == regularize_name(const.CONVERTER_OB)
    args = _get_parsed_args(f"file1.mmcif -t pdb -w '{const.CONVERTER_OB}'")
    assert args.converter == regularize_name(const.CONVERTER_OB)


def test_default_input_dir():
    """Test that the input directory is set to the current directory if not specified"""
    args = _get_parsed_args(f"file1.mmcif -t pdb -w {const.CONVERTER_OB}")
    assert args.input_dir == os.getcwd()


def test_default_output_dir():
    """Test that the output dir defaults to match input dir"""
    cwd = os.getcwd()
    args = _get_parsed_args(f"file1.mmcif -i {cwd}/.. -t pdb -w {const.CONVERTER_OB}")
    assert args.output_dir == f"{cwd}/.."


def test_default_coord_gen(subtests):
    """Test that we get the default coordinate generation options if they aren't explicitly specified"""
    args = _get_parsed_args(f"file1.mmcif -t pdb -w {const.CONVERTER_OB}")
    assert args.d_converter_args[COORD_GEN_KEY] == DEFAULT_COORD_GEN
    assert args.d_converter_args[COORD_GEN_QUAL_KEY] == DEFAULT_COORD_GEN_QUAL
    assert _get_parsed_args(f"file1.mmcif -t pdb -w {const.CONVERTER_OB} --coord-gen Gen3D"
                            ).d_converter_args[COORD_GEN_QUAL_KEY] == DEFAULT_COORD_GEN_QUAL


def test_fail_log_file_not_set(subtests):
    """Test that trying to get the log file raises an exception due to the test file not existing"""
    args = _get_parsed_args(f"file1.mmcif -t pdb -w {const.CONVERTER_OB}")
    with subtests.test("get input exception"), pytest.raises(FileConverterInputException) as e:
        _ = args.log_file
    with subtests.test("correct error message"):
        assert compressed_match(f"Input file '{os.getcwd()}/file1.mmcif' cannot be found", e.value)


def test_log_file_list_mode():
    """Test that the log file uses the expected default value in list mode"""
    args = _get_parsed_args("--list")
    assert args.log_file == const.DEFAULT_LISTING_LOG_FILE


def _check_no_errors(captured, subtests):
    """Check that no errors were produced in output"""
    with subtests.test("There's no output to stderr"):
        assert not captured.err
    with subtests.test("There aren't any tracebacks in stdout"):
        assert "Traceback" not in captured.out


@pytest.mark.parametrize("auto_str", ["", "-w auto", "-w Auto", "--with AUTO"])
def test_auto_converter(auto_str, subtests):
    """Ensure that a converter can be properly determined automatically
    """

    # Test that Open Babel is chosen when expected
    args = _get_parsed_args(f"file1.pdb -f pdb-0 -t inchi {auto_str}")
    assert args.converter == regularize_name(const.CONVERTER_OB)

    # Test that c2x is chosen when expected
    args = _get_parsed_args(f"file1.pdb -f pdb-0 -t xyz-0 {auto_str}")
    assert args.converter == regularize_name(const.CONVERTER_C2X)


@pytest.mark.parametrize("converter_name", [const.CONVERTER_AUTO, const.CONVERTER_AUTOCHAIN])
def test_auto_ambiguous_from_format(converter_name, subtests):
    """Ensure that the proper error is raised if the input format is ambiguous when using 'auto(chain)' converter
    """
    with subtests.test("get input exception"), pytest.raises(FileConverterInputException) as e:
        _get_parsed_args(f"file1 -f pdb -t xyz-0 -w {converter_name}")
    with subtests.test("correct error message"):
        assert compressed_match("the input format determined from the extension of the input file or specified "
                                "with `-f/--from` must unambiguously", e.value)


@pytest.mark.parametrize("converter_name", [const.CONVERTER_AUTO, const.CONVERTER_AUTOCHAIN])
def test_auto_ambiguous_ext(converter_name, subtests):
    """Ensure that the proper error is raised if the input format is ambiguous when using 'auto(chain)' converter
    """
    with subtests.test("get input exception"), pytest.raises(FileConverterInputException) as e:
        _get_parsed_args(f"file1.pdb -t xyz-0 -w {converter_name}")
    with subtests.test("correct error message"):
        assert compressed_match("the input format determined from the extension of the input file or specified "
                                "with `-f/--from` must unambiguously", e.value)


@pytest.mark.parametrize("converter_name", [const.CONVERTER_AUTO, const.CONVERTER_AUTOCHAIN])
def test_auto_archive_no_from(converter_name, subtests):
    """Ensure that the proper error is raised if the input format is ambiguous when using 'auto(chain)' converter
    """
    with subtests.test("get input exception"), pytest.raises(FileConverterInputException) as e:
        _get_parsed_args(f"file1.tar.gz -t xyz-0 -w {converter_name}")
    with subtests.test("correct error message"):
        assert compressed_match("Cannot determine input format for file", e.value)


def test_auto_multi_ambiguous_ext(subtests):
    """Ensure that the proper error is raised if the input format is ambiguous for one or more in a list of input files
    when using 'auto' converter
    """
    with subtests.test("get input exception"), pytest.raises(FileConverterInputException) as e:
        _get_parsed_args(f"file1.pdb file2.cif -t xyz-0 -w {const.CONVERTER_AUTO}")
    with subtests.test("correct error message"):
        assert compressed_match("input format must be uniquely identifiable for all input files.", e.value)


@pytest.mark.parametrize("converter_name", [const.CONVERTER_AUTO, const.CONVERTER_AUTOCHAIN])
def test_auto_invalid_to_format(converter_name, subtests):
    """Ensure that the proper error is raised if the output format is invalid when using 'auto(chain)' converter
    """
    with subtests.test("get input exception"), pytest.raises(FileConverterInputException) as e:
        _get_parsed_args(f"file1 -f pdb-0 -t invalid_format -w {converter_name}")
    with subtests.test("correct error message"):
        assert compressed_match("is not recognised as a valid output format. To see supported formats", e.value)


@pytest.mark.parametrize("converter_name", [const.CONVERTER_AUTO, const.CONVERTER_AUTOCHAIN])
def test_auto_ambiguous_to_format(converter_name, subtests):
    """Ensure that the proper error is raised if the to format is ambiguous when using 'auto(chain)' converter
    """
    with subtests.test("get input exception"), pytest.raises(FileConverterInputException) as e:
        _get_parsed_args(f"file1 -f pdb-0 -t xyz -w {converter_name}")
    with subtests.test("correct error message"):
        assert compressed_match("is ambiguous and can correspond to multiple possible output formats", e.value)


def test_auto_no_common_converter(subtests):
    """Ensure that the proper error is raised if no one converter can perform all conversions when using 'auto'
    converter.
    """
    with subtests.test("get input exception"), pytest.raises(FileConverterInputException) as e:
        _get_parsed_args(f"file1.abi file2.inchi -t pdb-0 -w {const.CONVERTER_AUTO}")
    with subtests.test("correct error message"):
        assert compressed_match("No converter is available which can perform a conversion of all input files", e.value)


@pytest.mark.parametrize("chain_str", const.L_CONVERTER_AUTOCHAIN)
def test_autochain_args(chain_str, subtests):
    """Test that the parsing correctly detects when an automatic chain is requested"""
    args = _get_parsed_args(f"file1.mmcif -t pdb-0 -w '{chain_str}'")
    assert args.converter == const.CONVERTER_AUTOCHAIN
    assert args.chain is True


l_arg_strs = [f"--{x[0]}-{x[1]} foo" for x in product(("to", "from"), ("flags", "options"))]
l_converter_and_arg = list(product(("auto", "autochain"), l_arg_strs))


@pytest.mark.parametrize("converter, arg_str", l_converter_and_arg)
def test_auto_no_converter_args(converter, arg_str, subtests):
    """Test that an error is raised if any converter-specific arguments are provided when an automatic converter or
    chain is requested"""
    with subtests.test("get input exception"), pytest.raises(FileConverterInputException) as e:
        _get_parsed_args(f"file1 -f pdb-0 -t xyz-0 -w {converter} {arg_str}")
    with subtests.test("correct error message"):
        assert compressed_match("Converter-specific arguments cannot be provided when the converter or chain is "
                                "automatically-determined.", e.value)


l_path_args = [(FORMAT_MOLDY, f"openbabel {FORMAT_PDB_0} atomsk", FORMAT_INCHI),
               (FORMAT_MOLDY, f"openbabel {FORMAT_PDB_0} atomsk {FORMAT_INCHI}", None),
               (None, f"{FORMAT_MOLDY} openbabel {FORMAT_PDB_0} atomsk", FORMAT_INCHI),
               (None, f"{FORMAT_MOLDY} openbabel {FORMAT_PDB_0} atomsk {FORMAT_INCHI}", None)]


@pytest.mark.parametrize("from_format, path_str, to_format", l_path_args)
def test_path_args(from_format: int | None, path_str: str, to_format: int | None):
    """Test that the parsing correctly interprets a path when one is requested"""
    arg_str = f"file1 --path {path_str}"
    if from_format:
        arg_str += f" -f {from_format}"
    if to_format:
        arg_str += f" -t {to_format}"
    args = _get_parsed_args(arg_str)
    assert args.path


def test_path_consistent_from_format():
    """Test that no error is raised if the `from_format` matches the first in the `path`"""
    args = _get_parsed_args(f"file1 -f {FORMAT_MOLDY} --path {FORMAT_MOLDY} openbabel {FORMAT_PDB_0} atomsk "
                            f"{FORMAT_INCHI}")
    assert args.path


def test_path_inconsistent_from_format(subtests):
    """Test that an error is raised if the `from_format` doesn't match the first in the `path`"""
    with subtests.test("get input exception"), pytest.raises(FileConverterInputException) as e:
        _get_parsed_args(f"file1 -f {FORMAT_PDB_0} --path {FORMAT_MOLDY} openbabel {FORMAT_PDB_0} atomsk "
                         f"{FORMAT_INCHI}")
    with subtests.test("correct error message"):
        assert compressed_match(f"The format '{FORMAT_PDB_0}' provided to `-f/--from` does not match '{FORMAT_MOLDY}"
                                f"', the first format provided to `--path`", e.value)


def test_path_consistent_to_format():
    """Test that no error is raised if the `to_format` matches the last in the `path`"""
    args = _get_parsed_args(f"file1 -t {FORMAT_INCHI} --path {FORMAT_MOLDY} openbabel {FORMAT_PDB_0} atomsk "
                            f"{FORMAT_INCHI}")
    assert args.path


def test_path_inconsistent_to_format(subtests):
    """Test that an error is raised if the `to_format` doesn't match the last in the `path`"""
    with subtests.test("get input exception"), pytest.raises(FileConverterInputException) as e:
        _get_parsed_args(f"file1 -t {FORMAT_PDB_0} --path {FORMAT_MOLDY} openbabel {FORMAT_PDB_0} atomsk "
                         f"{FORMAT_INCHI}")
    with subtests.test("correct error message"):
        assert compressed_match(f"The format '{FORMAT_PDB_0}' provided to `-t/--to` does not match '{FORMAT_INCHI}"
                                f"', the last format provided to `--path`", e.value)


def test_path_ambiguous_format(subtests):
    """Test that an error is raised if a format in `path` is ambiguous"""
    with subtests.test("get input exception"), pytest.raises(FileConverterInputException) as e:
        _get_parsed_args(f"file1 --path {FORMAT_MOLDY} openbabel pdb atomsk {FORMAT_INCHI}")
    with subtests.test("correct error message"):
        assert compressed_match("'pdb' is ambiguous and can correspond to multiple possible formats.", e.value)


@pytest.mark.parametrize("bad_format, bad_converter", [(True, False),
                                                       (False, True),
                                                       (True, True)])
def test_path_bad_format_converter(bad_format: bool, bad_converter: bool, subtests):
    """Test that an error is raised if a format and/or converter in `path` is invalid"""
    file_format = "INVALID" if bad_format else FORMAT_MOLDY
    converter = "INVALID" if bad_converter else "openbabel"
    with subtests.test("get input exception"), pytest.raises(FileConverterInputException) as e:
        _get_parsed_args(f"file1 --path {FORMAT_MOLDY} {converter} {file_format} atomsk {FORMAT_INCHI}")

    if bad_format:
        assert compressed_match("'INVALID' is not recognised as a valid format", e.value)

    if bad_converter:
        assert compressed_match("'INVALID' is not recognised as a valid converter", e.value)


def test_path_unsupported_converter(subtests):
    """Test that an error is raised if a converter in `path` is unsupported"""
    with subtests.test("get input exception"), pytest.raises(FileConverterInputException) as e:
        _get_parsed_args(f"file1 --path {FORMAT_MOLDY} molconverter {FORMAT_PDB_0} atomsk {FORMAT_INCHI}")
    with subtests.test("correct error message"):
        assert compressed_match("Molconverter is not supported", e.value)


def test_path_bad_length(subtests):
    """Test that an error is raised if the `path` has a bad length, e.g. due to spaces being used in a converter name"""
    with subtests.test("get input exception"), pytest.raises(FileConverterInputException) as e:
        _get_parsed_args(f"file1 --path {FORMAT_MOLDY} Open Babel {FORMAT_PDB_0} atomsk {FORMAT_INCHI}")
    with subtests.test("correct error message"):
        assert compressed_match("invalid due to an incorrect number of elements", e.value)


def test_list_converters(capsys, subtests):
    """Test the option to list available converters
    """
    run_with_arg_string("--list")
    captured = capsys.readouterr()
    assert "Available converters" in captured.out
    for converter_rname in L_REGISTERED_CONVERTERS:
        converter_name = get_registered_converter_class(converter_rname).meta.name
        assert converter_name in captured.out, converter_name

    _check_no_errors(captured, subtests)


@pytest.mark.parametrize("converter_name", L_REGISTERED_CONVERTERS)
def test_detail_converter(capsys, converter_name, subtests):
    """Test the option to provide detail on a converter
    """

    converter_info = get_converter_info(converter_name)

    run_with_arg_string(f"--list {converter_name}")
    captured = capsys.readouterr()

    with subtests.test("converter name displayed"):
        assert compressed_match(converter_info.pretty_name, captured.out)

    if not converter_info.description:
        assert compressed_match("available for this converter", captured.out)
    else:
        assert compressed_match(converter_info.description, captured.out)

    # Check for URL
    assert converter_info.url in captured.out

    # Check for list of allowed input/output formats
    assert "    INPUT    OUTPUT    DESCRIPTION" in captured.out

    l_allowed_in_formats, l_allowed_out_formats = get_possible_formats(converter_name)
    for in_format in l_allowed_in_formats:
        output_allowed = "yes" if in_format in l_allowed_out_formats else "no"
        assert compressed_match(f"{in_format.disambiguated_name}yes{output_allowed}{in_format.description}",
                                captured.out)
    for out_format in l_allowed_out_formats:
        input_allowed = "yes" if out_format in l_allowed_in_formats else "no"
        assert compressed_match(f"{out_format.disambiguated_name}{input_allowed}yes{out_format.description}",
                                captured.out)

    _check_no_errors(captured, subtests)


def test_detail_converter_bad_name(capsys, subtests):
    """Test we do get a simple error for a bad converter name
    """
    with pytest.raises(SystemExit):
        run_with_arg_string("--list bad_converter")
    captured = capsys.readouterr()
    assert "not recognised" in captured.err
    assert "Traceback" not in captured.out
    assert "Traceback" not in captured.err


def test_detail_converter_with(capsys, subtests):
    """Test that we can also provide the converter name with -w/--with
    """
    run_with_arg_string(f"-l -w {const.CONVERTER_C2X}")
    captured = capsys.readouterr()
    _check_no_errors(captured, subtests)
    assert const.CONVERTER_C2X in captured.out
    assert const.CONVERTER_OB not in captured.out


def test_get_conversions(capsys, subtests):
    """Test the option to get information on converters which can perform a desired conversion
    """
    in_format = "xyz-1"
    out_format = "inchi"
    l_conversions = get_possible_conversions(in_format, out_format)

    run_with_arg_string(f"-l -f {in_format} -t {out_format}")
    captured = capsys.readouterr()

    _check_no_errors(captured, subtests)

    assert bool(l_conversions) == compressed_match("The following registered converters can convert from "
                                                   f"{in_format} to {out_format}:", captured.out)

    for converter_info, _, _ in l_conversions:
        if converter_info.name in L_REGISTERED_CONVERTERS:
            assert compressed_match(converter_info.pretty_name, captured.out)
    for name in L_REGISTERED_CONVERTERS:
        converter_info = get_converter_info(name)
        if converter_info not in [x[0] for x in l_conversions]:
            assert not compressed_match(converter_info.pretty_name, captured.out)


def _check_step_details_present(step: Conversion, captured, subtests, **kwargs):
    with subtests.test("Display each step of the chain properly", **kwargs):
        assert compressed_match(MSG_CONVERSION_ONELINE.format(step.in_format.format_word(),
                                                              step.out_format.format_word(),
                                                              step.converter.format_word(),
                                                              format_weight(step.get_weight())),
                                captured.out)


def test_list_chain_if_not_direct(capsys, subtests):
    """Test the ability to get a pathway for a chained conversion if no direct conversion is possible
    """
    in_format = get_format_info(FORMAT_MOLDY)
    out_format = get_format_info(FORMAT_INCHI)
    pathway = get_conversion_pathway(in_format, out_format)
    assert len(pathway) > 1

    run_with_arg_string(f"-l -f {in_format.id} -t {out_format.id}")
    captured = capsys.readouterr()

    _check_no_errors(captured, subtests)

    with subtests.test("Indicate no direct conversions possible"):
        assert compressed_match(f"No direct conversions are possible from {in_format.format_word()} to "
                                f"{out_format.format_word()}", captured.out)

    with subtests.test("Indicate chained conversion is possible"):
        assert compressed_match(f"A chained conversion is possible from {in_format.format_word()} to "
                                f"{out_format.format_word()} using registered converters:", captured.out)

    for i, step in enumerate(pathway):
        _check_step_details_present(step, captured, subtests, i=i)


def test_list_conversion_impossible(capsys, subtests):
    """Test that we get the expected output from listing when a chained conversion is not possible
    """

    in_format = "cif"
    out_format = "abinit"

    run_with_arg_string(f"-l -f {in_format} -t {out_format}")
    captured = capsys.readouterr()

    with subtests.test("correct error message"):
        assert compressed_match(f"No chained conversions are possible from {in_format} to {out_format}.", captured.out)

    # Check that igraph's warning is suppressed
    assert not compressed_match("Couldn't reach some vertices", captured.out)

    _check_no_errors(captured, subtests)


@pytest.mark.parametrize("include", ("one", "best", "shortest"))
def test_list_possible_chains(include, capsys, subtests):
    """Test that we can get lowest-weight chains with `--lp`"""

    in_format = FORMAT_INCHI
    out_format = FORMAT_MOLDY
    if include == "one":
        l_pathways = [get_conversion_pathway(in_format, out_format)]
    else:
        l_pathways = get_possible_conversion_pathways(in_format, out_format, include=include)

    run_with_arg_string(f"--lp {include} -f {in_format} -t {out_format}")
    captured = capsys.readouterr()

    with subtests.test("Explanatory test present"):
        assert compressed_match("Conversion pathways are listed below", captured.out)

    _check_no_errors(captured, subtests)

    for i, pathway in enumerate(l_pathways):
        for j, step in enumerate(pathway):
            _check_step_details_present(step, captured, subtests, include=include, i=i, j=j)


def test_list_chain_impossible(capsys, subtests):
    """Test that we get the expected output from listing chains when no chained conversion is possible
    """

    in_format = "cif"
    out_format = "abinit"

    run_with_arg_string(f"--lp -f {in_format} -t {out_format}")
    captured = capsys.readouterr()

    with subtests.test("correct error message"):
        assert compressed_match(f"No conversion pathway is possible from {in_format} to {out_format}.", captured.out)

    # Check that igraph's warning is suppressed
    assert not compressed_match("Couldn't reach some vertices", captured.out)

    _check_no_errors(captured, subtests)


def test_conversion_info_open_babel(capsys, subtests):
    """Test that we get the expected information on the 'Open Babel' converter
    """

    converter_name = const.CONVERTER_OB

    in_format = "xyz-1"
    out_format = "inchi"
    qual = get_conversion_quality(converter_name, in_format, out_format)

    # Test a basic listing of arguments, checking with the converter name in lowercase to be sure that works
    run_with_arg_string(f"-l {converter_name.lower()} -f {in_format} -t {out_format}")
    captured = capsys.readouterr()

    _check_no_errors(captured, subtests)

    # Check that conversion quality details are in the output as expected
    with subtests.test("conversion listed as possible"):
        assert compressed_match(f"Conversion from {in_format} to {out_format} with {converter_name} is "
                                f"possible with {qual.qual_str} conversion quality", captured.out)
    with subtests.test("data loss warning displayed"):
        assert compressed_match("WARNING: Potential data loss or extrapolation issues with this conversion:",
                                captured.out)
    with subtests.test("correct 2D quality note"):
        assert compressed_match(const.QUAL_NOTE_OUT_MISSING.format(const.QUAL_2D_LABEL), captured.out)
    with subtests.test("correct 3D quality note"):
        assert compressed_match(const.QUAL_NOTE_OUT_MISSING.format(const.QUAL_3D_LABEL), captured.out)
    with subtests.test("correct connections quality note"):
        assert compressed_match(const.QUAL_NOTE_IN_MISSING.format(const.QUAL_CONN_LABEL), captured.out)

    l_in_flags, l_in_options = get_in_format_args(converter_name, in_format)
    l_out_flags, l_out_options = get_out_format_args(converter_name, out_format)

    # Check headings for input/output flags/options are present if and only if some of those flags/options exist
    assert bool(l_in_flags) == compressed_match(f"Allowed input flags for format '{in_format}'", captured.out)
    assert bool(l_out_flags) == compressed_match(f"Allowed output flags for format '{out_format}'", captured.out)
    assert bool(l_in_options) == compressed_match(f"Allowed input options for format '{in_format}'", captured.out)
    assert bool(l_out_options) == compressed_match(f"Allowed output options for format '{out_format}'", captured.out)

    # Check that info for each flag and option is printed as expected
    for flag_info in l_in_flags + l_out_flags:
        info = flag_info.info if flag_info.info and flag_info.info != "N/A" else ""
        assert compressed_match(f"{flag_info.name}{flag_info.description}{info}", captured.out)
    for option_info in l_in_options + l_out_options:
        info = option_info.info if option_info.info and option_info.info != "N/A" else ""
        assert compressed_match(f"{option_info.name}<{option_info.brief}>{option_info.description}{info}",
                                captured.out)


@pytest.mark.parametrize("converter_name", [const.CONVERTER_C2X, const.CONVERTER_ATO])
def test_conversion_info_others(capsys, converter_name, subtests):
    """Test that we get the expected information on other converters
    """

    in_format = "pdb-0"
    out_format = "cif"
    qual = get_conversion_quality(converter_name, in_format, out_format)

    run_with_arg_string(f"-l {converter_name} -f {in_format} -t {out_format}")

    captured = capsys.readouterr()
    _check_no_errors(captured, subtests)

    # Check that conversion quality details are in the output as expected
    with subtests.test("conversion listed as possible"):
        assert compressed_match(f"Conversion from {in_format} to {out_format} with {converter_name} is "
                                f"possible with {qual.qual_str} conversion quality", captured.out)
    with subtests.test("data loss warning displayed"):
        assert compressed_match("WARNING: Potential data loss or extrapolation issues with this conversion:",
                                captured.out)
    with subtests.test("correct connections quality note"):
        assert compressed_match(const.QUAL_NOTE_OUT_MISSING.format(const.QUAL_CONN_LABEL), captured.out)


def test_format_info(capsys, subtests):
    """Test that we can successfully get information on a file format"""

    # Try to get info on an unambiguous format

    in_format = "inchi"
    in_format_info = get_format_info(in_format)
    run_with_arg_string(f"-l -f {in_format}")

    captured = capsys.readouterr()

    _check_no_errors(captured, subtests)

    # Check for basic format information
    with subtests.test("basic format information is present"):
        assert compressed_match(f"{in_format_info.disambiguated_name} (ID {in_format_info.id}): " +
                                in_format_info.description, captured.out)

    # Check for property information
    for attr, label in D_FORMAT_PROPERTY_ATTRS.items():
        with subtests.test(f"{attr} support is correct"):
            support_status = getattr(in_format_info, attr)
            if support_status:
                assert compressed_match(label + " supported", captured.out)
            elif support_status is False:
                assert compressed_match(label + " not supported", captured.out)
            else:
                assert compressed_match(label + " unknown whether or not to be supported", captured.out)


@pytest.mark.parametrize("in_format", [FORMAT_PDB_0, FORMAT_MOL])
def test_format_aliases(in_format, capsys, subtests):
    """Test that when a format has aliases, they're printed in the info about it"""

    in_format_info = get_format_info(in_format)
    run_with_arg_string(f"-l -f {in_format}")

    captured = capsys.readouterr()

    _check_no_errors(captured, subtests)

    l_other_alias_ids = [x for x in in_format_info.format_common_info.d_alias_exts.keys()
                         if x != in_format_info.id]
    for other_alias_id in l_other_alias_ids:
        with subtests.test("info for alias is present", id=other_alias_id):
            assert compressed_match(get_format_info(other_alias_id).format_inline(), captured.out)


def test_format_info_ambiguous(capsys, subtests):
    """Test that we get expected information for an ambiguous format"""

    out_format = "pdb"
    l_out_format_info = get_format_info(out_format, which="all")
    run_with_arg_string(f"-l -t {out_format}")

    captured = capsys.readouterr()

    _check_no_errors(captured, subtests)

    with subtests.test("ambiguous warning is present"):
        assert compressed_match(f"WARNING: Format '{out_format}' is ambiguous", captured.out)

    for out_format_info in l_out_format_info:
        with subtests.test("each possible match is present", id=out_format_info.id):
            assert compressed_match(out_format_info.format_oneline(), captured.out)


def test_format_info_in_unrecognised(capsys, subtests):
    """Test we get expected errors for unrecognised input format"""

    in_format = 99999
    with subtests.test("CLI exit on unrecognised format"), pytest.raises(SystemExit):
        run_with_arg_string(f"-l -f {in_format}")

    with subtests.test("unrecognised error present"):
        assert compressed_match(f"ERROR: Format '{in_format}' not recognised", capsys.readouterr().err)


def test_format_info_out_unrecognised(capsys, subtests):
    """Test we get expected errors for unrecognised output format"""

    out_format = "not_a_format"

    with subtests.test("CLI exit on unrecognised format"), pytest.raises(SystemExit):
        run_with_arg_string(f"-l -t {out_format}")

    with subtests.test("unrecognised error present"):
        assert compressed_match(f"ERROR: Format '{out_format}' not recognised", capsys.readouterr().err)
