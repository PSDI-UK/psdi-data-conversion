#!/usr/bin/env python3

"""@file psdi_data_conversion/main.py

Created 2025-01-14 by Bryan Gillis.

Entry-point file for the command-line interface for data conversion.
"""

import logging
import os
import sys
import textwrap
from argparse import ArgumentParser
from itertools import product

from psdi_data_conversion import constants as const
from psdi_data_conversion.constants import CL_SCRIPT_NAME, CONVERTER_DEFAULT, TERM_WIDTH
from psdi_data_conversion.converter import (D_CONVERTER_ARGS, L_REGISTERED_CONVERTERS, L_SUPPORTED_CONVERTERS,
                                            converter_is_registered, converter_is_supported,
                                            get_registered_converter_class, get_supported_converter_class,
                                            run_converter)
from psdi_data_conversion.converters.base import (FileConverterAbortException, FileConverterException,
                                                  FileConverterInputException)
from psdi_data_conversion.database import (FormatInfo, get_conversion_quality, get_converter_info, get_format_info,
                                           get_in_format_args, get_out_format_args, get_possible_conversions,
                                           get_possible_formats)
from psdi_data_conversion.file_io import split_archive_ext
from psdi_data_conversion.log_utility import get_log_level_from_str
from psdi_data_conversion.utils import regularize_name


def print_wrap(s: str, newline=False, err=False, **kwargs):
    """Print a string wrapped to the terminal width
    """
    if err:
        file = sys.stderr
    else:
        file = sys.stdout
    for line in s.split("\n"):
        print(textwrap.fill(line, width=TERM_WIDTH, **kwargs), file=file)
    if newline:
        print("")


class ConvertArgs:
    """Class storing arguments for data conversion, processed and determined from the input arguments.
    """

    def __init__(self, args):

        # Start by copying over arguments. Some share names with reserved words, so we have to use `getattr` for them

        # Positional arguments
        self.l_args: list[str] = args.l_args

        # Keyword arguments for standard conversion
        self.from_format: str | None = getattr(args, "from")
        self._input_dir: str | None = getattr(args, "in")
        self.to_format: str | None = args.to
        self._output_dir: str | None = args.out
        converter_name = getattr(args, "with")
        if isinstance(converter_name, str):
            self.name = regularize_name(converter_name)
        elif converter_name:
            self.name = regularize_name(" ".join(converter_name))
        else:
            self.name = None
        self.delete_input = args.delete_input
        self.from_flags: str = args.from_flags.replace(r"\-", "-")
        self.to_flags: str = args.to_flags.replace(r"\-", "-")
        self.from_options: str = args.from_options.replace(r"\-", "-")
        self.to_options: str = args.to_options.replace(r"\-", "-")
        self.no_check: bool = args.nc
        self.strict: bool = args.strict

        # Keyword arguments for alternative functionality
        self.list: bool = args.list

        # Logging/stdout arguments
        self.log_mode: bool = args.log_mode
        self.quiet = args.quiet
        self._log_file: str | None = args.log_file

        if not args.log_level:
            self.log_level = None
        else:
            try:
                self.log_level = get_log_level_from_str(args.log_level)
            except ValueError as e:
                # A ValueError indicates an unrecognised logging level, so we reraise this with the help flag to
                # indicate we want to provide this as feedback to the user so they can correct their command
                raise FileConverterInputException(str(e), help=True)

        # If formats were provided as ints, convert them to the int type now
        try:
            if self.from_format:
                self.from_format = int(self.from_format)
        except ValueError:
            pass
        try:
            if self.to_format:
                self.to_format = int(self.to_format)
        except ValueError:
            pass

        # Special handling for listing converters
        if self.list:
            # Force log mode to stdout and turn off quiet
            self.log_mode = const.LOG_STDOUT
            self.quiet = False

            # Get the converter name from the arguments if it wasn't provided by -w/--with
            if not self.name:
                self.name = regularize_name(" ".join(self.l_args))

            # For this operation, any other arguments can be ignored
            return

        # If not listing and a converter name wasn't supplied, use the default converter
        if not self.name:
            self.name = regularize_name(CONVERTER_DEFAULT)

        # Quiet mode is equivalent to logging mode == LOGGING_NONE, so normalize them if either is set
        if self.quiet:
            self.log_mode = const.LOG_NONE
        elif self.log_mode == const.LOG_NONE:
            self.quiet = True

        # Check validity of input

        if len(self.l_args) == 0:
            raise FileConverterInputException("One or more names of files to convert must be provided", help=True)

        if self._input_dir is not None and not os.path.isdir(self._input_dir):
            raise FileConverterInputException(f"The provided input directory '{self._input_dir}' does not exist as a "
                                              "directory", help=True)

        if self.to_format is None:
            msg = textwrap.fill("ERROR Output format (-t or --to) must be provided. For information on supported "
                                "formats and converters, call:\n")
            msg += f"{CL_SCRIPT_NAME} -l"
            raise FileConverterInputException(msg, msg_preformatted=True, help=True)

        # If the output directory doesn't exist, silently create it
        if self._output_dir is not None and not os.path.isdir(self._output_dir):
            if os.path.exists(self._output_dir):
                raise FileConverterInputException(f"Output directory '{self._output_dir}' exists but is not a "
                                                  "directory", help=True)
            os.makedirs(self._output_dir, exist_ok=True)

        # Check the converter is recognized
        if not converter_is_supported(self.name):
            msg = textwrap.fill(f"ERROR: Converter '{self.name}' not recognised", width=TERM_WIDTH)
            msg += f"\n\n{get_supported_converters()}"
            raise FileConverterInputException(msg, help=True, msg_preformatted=True)
        elif not converter_is_registered(self.name):
            converter_name = get_supported_converter_class(self.name).name
            msg = textwrap.fill(f"ERROR: Converter '{converter_name}' is not registered. It may be possible to "
                                "register it by installing an appropriate binary for your platform.", width=TERM_WIDTH)
            msg += f"\n\n{get_supported_converters()}"
            raise FileConverterInputException(msg, help=True, msg_preformatted=True)

        # Logging mode is valid
        if self.log_mode not in const.L_ALLOWED_LOG_MODES:
            raise FileConverterInputException(f"Unrecognised logging mode: {self.log_mode}. Allowed "
                                              f"modes are: {const.L_ALLOWED_LOG_MODES}", help=True)

        # Arguments specific to this converter
        self.d_converter_args = {}
        l_converter_args = D_CONVERTER_ARGS[self.name]
        if not l_converter_args:
            l_converter_args = []
        for arg_name, _, get_data in l_converter_args:
            # Convert the argument name to how it will be represented in the parsed_args object
            while arg_name.startswith("-"):
                arg_name = arg_name[1:]
            arg_name = arg_name.replace("-", "_")
            self.d_converter_args.update(get_data(getattr(args, arg_name)))

    @property
    def input_dir(self):
        """If the input directory isn't provided, use the current directory.
        """
        if self._input_dir is None:
            self._input_dir = os.getcwd()
        return self._input_dir

    @property
    def output_dir(self):
        """If the output directory isn't provided, use the input directory.
        """
        if self._output_dir is None:
            self._output_dir = self.input_dir
        return self._output_dir

    @property
    def log_file(self):
        """Determine a name for the log file if one is not provided.
        """
        if self._log_file is None:
            if self.list:
                self._log_file = const.DEFAULT_LISTING_LOG_FILE
            else:
                first_filename = os.path.join(self.input_dir, self.l_args[0])

                # Find the path to this file
                if not os.path.isfile(first_filename):
                    if self.from_format:
                        test_filename = first_filename + f".{self.from_format}"
                        if os.path.isfile(test_filename):
                            first_filename = test_filename
                        else:
                            raise FileConverterInputException(f"Input file {first_filename} cannot be found. Also "
                                                              f"checked for {test_filename}.", help=True)
                    else:
                        raise FileConverterInputException(f"Input file {first_filename} cannot be found.", help=True)

                filename_base = os.path.split(split_archive_ext(first_filename)[0])[1]
                if self.log_mode == const.LOG_FULL:
                    # For server-style logging, other files will be created and used for logs
                    self._log_file = None
                else:
                    self._log_file = os.path.join(self.output_dir, filename_base + const.LOG_EXT)
        return self._log_file


def get_argument_parser():
    """Get an argument parser for this script.

    Returns
    -------
    parser : ArgumentParser
        An argument parser set up with the allowed command-line arguments for this script.
    """

    parser = ArgumentParser()

    # Positional arguments
    parser.add_argument("l_args", type=str, nargs="*",
                        help="Normally, file(s) to be converted or zip/tar archives thereof. If an archive or archives "
                        "are provided, the output will be packed into an archive of the same type. Filenames should be "
                        "provided as either relative to the input directory (default current directory) or absolute. "
                        "If the '-l' or '--list' flag is set, instead the name of a converter can be used here to get "
                        "information on it.")

    # Keyword arguments for standard conversion
    parser.add_argument("-f", "--from", type=str, default=None,
                        help="The input (convert from) file extension (e.g., smi). If not provided, will attempt to "
                             "auto-detect format.")
    parser.add_argument("-i", "--in", type=str, default=None,
                        help="The directory containing the input file(s), default current directory.")
    parser.add_argument("-t", "--to", type=str, default=None,
                        help="The output (convert to) file extension (e.g., cmi).")
    parser.add_argument("-o", "--out", type=str, default=None,
                        help="The directory where output files should be created. If not provided, output files will "
                        "be created in -i/--in directory if that was provided, or else in the directory containing the "
                        "first input file.")
    parser.add_argument("-w", "--with", type=str, nargs="+",
                        help="The converter to be used (default 'Open Babel').")
    parser.add_argument("--delete-input", action="store_true",
                        help="If set, input files will be deleted after conversion, default they will be kept")
    parser.add_argument("--from-flags", type=str, default="",
                        help="Any command-line flags to be provided to the converter for reading in the input file(s). "
                             "For information on the flags accepted by a converter and its required format for them, "
                             "call this script with '-l <converter name>'. If the set of flags includes any spaces, it "
                             "must be quoted, and if hyphens are used, the first preceding hyphen for each flag must "
                             "be backslash-escaped, e.g. '--from-flags \"\\-a \\-bc \\--example\"'")
    parser.add_argument("--to-flags", type=str, default="",
                        help="Any command-line flags to be provided to the converter for writing the output file(s). "
                             "For information on the flags accepted by a converter and its required format for them, "
                             "call this script with '-l <converter name>'. If the set of flags includes any spaces, it "
                             "must be quoted, and if hyphens are used, the first preceding hyphen for each flag must "
                             "be backslash-escaped, e.g. '--to-flags \"\\-a \\-bc \\--example\"'")
    parser.add_argument("--from-options", type=str, default="",
                        help="Any command-line options to be provided to the converter for reading in the input "
                             "file(s). For information on the options accepted by a converter and its required format "
                             "for them, call this script with '-l <converter name>'. If the set of options includes "
                             "any spaces, it must be quoted, and the first preceding hyphen for each option must be "
                             "backslash-escaped, e.g. '--from-options \"\\-x xval --opt optval\"'")
    parser.add_argument("--to-options", type=str, default="",
                        help="Any command-line options to be provided to the converter for writing the output "
                             "file(s). For information on the options accepted by a converter and its required format "
                             "for them, call this script with '-l <converter name>'. If the set of options includes "
                             "any spaces, it must be quoted, and the first preceding hyphen for each option must be "
                             "backslash-escaped, e.g. '--to-options \"\\-x xval --opt optval\"'")
    parser.add_argument("-s", "--strict", action="store_true",
                        help="If set, will fail if one of the input files has the wrong extension (including those "
                        "contained in archives, but not the archive files themselves). Otherwise, will only print a "
                        "warning in this case.")
    parser.add_argument("--nc", "--no-check", action="store_true",
                        help="If set, will not perform a pre-check in the database on the validity of a conversion. "
                        "Setting this will result in a less human-friendly error message (or may even falsely indicate "
                        "success) if the conversion is not supported, but will save some execution time. Recommended "
                        "only for automated execution after the user has confirmed a conversion is supported")

    # Keyword arguments specific to converters
    for converter_name in L_REGISTERED_CONVERTERS:
        l_converter_args = D_CONVERTER_ARGS[converter_name]
        if l_converter_args:
            for arg_name, kwargs, _ in l_converter_args:
                parser.add_argument(arg_name, **kwargs)

    # Keyword arguments for alternative functionality
    parser.add_argument("-l", "--list", action="store_true",
                        help="If provided alone, lists all available converters. If the name of a converter is "
                             "provided, gives information on the converter and any command-line flags it accepts.")

    # Logging/stdout arguments
    parser.add_argument("-g", "--log-file", type=str, default=None,
                        help="The name of the file to log to. This can be provided relative to the current directory "
                        "(e.g. '-g ../logs/log-file.txt') or fully qualified (e.g. /path/to/log-file.txt). "
                        "If not provided, the log file will be named after the =first input file (+'.log') and placed "
                        "in the output directory (specified with -o/--out).\n"
                        "In 'full' logging mode (not recommended with this interface), this will apply only to logs "
                        "from the outermost level of the script if explicitly specified. If not explicitly specified, "
                        "those logs will be sent to stderr.")
    parser.add_argument("--log-mode", type=str, default=const.LOG_SIMPLE,
                        help="How logs should be stored. Allowed values are: \n"
                        "- 'full' - Multi-file logging, not recommended for the CLI, but allowed for a compatible "
                        "interface with the public web app"
                        "- 'simple' - Logs saved to one file"
                        "- 'stdout' - Output logs and errors only to stdout"
                        "- 'none' - Output only errors to stdout")
    parser.add_argument("-q", "--quiet", action="store_true",
                        help="If set, all terminal output aside from errors will be suppressed and no log file will be "
                             "generated.")
    parser.add_argument("--log-level", type=str, default=None,
                        help="The desired level to log at. Allowed values are: 'DEBUG', 'INFO', 'WARNING', 'ERROR, "
                             "'CRITICAL'. Default: 'INFO' for logging to file, 'WARNING' for logging to stdout")

    return parser


def parse_args():
    """Parses arguments for this executable.

    Returns
    -------
    args : Namespace
        The parsed arguments.
    """

    parser = get_argument_parser()

    args = ConvertArgs(parser.parse_args())

    return args


def detail_converter_use(args: ConvertArgs):
    """Prints output providing information on a specific converter, including the flags and options it allows
    """

    converter_info = get_converter_info(args.name)
    converter_class = get_supported_converter_class(args.name)
    converter_name = converter_class.name

    print_wrap(f"{converter_name}: {converter_info.description} ({converter_info.url})", break_long_words=False,
               break_on_hyphens=False, newline=True)

    if converter_class.info:
        print_wrap(converter_class.info, break_long_words=False, break_on_hyphens=False, newline=True)

    # If both an input and output format are specified, provide the degree of success for this conversion. Otherwise
    # list possible input/output formats
    if args.from_format is not None and args.to_format is not None:
        qual = get_conversion_quality(args.name, args.from_format, args.to_format)
        if qual is None:
            print_wrap(f"Conversion from '{args.from_format}' to '{args.to_format}' with {converter_name} is not "
                       "supported.", newline=True)
        else:
            print_wrap(f"Conversion from '{args.from_format}' to '{args.to_format}' with {converter_name} is "
                       f"possible with {qual.qual_str} conversion quality", newline=True)
            # If there are any potential issues with the conversion, print them out
            if qual.details:
                print_wrap("WARNING: Potential data loss or extrapolation issues with this conversion:")
                for detail_line in qual.details.split("\n"):
                    print_wrap(f"- {detail_line}")
                print("")
    else:
        l_input_formats, l_output_formats = get_possible_formats(args.name)

        # If one format was supplied, check if it's supported
        for (format_name, l_formats, to_or_from) in ((args.from_format, l_input_formats, "from"),
                                                     (args.to_format, l_output_formats, "to")):
            if format_name is None:
                continue
            if format_name in l_formats:
                optional_not: str = ""
            else:
                optional_not: str = "not "
            print_wrap(f"Conversion {to_or_from} {format_name} is {optional_not}supported by {converter_name}.\n")

        # List all possible formats, and which can be used for input and which for output
        s_all_formats: set[FormatInfo] = set(l_input_formats)
        s_all_formats.update(l_output_formats)
        l_all_formats: list[FormatInfo] = list(s_all_formats)
        l_all_formats.sort(key=lambda x: x.disambiguated_name.lower())

        print_wrap(f"File formats supported by {converter_name}:", newline=True)
        max_format_length = max([len(x.disambiguated_name) for x in l_all_formats])
        print(" "*(max_format_length+4) + "   INPUT  OUTPUT")
        print(" "*(max_format_length+4) + "   -----  ------")
        for file_format in l_all_formats:
            in_yes_or_no = "yes" if file_format in l_input_formats else "no"
            out_yes_or_no = "yes" if file_format in l_output_formats else "no"
            print(f"    {file_format.disambiguated_name:>{max_format_length}}{in_yes_or_no:>8}{out_yes_or_no:>8}")
        print("")

    if converter_class.allowed_flags is None:
        print_wrap("Information has not been provided about general flags accepted by this converter.", newline=True)
    elif len(converter_class.allowed_flags) > 0:
        print_wrap("Allowed general flags:")
        for flag, d_data, _ in converter_class.allowed_flags:
            help = d_data.get("help", "(No information provided)")
            print(f"  {flag}")
            print_wrap(help, width=TERM_WIDTH, initial_indent=" "*4, subsequent_indent=" "*4)
        print("")

    if converter_class.allowed_options is None:
        print_wrap("Information has not been provided about general options accepted by this converter.", newline=True)
    elif len(converter_class.allowed_options) > 0:
        print_wrap("Allowed general options:")
        for option, d_data, _ in converter_class.allowed_options:
            help = d_data.get("help", "(No information provided)")
            print(f"  {option} <val(s)>")
            print(textwrap.fill(help, initial_indent=" "*4, subsequent_indent=" "*4))
        print("")

    # If input/output-format specific flags or options are available for the converter but a format isn't available,
    # we'll want to take note of that and mention that at the end of the output
    mention_input_format = False
    mention_output_format = False

    if args.from_format is not None:
        from_format = args.from_format
        in_flags, in_options = get_in_format_args(args.name, from_format)
    else:
        in_flags, in_options = [], []
        from_format = "N/A"
        if converter_class.has_in_format_flags_or_options:
            mention_input_format = True

    if args.to_format is not None:
        to_format = args.to_format
        out_flags, out_options = get_out_format_args(args.name, to_format)
    else:
        out_flags, out_options = [], []
        to_format = "N/A"
        if converter_class.has_out_format_flags_or_options:
            mention_output_format = True

    # Number of character spaces allocated for flags/options when printing them out
    ARG_LEN = 20

    for l_args, flag_or_option, input_or_output, format_name in ((in_flags, "flag", "input", from_format),
                                                                 (in_options, "option", "input", from_format),
                                                                 (out_flags, "flag", "output", to_format),
                                                                 (out_options, "option", "output", to_format)):
        if len(l_args) == 0:
            continue
        print_wrap(f"Allowed {input_or_output} {flag_or_option}s for format '{format_name}':")
        for arg_info in l_args:
            if flag_or_option == "flag":
                optional_brief = ""
            else:
                optional_brief = f" <{arg_info.brief}>"
            print_wrap(f"{arg_info.flag+optional_brief:>{ARG_LEN}}  {arg_info.description}",
                       subsequent_indent=" "*(ARG_LEN+2))
            if arg_info.info and arg_info.info != "N/A":
                print_wrap(arg_info.info,
                           initial_indent=" "*(ARG_LEN+2),
                           subsequent_indent=" "*(ARG_LEN+2))
        print("")

    # Now at the end, bring up input/output-format-specific flags and options
    if mention_input_format and mention_output_format:
        print_wrap("For details on input/output flags and options allowed for specific formats, call:\n"
                   f"{CL_SCRIPT_NAME} -l {converter_name} -f <input_format> -t <output_format>")
    elif mention_input_format:
        print_wrap("For details on input flags and options allowed for a specific format, call:\n"
                   f"{CL_SCRIPT_NAME} -l {converter_name} -f <input_format> [-t <output_format>]")
    elif mention_output_format:
        print_wrap("For details on output flags and options allowed for a specific format, call:\n"
                   f"{CL_SCRIPT_NAME} -l {converter_name} -t <output_format> [-f <input_format>]")


def list_supported_formats(err=False):
    """Prints a list of all formats recognised by at least one registered converter
    """
    # Make a list of all formats recognised by at least one registered converter
    s_all_formats: set[FormatInfo] = set()
    s_registered_formats: set[FormatInfo] = set()
    for converter_name in L_SUPPORTED_CONVERTERS:
        l_in_formats, l_out_formats = get_possible_formats(converter_name)

        # To make sure we don't see any unexpected duplicates in the set due to cached/uncached values, get the
        # disambiguated name of each format first
        [x.disambiguated_name for x in l_in_formats]
        [x.disambiguated_name for x in l_out_formats]

        s_all_formats.update(l_in_formats)
        s_all_formats.update(l_out_formats)
        if converter_name in L_REGISTERED_CONVERTERS:
            s_registered_formats.update(l_in_formats)
            s_registered_formats.update(l_out_formats)

    s_unregistered_formats = s_all_formats.difference(s_registered_formats)

    # Convert the sets to lists and alphabetise them
    l_registered_formats = list(s_registered_formats)
    l_registered_formats.sort(key=lambda x: x.disambiguated_name.lower())
    l_unregistered_formats = list(s_unregistered_formats)
    l_unregistered_formats.sort(key=lambda x: x.disambiguated_name.lower())

    # Pad the format strings to all be the same length. To keep columns aligned, all padding is done with non-
    # breaking spaces (\xa0), and each format is followed by a single normal space
    longest_format_len = max([len(x.disambiguated_name) for x in l_registered_formats])
    l_padded_formats = [f"{x.disambiguated_name:\xa0<{longest_format_len}} " for x in l_registered_formats]

    print_wrap("Formats supported by registered converters: ", err=err, newline=True)
    print_wrap("".join(l_padded_formats), err=err, initial_indent="  ", subsequent_indent="  ", newline=True)

    if l_unregistered_formats:
        longest_unregistered_format_len = max([len(x) for x in l_unregistered_formats])
        l_padded_unregistered_formats = [f"{x:\xa0<{longest_unregistered_format_len}} "
                                         for x in l_unregistered_formats]
        print_wrap("Formats supported by unregistered converters which are supported by this package: ", err=err,
                   newline=True)
        print_wrap("".join(l_padded_unregistered_formats), err=err,
                   initial_indent="  ", subsequent_indent="  ", newline=True)

    print_wrap("Note that not all formats are supported with all converters, or both as input and as output.")
    if err:
        print("")
        print_wrap("For more details on a format, call:")
        print(f"{CL_SCRIPT_NAME} -l -f <format>")


def detail_format(format_name: str):
    """Prints details on a format
    """

    l_format_info = get_format_info(format_name, which="all")

    if len(l_format_info) == 0:
        print_wrap(f"ERROR: Format '{format_name}' not recognised", err=True, newline=True)
        list_supported_formats(err=True)
        exit(1)

    if len(l_format_info) > 1:
        print_wrap(f"WARNING: Format '{format_name}' is ambiguous and could refer to multiple formats. It may be "
                   "necessary to explicitly specify which you want to use when calling this script, e.g. with "
                   f"'-f {format_name}-0' - see the disambiguated names in the list below:", newline=True)

    for format_info in l_format_info:
        print_wrap(f"{format_info.id}: {format_info.disambiguated_name} ({format_info.note})")


def detail_formats_and_possible_converters(from_format: str, to_format: str):
    """Prints details on converters that can perform a conversion from one format to another
    """

    # Check that both formats are valid, and print an error if not
    either_format_failed = False

    try:
        get_format_info(from_format, which=0)
    except KeyError:
        either_format_failed = True
        print_wrap(f"ERROR: Input format '{from_format}' not recognised", newline=True, err=True)

    try:
        get_format_info(to_format, which=0)
    except KeyError:
        either_format_failed = True
        print_wrap(f"ERROR: Output format '{from_format}' not recognised", newline=True, err=True)

    if either_format_failed:
        # Let the user know about formats which are allowed
        list_supported_formats(err=True)
        exit(1)

    # Provide details on both the input and output formats
    detail_format(from_format)
    print()
    detail_format(to_format)

    l_possible_conversions = get_possible_conversions(from_format, to_format)

    # Get a list of all different formats which share the provided name, cutting out duplicates
    l_from_formats = list(set([x[1] for x in l_possible_conversions]))
    l_from_formats.sort(key=lambda x: x.disambiguated_name)
    l_to_formats = list(set([x[2] for x in l_possible_conversions]))
    l_to_formats.sort(key=lambda x: x.disambiguated_name)

    # Loop over all possible combinations of formats

    for possible_from_format, possible_to_format in product(l_from_formats, l_to_formats):

        from_name = possible_from_format.disambiguated_name
        to_name = possible_to_format.disambiguated_name

        l_conversions_matching_formats = [x for x in l_possible_conversions
                                          if x[1] == possible_from_format and x[2] == possible_to_format]

        l_possible_registered_converters = [get_registered_converter_class(x[0]).name
                                            for x in l_conversions_matching_formats
                                            if x[0] in L_REGISTERED_CONVERTERS]
        l_possible_unregistered_converters = [get_supported_converter_class(x[0]).name
                                              for x in l_conversions_matching_formats
                                              if x[0] in L_SUPPORTED_CONVERTERS and x[0] not in L_REGISTERED_CONVERTERS]

        print()

        if len(l_possible_registered_converters)+len(l_possible_unregistered_converters) == 0:
            print_wrap(f"No converters are available which can perform a conversion from {from_name} to "
                       f"{to_name}")
            continue
        elif len(l_possible_registered_converters) == 0:
            print_wrap(f"No registered converters can perform a conversion from {from_name} to "
                       f"{to_name}, however the following converters are supported by this package on other "
                       "platforms and can perform this conversion:", newline=True)
            print("\n    ".join(l_possible_unregistered_converters))
            continue

        print_wrap(f"The following registered converters can convert from {from_name} to "
                   f"{to_name}:", newline=True)
        print("    " + "\n    ".join(l_possible_registered_converters) + "\n")
        if l_possible_unregistered_converters:
            print("")
            print_wrap("Additionally, the following converters are supported by this package on other platforms and "
                       "can perform this conversion:", newline=True)
            print("    " + "\n    ".join(l_possible_unregistered_converters) + "\n")

        print_wrap("For details on input/output flags and options allowed by a converter for this conversion, call:")
        print(f"{CL_SCRIPT_NAME} -l <converter name> -f {from_name} -t {to_name}")


def get_supported_converters():
    """Gets a string containing a list of supported converters
    """

    MSG_NOT_REGISTERED = "(supported but not registered)"

    l_converters: list[str] = []
    any_not_registered = False
    for converter_name in L_SUPPORTED_CONVERTERS:
        converter_text = get_supported_converter_class(converter_name).name
        if converter_name not in L_REGISTERED_CONVERTERS:
            converter_text += f" {MSG_NOT_REGISTERED}"
            any_not_registered = True
        l_converters.append(converter_text)

    output_str = "Available converters: \n\n    " + "\n    ".join(l_converters)

    if any_not_registered:
        output_str += (f"\n\nConverters marked as \"{MSG_NOT_REGISTERED}\" are supported by this package, but no "
                       "appropriate binary for your platform was either distributed with this package or "
                       "found on your system")

    return output_str


def list_supported_converters(err=False):
    """Prints a list of supported converters for the user
    """
    if err:
        file = sys.stderr
    else:
        file = sys.stdout
    print(get_supported_converters() + "\n", file=file)


def detail_converters_and_formats(args: ConvertArgs):
    """Prints details on available converters and formats for the user.
    """
    if args.name in L_SUPPORTED_CONVERTERS:
        detail_converter_use(args)
        if args.name not in L_REGISTERED_CONVERTERS:
            print_wrap("WARNING: This converter is supported by this package but is not registered. It may be possible "
                       "to register it by installing an appropriate binary on your system.", err=True)
        return

    elif args.name != "":
        print_wrap(f"ERROR: Converter '{args.name}' not recognized.", err=True, newline=True)
        list_supported_converters(err=True)
        exit(1)
    elif args.from_format and args.to_format:
        detail_formats_and_possible_converters(args.from_format, args.to_format)
        return
    elif args.from_format:
        detail_format(args.from_format)
        return
    elif args.to_format:
        detail_format(args.to_format)
        return

    list_supported_converters()
    list_supported_formats()

    print("")

    print_wrap("For more details on a converter, call:")
    print(f"{CL_SCRIPT_NAME} -l <converter name>\n")

    print_wrap("For more details on a format, call:")
    print(f"{CL_SCRIPT_NAME} -l -f <format>\n")

    print_wrap("For a list of converters that can perform a desired conversion, call:")
    print(f"{CL_SCRIPT_NAME} -l -f <input format> -t <output format>\n")

    print_wrap("For a list of options provided by a converter for a desired conversion, call:")
    print(f"{CL_SCRIPT_NAME} -l <converter name> -f <input format> -t <output format>")


def run_from_args(args: ConvertArgs):
    """Workhorse function to perform primary execution of this script, using the provided parsed arguments.

    Parameters
    ----------
    args : ConvertArgs
        The parsed arguments for this script.
    """

    # Check if we've been asked to list options
    if args.list:
        return detail_converters_and_formats(args)

    data = {'success': 'unknown',
            'from_flags': args.from_flags,
            'to_flags': args.to_flags,
            'from_options': args.from_options,
            'to_options': args.to_options,
            'from_arg_flags': '',
            'from_args': '',
            'to_arg_flags': '',
            'to_args': '',
            'upload_file': 'true'}
    data.update(args.d_converter_args)

    success = True

    for filename in args.l_args:

        # Search for the file in the input directory
        qualified_filename = os.path.join(args.input_dir, filename)
        if not os.path.isfile(qualified_filename):
            # Check if we can add the format to it as an extension to find it
            ex_extension = f".{args.from_format}"
            if not qualified_filename.endswith(ex_extension):
                qualified_filename += ex_extension
                if not os.path.isfile(qualified_filename):
                    print_wrap(f"ERROR: Cannot find file {filename+ex_extension} in directory {args.input_dir}",
                               err=True)
                    continue
            else:
                print_wrap(f"ERROR: Cannot find file {filename} in directory {args.input_dir}", err=True)
                continue

        if not args.quiet:
            print_wrap(f"Converting {filename} to {args.to_format}...", newline=True)

        try:
            conversion_result = run_converter(filename=qualified_filename,
                                              to_format=args.to_format,
                                              from_format=args.from_format,
                                              name=args.name,
                                              data=data,
                                              use_envvars=False,
                                              upload_dir=args.input_dir,
                                              download_dir=args.output_dir,
                                              no_check=args.no_check,
                                              strict=args.strict,
                                              log_file=args.log_file,
                                              log_mode=args.log_mode,
                                              log_level=args.log_level,
                                              delete_input=args.delete_input,
                                              refresh_local_log=False)
        except FileConverterAbortException as e:
            if not e.logged:
                print_wrap(f"ERROR: Attempt to convert file {filename} aborted with status code {e.status_code} and "
                           f"message:\n{e}\n", err=True)
                e.logged = True
            success = False
            continue
        except FileConverterException as e:
            if e.help and not e.logged:
                print_wrap(f"ERROR: {e}", err=True)
                e.logged = True
            elif "Conversion from" in str(e) and "is not supported" in str(e):
                if not e.logged:
                    print_wrap(f"ERROR: {e}", err=True, newline=True)
                detail_formats_and_possible_converters(args.from_format, args.to_format)
            elif e.help and not e.logged:
                if e.msg_preformatted:
                    print(e, file=sys.stderr)
                else:
                    print_wrap(f"ERROR: {e}", err=True)
            elif not e.logged:
                print_wrap(f"ERROR: Attempt to convert file {filename} failed at converter initialization with "
                           f"exception type {type(e)} and message: \n{e}\n", err=True)
            e.logged = True
            success = False
            continue
        except Exception as e:
            if not hasattr(e, "logged") or e.logged is False:
                print_wrap(f"ERROR: Attempt to convert file {filename} failed with exception type {type(e)} and "
                           f"message: \n{e}\n", err=True)
                e.logged = True
            success = False
            continue

        if not args.quiet:
            print_wrap("Success! The converted file can be found at:",)
            print(f"  {conversion_result.output_filename}")
            print_wrap("The log can be found at:")
            print(f"  {conversion_result.log_filename}")

    if not success:
        exit(1)


def main():
    """Standard entry-point function for this script.
    """

    # If no inputs were provided, print a message about usage
    if len(sys.argv) == 1:
        print_wrap("See the README.md file for information on using this utility and examples of basic usage, or for "
                   "detailed explanation of arguments call:")
        print(f"{CL_SCRIPT_NAME} -h")
        exit(1)

    try:
        args = parse_args()
    except FileConverterInputException as e:
        if not e.help:
            raise
        # If we get an exception with the help flag set, it's likely due to user error, so don't bother them with a
        # traceback and simply print the message to stderr
        if e.msg_preformatted:
            print(e, file=sys.stderr)
        else:
            print_wrap(f"ERROR: {e}", err=True)
        exit(1)

    if (args.log_mode == const.LOG_SIMPLE or args.log_mode == const.LOG_FULL) and args.log_file:
        # Delete any previous local log if it exists
        try:
            os.remove(args.log_file)
        except FileNotFoundError:
            pass
        logging.basicConfig(filename=args.log_file, level=args.log_level,
                            format=const.LOG_FORMAT, datefmt=const.TIMESTAMP_FORMAT)
    else:
        logging.basicConfig(level=args.log_level, format=const.LOG_FORMAT)

    logging.debug("#")
    logging.debug("# Beginning execution of script `%s`", __file__)
    logging.debug("#")

    run_from_args(args)

    logging.debug("#")
    logging.debug("# Finished execution of script `%s`", __file__)
    logging.debug("#")


if __name__ == "__main__":

    main()
