#!/usr/bin/env python3

"""@file psdi_data_conversion/main.py

Created 2025-01-14 by Bryan Gillis.

Entry-point file for the command-line interface for data conversion.
"""

import logging
from argparse import ArgumentParser
import os
import sys
import textwrap

from psdi_data_conversion import constants as const
from psdi_data_conversion.constants import ARG_LEN, TERM_WIDTH
from psdi_data_conversion.converter import D_REGISTERED_CONVERTERS, L_REGISTERED_CONVERTERS, run_converter
from psdi_data_conversion.converters.base import FileConverterAbortException, FileConverterInputException
from psdi_data_conversion.database import (get_converter_info, get_degree_of_success, get_in_format_args,
                                           get_out_format_args, get_possible_converters, get_possible_formats)


class FileConverterHelpException(FileConverterInputException):
    """An exception class which indicates an error where we will likely want to help the user figure out how to
    correctly use the CLI instead of simply printing a traceback
    """
    pass


def print_wrap(s, newline=False, **kwargs):
    """Print a string wrapped to the terminal width
    """
    print(textwrap.fill(s, width=TERM_WIDTH, **kwargs))
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
        self._from_format: str | None = getattr(args, "from")
        self._input_dir: str | None = getattr(args, "in")
        self.to_format: str | None = args.to
        self._output_dir: str | None = args.out
        converter_name = getattr(args, "with")
        if isinstance(converter_name, str):
            self.name = converter_name
        else:
            self.name: str = " ".join(converter_name)
        self.delete_input = args.delete_input
        self.from_flags: str = args.from_flags.replace(r"\-", "-")
        self.to_flags: str = args.to_flags.replace(r"\-", "-")

        # Keyword arguments specific to OpenBabel conversion
        self.coord_gen: str
        if args.coord_gen is None:
            self.coord_gen = const.DEFAULT_COORD_GEN
        else:
            self.coord_gen = args.coord_gen[0]

        self.coord_gen_qual: str
        if args.coord_gen is None or len(args.coord_gen) == 1:
            self.coord_gen_qual = const.DEFAULT_COORD_GEN_QUAL
        else:
            self.coord_gen_qual = args.coord_gen[1]

        # Keyword arguments for alternative functionality
        self.list: bool = args.list

        # Logging/stdout arguments
        self.log_mode: bool = args.log_mode
        self.quiet = args.quiet
        self._log_file: str | None = args.log_file

        if args.log_level.lower() == "debug":
            self.log_level = logging.DEBUG
        elif args.log_level.lower() == "info":
            self.log_level = logging.INFO
        elif args.log_level.lower().startswith("warn"):
            self.log_level = logging.WARNING
        elif args.log_level.lower() == "error":
            self.log_level = logging.ERROR
        elif args.log_level.lower() == "critical":
            self.log_level = logging.CRITICAL
        elif args.log_level:
            raise FileConverterHelpException(f"Unrecognised logging level: {args.log_level}")

        # Special handling for listing converters
        if self.list:
            # Force log mode to stdout and turn off quiet
            self.log_mode = const.LOG_STDOUT
            self.quiet = False

            # Get the converter name from the arguments
            self.name = " ".join(self.l_args)

            # For this operation, any other arguments can be ignored
            return

        # Quiet mode is equivalent to logging mode == LOGGING_NONE, so normalize them if either is set
        if self.quiet:
            self.log_mode = const.LOG_NONE
        elif self.log_mode == const.LOG_NONE:
            self.quiet = True

        # Check validity of input

        if len(self.l_args) == 0:
            raise FileConverterHelpException("One or more names of files to convert must be provided")

        if self._input_dir is not None and not os.path.isdir(self._input_dir):
            raise FileConverterHelpException(f"The provided input directory '{self._input_dir}' does not exist as a "
                                             "directory")

        if self.to_format is None:
            raise FileConverterHelpException("Output format (-t or --to) must be provided")

        # If the output directory doesn't exist, silently create it
        if self._output_dir is not None and not os.path.isdir(self._output_dir):
            if os.path.exists(self._output_dir):
                raise FileConverterHelpException(
                    f"Output directory '{self._output_dir}' exists but is not a directory")
            os.makedirs(self._output_dir, exist_ok=True)

        # Check the converter is recognized
        if self.name not in L_REGISTERED_CONVERTERS:
            raise FileConverterHelpException(f"Converter '{self.name}' not recognised")

        # No more than two arguments supplied to --coord-gen
        if args.coord_gen is not None and len(args.coord_gen) > 2:
            raise FileConverterHelpException("At most two arguments may be provided to --coord-gen, the mode and "
                                             "quality, e.g. '--coord-gen Gen3D best'")

        # Coordinate generation options are valid
        if self.coord_gen not in const.L_ALLOWED_COORD_GENS:
            raise FileConverterHelpException(f"Coordinate generation type '{self.coord_gen}' not recognised. Allowed "
                                             f"types are: {const.L_ALLOWED_COORD_GENS}")
        if self.coord_gen_qual not in const.L_ALLOWED_COORD_GEN_QUALS:
            raise FileConverterHelpException(f"Coordinate generation quality '{self.coord_gen_qual}' not recognised. "
                                             f"Allowed qualities are: {const.L_ALLOWED_COORD_GEN_QUALS}")

        # Logging mode is valid
        if self.log_mode not in const.L_ALLOWED_LOG_MODES:
            raise FileConverterHelpException(f"Unrecognised logging mode: {self.log_mode}. Allowed "
                                             f"modes are: {const.L_ALLOWED_LOG_MODES}")

    @property
    def from_format(self):
        """If the input file format isn't provided, determine it from the first file in the list.
        """
        if self._from_format is None:
            if self.list:
                # from_format isn't required in list mode, so don't raise an exception
                return None
            first_filename = self.l_args[0]
            ext = os.path.splitext(first_filename)[1]
            if len(ext) == 0:
                raise FileConverterHelpException("Input file format (-f or --from) was not provided, and cannot "
                                                 f"determine it automatically from filename '{first_filename}'")
            # Format will be the extension, minus the leading period
            self._from_format = ext[1:]
        return self._from_format

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
                    test_filename = first_filename + f".{self.from_format}"
                    if os.path.isfile(test_filename):
                        first_filename = test_filename
                    else:
                        raise FileConverterHelpException(f"Input file {first_filename} cannot be found. Also "
                                                         f"checked for {test_filename}.")

                filename_base = os.path.split(os.path.splitext(first_filename)[0])[1]
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
                        help="Normally, file(s) to be converted. Filenames should be provided as either relative to "
                             "the input directory (default current directory) or absolute. If the '-l' or '--list' "
                             "flag is set, instead the name of a converter can be used here to get information on it.")

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
    parser.add_argument("-w", "--with", type=str, nargs="+", default="Open Babel",
                        help="The converter to be used (default 'Open Babel').")
    parser.add_argument("--delete-input", action="store_true",
                        help="If set, input files will be deleted after conversion, default they will be kept")
    parser.add_argument("--from-flags", type=str, default="",
                        help="Any command-line flags to be provided to the converter for reading in the input file(s). "
                             "For information on the flags accepted by a converter, call this script with '-l "
                             "<converter name>'. The first preceding hyphen for each flag must be backslash-escaped, "
                             "e.g. '--from-flags \"\\-a \\-bc \\--example\"'")
    parser.add_argument("--to-flags", type=str, default="",
                        help="Any command-line flags to be provided to the converter for writing the output file(s). "
                             "For information on the flags accepted by a converter, call this script with '-l "
                             "<converter name>'. The first preceding hyphen for each flag must be backslash-escaped, "
                             "e.g. '--to-flags \"\\-a \\-bc \\--example\"'")

    # Keyword arguments specific to OpenBabel conversion
    parser.add_argument("--coord-gen", type=str, default=None, nargs="+",
                        help="(Open Babel converter only). The mode to be used for Open Babel calculation of atomic "
                             "coordinates, and optionally the quality of the conversion. The mode should be one of "
                             "'Gen2D', 'Gen3D', or 'neither' (default 'neither'). The quality, if supplied, should be "
                             "one of 'fastest', 'fast', 'medium', 'better' or 'best' (default 'medium'). E.g. "
                             "'--coord-gen Gen2D' (quality defaults to 'medium'), '--coord-gen Gen3D best'")

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
    parser.add_argument("--log-level", type=str, default="WARNING",
                        help="The desired level to log at. Allowed values are: 'DEBUG', 'INFO', 'WARNING', 'ERROR, "
                             "'CRITICAL'. Default: 'INFO'")

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
    converter_class = D_REGISTERED_CONVERTERS[args.name]

    print_wrap(f"{converter_info.name}: {converter_info.description} ({converter_info.url})", break_long_words=False,
               break_on_hyphens=False, newline=True)

    # If both an input and output format are specified, provide the degree of success for this conversion. Otherwise
    # list possible input output formats
    if args.from_format is not None and args.to_format is not None:
        dos = get_degree_of_success(args.name, args.from_format, args.to_format)
        if dos is None:
            print_wrap(f"Conversion from '{args.from_format}' to '{args.to_format}' with {args.name} is not "
                       "supported.", newline=True)
        else:
            print_wrap(f"Conversion from '{args.from_format}' to '{args.to_format}' with {args.name} is "
                       f"possible with the following note on degree of success: {dos}", newline=True)
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
                optional_not: str = " not"
            print_wrap(f"Conversion {to_or_from} {format_name} is {optional_not}supported by {args.name}.\n")

        # List all possible formats, and which can be used for input and which for output
        s_all_formats: set[str] = set(l_input_formats)
        s_all_formats.update(l_output_formats)
        l_all_formats: list[str] = list(s_all_formats)
        l_all_formats.sort(key=lambda s: s.lower())

        print_wrap(f"File formats supported by {args.name}:", newline=True)
        max_format_length = max([len(x) for x in l_all_formats])
        print(" "*(max_format_length+4) + "   INPUT  OUTPUT")
        print(" "*(max_format_length+4) + "   -----  ------")
        for file_format in l_all_formats:
            in_yes_or_no = "yes" if file_format in l_input_formats else "no"
            out_yes_or_no = "yes" if file_format in l_output_formats else "no"
            print(f"    {file_format:>{max_format_length}}{in_yes_or_no:>8}{out_yes_or_no:>8}")
        print("")

    if converter_class.allowed_flags is None:
        print_wrap("Information has not been provided about general flags accepted by this converter.", newline=True)
    elif len(converter_class.allowed_flags) == 0:
        print_wrap("This converter does not accept any general flags.", newline=True)
    else:
        print_wrap("Allowed general flags:")
        for flag, help in converter_class.allowed_flags:
            print(f"  {flag}")
            print_wrap(help, width=TERM_WIDTH, initial_indent=" "*4, subsequent_indent=" "*4)
        print("")

    if converter_class.allowed_options is None:
        print_wrap("Information has not been provided about general options accepted by this converter.", newline=True)
    elif len(converter_class.allowed_options) == 0:
        print_wrap("This converter does not accept any general options.", newline=True)
    else:
        print_wrap("Allowed general options:")
        for option, help in converter_class.allowed_options:
            print(f"  {option} <val>")
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
                   f"psdi-data-convert -l {args.name} -f <input_format> -t <output_format>")
    elif mention_input_format:
        print_wrap("For details on input flags and options allowed for a specific format, call:\n"
                   f"psdi-data-convert -l {args.name} -f <input_format> [-t <output_format>]")
    elif mention_output_format:
        print_wrap("For details on output flags and options allowed for a specific format, call:\n"
                   f"psdi-data-convert -l {args.name} -t <output_format> [-f <input_format>]")


def detail_possible_converters(from_format: str, to_format: str):
    """Prints details on converters that can perform a conversion from one format to another
    """
    l_possible_converters = [x for x in get_possible_converters(from_format, to_format)
                             if x[0] in L_REGISTERED_CONVERTERS]

    if len(l_possible_converters) == 0:
        print_wrap(f"No converters are available which can perform a conversion from {from_format} to {to_format}")
        return

    print_wrap(f"The following converters can convert from {from_format} to {to_format}:", newline=True)

    max_converter_len = max([len(name) for name, dos in l_possible_converters])
    print(f"{'    CONVERTER':<{max_converter_len+4}}  DEGREE OF SUCCESS")
    print(f"{'    ---------':<{max_converter_len+4}}  -----------------")
    for name, dos in l_possible_converters:
        print_wrap(f"    {name:<{max_converter_len}}  {dos}", subsequent_indent=" "*(max_converter_len+6))
    print("")

    print_wrap("For details on input/output flags and options allowed by a converter for this conversion, call:")
    print(f"psdi-data-convert -l <converter name> -f {from_format} -t {to_format}")


def detail_converters(args: ConvertArgs):
    """Prints details on available converters for the user.
    """
    if args.name in L_REGISTERED_CONVERTERS:
        return detail_converter_use(args)
    elif args.name != "":
        print(f"ERROR: Converter '{args.name}' not recognized.", file=sys.stderr)
    elif args.from_format and args.to_format:
        return detail_possible_converters(args.from_format, args.to_format)
    print("Available converters: \n\n    " + "\n    ".join(L_REGISTERED_CONVERTERS) + "\n")

    print_wrap("For more details on a converter, call:")
    print("psdi-data-convert -l <converter name>\n")

    print_wrap("For a list of converters that can perform a desired conversion, call:")
    print("psdi-data-convert -l -f <input format> -t <output format>\n")

    print_wrap("For a list of options provided by a converter for a desired conversion, call:")
    print("psdi-data-convert -l <converter name> -f <input format> -t <output format>")


def run_from_args(args: ConvertArgs):
    """Workhorse function to perform primary execution of this script, using the provided parsed arguments.

    Parameters
    ----------
    args : ConvertArgs
        The parsed arguments for this script.
    """

    # Check if we've been asked to list options
    if args.list:
        return detail_converters(args)

    data = {'success': 'unknown',
            'from_flags': args.from_flags,
            'to_flags': args.to_flags,
            'from_arg_flags': '',
            'from_args': '',
            'to_arg_flags': '',
            'to_args': '',
            'coordinates': args.coord_gen,
            'coordOption': args.coord_gen_qual,
            'upload_file': 'true'}

    for filename in args.l_args:

        # Search for the file in the input directory
        qualified_filename = os.path.join(args.input_dir, filename)
        if not os.path.isfile(qualified_filename):
            # Check if we can add the format to it as an extension to find it
            ex_extension = f".{args.from_format}"
            if not qualified_filename.endswith(ex_extension):
                qualified_filename += ex_extension
                if not os.path.isfile(qualified_filename):
                    print(f"ERROR: Cannot find file {filename+ex_extension} in directory {args.input_dir}",
                          file=sys.stderr)
                    continue
            else:
                print(f"ERROR: Cannot find file {filename} in directory {args.input_dir}", file=sys.stderr)
                continue

        if not args.quiet:
            print(f"Converting {filename} to {args.to_format}...")

        try:
            run_converter(filename=qualified_filename,
                          to_format=args.to_format,
                          from_format=args.from_format,
                          name=args.name,
                          data=data,
                          use_envvars=False,
                          upload_dir=args.input_dir,
                          download_dir=args.output_dir,
                          log_file=args.log_file,
                          log_mode=args.log_mode,
                          log_level=args.log_level,
                          delete_input=args.delete_input,
                          refresh_local_log=False)
        except FileConverterAbortException as e:
            print(f"ERROR: Attempt to convert file {filename} aborted with status code {e.status_code} and message: " +
                  f"\n{e}\n", file=sys.stderr)
            continue
        except Exception as e:
            print(f"ERROR: Attempt to convert file {filename} failed with exception type {type(e)} and message: " +
                  f"\n{e}\n", file=sys.stderr)
            continue

        if not args.quiet:
            sys.stdout.write("Success!\n")


def main():
    """Standard entry-point function for this script.
    """

    # If no inputs were provided, print a message about usage
    if len(sys.argv) == 1:
        print_wrap("See the README.md file for information on using this utility and examples of basic usage, or for "
                   "detailed explanation of arguments call:")
        print("psdi-data-convert -h")
        exit(1)

    try:
        args = parse_args()
    except FileConverterHelpException as e:
        # If we get a Help exception, it's likely due to user error, so don't bother them with a traceback and simply
        # print the message to stderr
        print(f"ERROR: {e}", file=sys.stderr)
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
