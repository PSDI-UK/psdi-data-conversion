#!/usr/bin/env python3

"""@file psdi_data_conversion/main.py

Created 2025-01-14 by Bryan Gillis.

Entry-point file for the command-line interface for data conversion.
"""

import logging
from argparse import ArgumentParser
import os
import sys

from psdi_data_conversion import constants as const
from psdi_data_conversion.converter import L_REGISTERED_CONVERTERS, run_converter
from psdi_data_conversion.converters.base import (FileConverterAbortException, FileConverterInputException,
                                                  get_file_storage)

logger = logging.getLogger(__name__)


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
        self._output_dir: str | None = args.at
        self.name: str = getattr(args, "with")
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
        self.log_level: str = args.log_level

        # Quiet mode is equivalent to logging mode == LOGGING_NONE, so normalize them if either is set
        if self.quiet:
            self.log_mode = const.LOG_NONE
        elif self.log_mode == const.LOG_NONE:
            self.quiet = True

        # Check validity of input

        if self.list:
            # If requesting to list converters, any other arguments can be ignored
            return

        if len(self.l_args) == 0:
            raise FileConverterInputException("One or more names of files to convert must be provided")

        if self._input_dir is not None and not os.path.isdir(self._input_dir):
            raise FileConverterInputException(f"The provided input directory '{self._input_dir}' does not exist as a "
                                              "directory")

        if self.to_format is None:
            raise FileConverterInputException("Output format (-t or --to) must be provided")

        # If the output directory doesn't exist, silently create it
        if self._output_dir is not None and not os.path.isdir(self._output_dir):
            if os.path.exists(self._output_dir):
                raise FileConverterInputException(
                    f"Output directory '{self._output_dir}' exists but is not a directory")
            os.makedirs(self._output_dir, exist_ok=True)

        # Check the converter is recognized
        if self.name not in L_REGISTERED_CONVERTERS:
            raise FileConverterInputException(f"Converter '{self.name}' not recognised")

        # No more than two arguments supplied to --coord-gen
        if args.coord_gen is not None and len(args.coord_gen) > 2:
            raise FileConverterInputException("At most two arguments may be provided to --coord-gen, the mode and "
                                              "quality, e.g. '--coord-gen Gen3D best'")

        # Coordinate generation options are valid
        if self.coord_gen not in const.L_ALLOWED_COORD_GENS:
            raise FileConverterInputException(f"Coordinate generation type '{self.coord_gen}' not recognised. Allowed "
                                              f"types are: {const.L_ALLOWED_COORD_GENS}")
        if self.coord_gen_qual not in const.L_ALLOWED_COORD_GEN_QUALS:
            raise FileConverterInputException(f"Coordinate generation quality '{self.coord_gen_qual}' not recognised. "
                                              f"Allowed qualities are: {const.L_ALLOWED_COORD_GEN_QUALS}")

        # Logging mode is valid
        if self.log_mode not in const.L_ALLOWED_LOG_MODES:
            raise FileConverterInputException(f"ERROR: Unrecognised logging option: {self.log_mode}. Allowed "
                                              f"options are: {const.L_ALLOWED_LOG_MODES}")

    @property
    def from_format(self):
        """If the input file format isn't provided, determine it from the first file in the list.
        """
        if self._from_format is None:
            first_filename = self.l_args[0]
            ext = os.path.splitext(first_filename)[1]
            if len(ext) == 0:
                raise FileConverterInputException("Input file format (-f or --from) was not provided, and cannot "
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
                        raise FileConverterInputException(f"ERROR: Input file {first_filename} cannot be found. Also "
                                                          f"checked for {test_filename}.")

                base = os.path.splitext(first_filename)[0]
                self._log_file = base + const.LOG_EXT
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
    parser.add_argument("-a", "--at", type=str, default=None,
                        help="The directory where output files should be created, default same as input directory.")
    parser.add_argument("-w", "--with", type=str, default="Open Babel",
                        help="The converter to be used (default 'Open Babel').")
    parser.add_argument("-d", "--delete-input", action="store_true",
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
    parser.add_argument("--list", action="store_true",
                        help="If provided alone, lists all available converters. If the name of a converter is "
                             "provided, gives information on the converter and any command-line flags it accepts.")

    # Logging/stdout arguments
    parser.add_argument("-l", "--log-file", type=str, default=None,
                        help="The name of the file to log to. If not provided, the log file will be named after the "
                             "first input file (+'.log') and placed in the current directory.")
    parser.add_argument("--log-mode", type=str, default=const.LOG_SIMPLE,
                        help="How logs should be stores. Allowed values are: \n"
                        "- 'full' - Multi-file logging, only recommended when running as a public web app"
                        "- 'simple' - Logs saved to one file"
                        "- 'stdout' - Output logs and errors only to stdout"
                        "- 'none' - Output only errors to stdout")
    parser.add_argument("-q", "--quiet", action="store_true",
                        help="If set, all output aside from errors will be suppressed and no log file will be "
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


def detail_converter_use(converter_name: str):
    """TODO
    """
    print("Converter use detailing is still TBD")


def detail_converters(l_args: list[str]):
    """Prints details on available converters for the user.
    """
    converter_name = " ".join(l_args)
    if converter_name in L_REGISTERED_CONVERTERS:
        return detail_converter_use(converter_name)
    elif converter_name != "":
        print(f"ERROR: Converter {converter_name} not recognized.", file=sys.stderr)
    print("Available converters are: \n" + "\n".join(L_REGISTERED_CONVERTERS) + "\n" +
          "For more details on a converter, call: \n" +
          "psdi-data-convert --list <Converter name>")


def run_from_args(args: ConvertArgs):
    """Workhorse function to perform primary execution of this script, using the provided parsed arguments.

    Parameters
    ----------
    args : ConvertArgs
        The parsed arguments for this script.
    """

    logger.debug("# Entering function `run_from_args`")

    # Check if we've been asked to list options
    if args.list:
        return detail_converters(args.l_args)

    form = {'token': '1041c0a661d118d5f28e7c6830375dd0',
            'from': args.from_format,
            'to': args.to_format,
            'from_full': args.from_format,
            'to_full': args.to_format,
            'success': 'unknown',
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

        file_storage = get_file_storage(qualified_filename)

        try:
            run_converter(name=args.name,
                          files=file_storage,
                          form=form,
                          file_to_convert=const.FILE_TO_UPLOAD_KEY,
                          use_envvars=False,
                          upload_dir=args.input_dir,
                          download_dir=args.output_dir,
                          log_file=args.log_file,
                          log_mode=args.log_mode,
                          delete_input=args.delete_input,
                          max_file_size=0)
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

    logger.debug("# Exiting function `run_from_args`")


def main():
    """Standard entry-point function for this script.
    """

    args = parse_args()

    if args.log_mode == const.LOG_SIMPLE or args.log_mode == const.LOG_FULL:
        logging.basicConfig(filename=args.log_file)

    logger.info("#")
    logger.info("# Beginning execution of script `%s`", __file__)
    logger.info("#")

    run_from_args(args)

    logger.info("#")
    logger.info("# Finished execution of script `%s`", __file__)
    logger.info("#")


if __name__ == "__main__":

    main()
