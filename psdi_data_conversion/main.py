#!/usr/bin/env python3

"""@file psdi_data_conversion/main.py

Created 2025-01-14 by Bryan Gillis.

Entry-point file for the command-line interface for data conversion.
"""

import logging
from argparse import ArgumentParser

logger = logging.getLogger(__name__)


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
    parser.add_argument("-t", "--to", type=str,
                        help="The output (convert to) file extension (e.g., cmi).")
    parser.add_argument("-a", "--at", type=str, default=None,
                        help="The directory where output files should be created, default same as input directory.")
    parser.add_argument("-w", "--with", type=str, default="Open Babel",
                        help="The converter to be used (default 'Open Babel').")
    parser.add_argument("--flags", type=str, default="",
                        help="Any command-line flags to be provided to the converter. For information on the flags "
                             "accepted by a converter, call this script with '-l <converter name>'.")

    # Keyword arguments for alternative functionality
    parser.add_argument("-l", "--list", type=bool, action="store_true",
                        help="If provided alone, lists all available converters. If the name of a converter is "
                             "provided, gives information on the converter and any command-line flags it accepts.")

    # Logging/stdout arguments
    parser.add_argument("-q", "--quiet", type=bool, action="store_true",
                        help="If set, all output aside from errors will be suppressed and no log file will be "
                             "generated.")
    parser.add_argument("-l", "--log-file", type=str, default=None,
                        help="The name of the file to log to. If not provided, the log file will be named after the "
                             "first input file (+'.log') and placed in the current directory.")
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

    args = parser.parse_args()

    return args


def run_from_args(args):
    """Workhorse function to perform primary execution of this script, using the provided parsed arguments.

    Parameters
    ----------
    args : Namespace
        The parsed arguments for this script.
    """

    logger.debug("# Entering function `run_from_args`")

    print("This is currently a dummy executable, with functionality TBD.")

    logger.debug("# Exiting function `run_from_args`")


def main():
    """Standard entry-point function for this script.
    """

    args = parse_args()

    logging.basicConfig(level=args.log_level)

    logger.info("#")
    logger.info("# Beginning execution of script `%s`", __file__)
    logger.info("#")

    run_from_args(args)

    logger.info("#")
    logger.info("# Finished execution of script `%s`", __file__)
    logger.info("#")


if __name__ == "__main__":

    main()
