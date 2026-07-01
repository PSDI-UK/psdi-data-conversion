#!/usr/bin/env python3

"""psdi_data_conversion/main_create_plugin.py
=============

Created 2026-07-01 by Bryan Gillis.

Entry-point file for the script to create a new converter plugin.
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

    logger.debug("#")
    logger.debug("# Beginning execution of script `%s`", __file__)
    logger.debug("#")

    run_from_args(args)

    logger.debug("#")
    logger.debug("# Finished execution of script `%s`", __file__)
    logger.debug("#")


if __name__ == "__main__":

    main()
