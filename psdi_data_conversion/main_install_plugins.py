#!/usr/bin/env python3

"""psdi_data_conversion/main_install_plugins.py
=============

Created 2026-07-07 by Bryan Gillis.

Entry-point file for the script to install converter plugins.
"""

import os
from argparse import ArgumentParser

from psdi_data_conversion.converters import base as converters_base

PLUGIN_DATAFILE = "data.json"


def get_argument_parser():
    """Get an argument parser for this script.

    Returns
    -------
    parser : ArgumentParser
        An argument parser set up with the allowed command-line arguments for this script.
    """

    parser = ArgumentParser()

    parser.add_argument("-f", "--force", type=str, default=None,
                        help="Assume that all provided formats are new and don't ask for confirmation if they resemble "
                        "any existing formats")

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

    # Get the directory containing converter plugins
    conv_path = os.path.split(os.path.realpath(converters_base.__file__))[0]

    for dir in os.listdir(conv_path):
        qual_dir = os.path.join(conv_path, dir)
        if not os.path.isdir(qual_dir) or not os.path.isfile(os.path.join(qual_dir, "__init__.py")):
            continue
        # TODO
        pass


def main():
    """Standard entry-point function for this script.
    """

    args = parse_args()

    run_from_args(args)


if __name__ == "__main__":

    main()
