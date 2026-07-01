#!/usr/bin/env python3

"""psdi_data_conversion/main_create_plugin.py
=============

Created 2026-07-01 by Bryan Gillis.

Entry-point file for the script to create a new converter plugin.
"""

import logging
import os
import sys
from argparse import ArgumentParser

from psdi_data_conversion.converters import base as converters_base

PLUGIN_TEMPLATE = "template"
PLUGIN_SCRIPT_TEMPLATE = "script_template"
PLUGIN_PYFILE = "converter.py"
PLUGIN_DATAFILE = "data.json"

logger = logging.getLogger(__name__)


def get_argument_parser():
    """Get an argument parser for this script.

    Returns
    -------
    parser : ArgumentParser
        An argument parser set up with the allowed command-line arguments for this script.
    """

    parser = ArgumentParser()

    parser.add_argument("name", type=str, nargs="+",
                        help="The name of the plugin to be created, e.g. 'Open Babel'")

    parser.add_argument("--label", type=str, default=None,
                        help="The label for the package (i.e. Python-compatible package name). By default, will "
                        "convert `name` to snake_case (e.g. 'Open Babel' -> 'open_babel')")

    parser.add_argument("--script", action="store_true",
                        help="If set, will create the plugin using the 'ScriptFileConverter' base class, which uses "
                        "a script to run the conversion. The script will by default be named `{label}.sh`")

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

    # Get the name and label from input arguments
    name: str = args.name
    label: str | None = args.label
    if label:
        # Check that the label appears to be properly in snake_case
        if label != label.lower().replace(" ", "_"):
            print(f"ERROR: Label '{label}' is invalid. The label should be in snake_case: All lower-case with "
                  "underscores in place of spaces", file=sys.stderr)
            exit(1)
    else:
        # Create the label by converting the name to snake_case
        label = name.lower().replace(" ", "_")

    # Determine the PascalCase name of the converter (to be used for classes)
    l_name_words = name.split(" ")
    pascal_name = "".join(map(lambda x: x.capitalize(), l_name_words))

    # Get the directory containing converter plugins
    conv_path = os.path.realpath(converters_base.__file__)

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
