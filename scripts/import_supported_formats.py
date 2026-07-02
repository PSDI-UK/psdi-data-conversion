#!/usr/bin/env python3

"""scripts/import_supported_formats.py
=============

Created 2026-07-02 by Bryan Gillis.

Import supported formats into new database format
"""

import json
import os
from argparse import ArgumentParser
from copy import deepcopy

from psdi_data_conversion.converters import base as converters_base
from psdi_data_conversion.database import get_possible_formats


def get_argument_parser():
    """Get an argument parser for this script.

    Returns
    -------
    parser : ArgumentParser
        An argument parser set up with the allowed command-line arguments for this script.
    """

    parser = ArgumentParser()

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

    # Loop through plugins
    for dir in os.listdir(conv_path):
        qual_dir = os.path.join(conv_path, dir)
        if (not os.path.isdir(qual_dir) or not os.path.isfile(os.path.join(qual_dir, "__init__.py"))
                or dir in ("example", "template", "script_template")):
            continue
        datafile = os.path.join(qual_dir, "data.json")
        d_in_data = json.load(open(datafile))
        d_out_data = deepcopy(d_in_data)

        # Get the formats supported by the converter as input-output, input-only, and output-only
        in_formats, out_formats = get_possible_formats(d_in_data["converter"]["name"])
        s_in_formats, s_out_formats = {*in_formats}, {*out_formats}
        l_supported_formats = [*s_in_formats.union(s_out_formats)]
        l_in_only_formats = [*s_in_formats.difference(s_out_formats)]
        l_out_only_formats = [*s_out_formats.difference(s_in_formats)]

        l_supported_formats.sort(), l_in_only_formats.sort(), l_out_only_formats.sort()

        d_out_data["supported_formats"] = [x.id for x in l_supported_formats]
        d_out_data["in_only_formats"] = [x.id for x in l_in_only_formats]
        d_out_data["out_only_formats"] = [x.id for x in l_out_only_formats]

        json.dump(open(datafile, "w"))


def main():
    """Standard entry-point function for this script.
    """

    args = parse_args()

    run_from_args(args)


if __name__ == "__main__":

    main()
