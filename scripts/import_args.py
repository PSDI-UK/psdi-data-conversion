#!/usr/bin/env python3

"""scripts/import_args.py
=============

Created 2026-07-02 by Bryan Gillis.

Import supported arguments and formats they apply to into new database format
"""

import json
import os
from argparse import ArgumentParser
from copy import deepcopy
from inspect import ArgInfo
from itertools import product

from psdi_data_conversion.converters import base as converters_base
from psdi_data_conversion.database import FlagInfo, get_in_format_args, get_out_format_args, get_possible_formats


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

        # Get the database key prefix, if one exists - if not, skip this plugic
        db_key_prefix: str | None = d_in_data["converter"]["database_key_prefix"]
        if db_key_prefix is None:
            continue
        conv_name = d_in_data["converter"]["name"]

        l_in_formats, l_out_formats = get_possible_formats(conv_name)
        l_arg_types = (("flag", 0), ("option", 1))
        l_in_out = (("in", get_in_format_args, 0), ("out", get_out_format_args, 1))

        d_arg_applies_to: dict[tuple[ArgInfo, str], list[int]] = {}

        for (l_format_info,
             (arg_type_str, arg_type_index),
             (in_out_str, in_out_func, in_out_index)) in product(zip(l_in_formats, l_out_formats),
                                                                 l_arg_types,
                                                                 l_in_out):
            format_info = l_format_info[in_out_index]
            l_args = in_out_func(conv_name, format_info)[arg_type_index]
            for arg in l_args:
                if not d_arg_applies_to.get((arg, in_out_str)):
                    d_arg_applies_to[(arg, in_out_str)] = []
                d_arg_applies_to[(arg, in_out_str)].append(format_info.id)

        for (arg, in_out_str), l_applies_to in d_arg_applies_to.items():
            arg_str: str
            if isinstance(arg, FlagInfo):
                arg_str = "flag"
            else:
                arg_str = "option"

            l_applies_to.sort()

            d_out_data[f"{in_out_str}_{arg_str}"] = d_in_data[f"{db_key_prefix}{arg_str}s_{in_out_str}"]

            d_out_data[f"{in_out_str}_{arg_str}"]["format_ids"] = l_applies_to

        json.dump(d_out_data, open(datafile, "w"))


def main():
    """Standard entry-point function for this script.
    """

    args = parse_args()

    run_from_args(args)


if __name__ == "__main__":

    main()
