#!/usr/bin/env python3

"""scripts/import_args.py
=============

Created 2026-07-02 by Bryan Gillis.

Import supported arguments and formats they apply to into new database format
"""

import json
import os
import sys
from argparse import ArgumentParser
from copy import deepcopy
from itertools import product

from psdi_data_conversion.converters import base as converters_base
from psdi_data_conversion.database import (get_database_path, get_in_format_args, get_out_format_args,
                                           get_possible_formats)


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
    d_db_data = json.load(open(get_database_path()))

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
        l_formats = list({*l_in_formats}.union({*l_out_formats}))
        l_arg_types = (("flag", 0), ("option", 1))
        l_in_out = (("in", get_in_format_args), ("out", get_out_format_args))

        d_arg_applies_to: dict[tuple[str, int, int, str], list[int]] = {}

        for (format_info,
             (arg_type_str, arg_type_index),
             (in_out_str, in_out_func)) in product(l_formats,
                                                   l_arg_types,
                                                   l_in_out):
            l_args = in_out_func(conv_name, format_info)[arg_type_index]
            for arg in l_args:
                if not d_arg_applies_to.get((arg_type_str, arg_type_index, arg.id, in_out_str)):
                    d_arg_applies_to[(arg_type_str, arg_type_index, arg.id, in_out_str)] = []
                d_arg_applies_to[(arg_type_str, arg_type_index, arg.id, in_out_str)].append(format_info.id)

        for (arg_type_str, arg_type_index, arg_id, in_out_str), l_applies_to in d_arg_applies_to.items():
            l_applies_to.sort()

            db_arg_type_str: str = "flag"
            if arg_type_str == "option":
                db_arg_type_str = "argflag"
            if not d_out_data.get(f"{in_out_str}_{arg_type_str}"):
                d_out_data[f"{in_out_str}_{arg_type_str}"] = d_db_data[
                    f"{db_key_prefix}{db_arg_type_str}s_{in_out_str}"]

            for d_out_arg_info in d_out_data[f"{in_out_str}_{arg_type_str}"]:
                if not d_out_arg_info["id"] == arg_id:
                    continue
                d_out_arg_info["format_ids"] = l_applies_to
                break
            else:
                print(f"Arg ID {arg_id} not found. {arg_type_str=}, {arg_type_index=}, {in_out_str=}", file=sys.stderr)
                exit(1)

        json.dump(d_out_data, open(datafile, "w"))


def main():
    """Standard entry-point function for this script.
    """

    args = parse_args()

    run_from_args(args)


if __name__ == "__main__":

    main()
