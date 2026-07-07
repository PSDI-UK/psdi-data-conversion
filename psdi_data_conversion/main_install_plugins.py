#!/usr/bin/env python3

"""psdi_data_conversion/main_install_plugins.py
=============

Created 2026-07-07 by Bryan Gillis.

Entry-point file for the script to install converter plugins.
"""

import json
import os
from argparse import ArgumentParser
from itertools import product

from psdi_data_conversion import database as db
from psdi_data_conversion.converters import base as converters_base

PLUGIN_DATAFILE = "data.json"
FORMATS_DATAFILE = "formats.json"
UNSUPPORTED_CONVERTERS_DATAFILE = "unsupported_converters.json"

MSG_DOS_NOT_TESTED = "not tested"

JsonDict = dict[str, None | int | str | bool | dict | list]
JsonMainDict = dict[str, None | int | str | bool | JsonDict | list[JsonDict]]


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

    db_out: JsonMainDict = {}

    # Get the main database path and directory
    db_path = db.get_database_path()
    db_dir = os.path.split(db_path)[0]

    # Load the formats data into the output dict
    db_out[db.DB_FORMATS_KEY] = json.load(open(os.path.join(db_dir, FORMATS_DATAFILE)))[db.DB_FORMATS_KEY]

    # Load data on unsupported converters into the output dict, and fill in missing items
    l_db_converters = json.load(
        open(os.path.join(db_dir, UNSUPPORTED_CONVERTERS_DATAFILE)))[db.DB_CONVERTERS_KEY]
    db_out[db.DB_CONVERTERS_KEY] = l_db_converters
    for d_conv_info in l_db_converters:
        if db.DB_DESC_KEY in d_conv_info:
            d_conv_info[db.DB_DESCRIPTION_KEY] = d_conv_info.pop(db.DB_DESC_KEY)
        if db.DB_INFO_KEY in d_conv_info:
            d_conv_info[db.DB_FURTHER_INFO_KEY] = d_conv_info.pop(db.DB_INFO_KEY)
        for key, default in ((db.DB_KEY_PREFIX_KEY, None),
                             (db.DB_DESCRIPTION_KEY, ""),
                             (db.DB_FURTHER_INFO_KEY, ""),
                             (db.DB_SUPPORT_AMBIG_EXT_KEY, False),
                             (db.DB_URL_KEY, "")):
            if key not in d_conv_info:
                d_conv_info[key] = default

    # Create initial entries for the rest of the output dict
    db_out[db.DB_CONVERTERS_KEY] = l_db_converters
    l_db_converts_to: list[JsonDict] = []
    db_out[db.DB_CONVERTS_TO_KEY] = l_db_converts_to

    # Get the directory containing converter plugins
    conv_path = os.path.split(os.path.realpath(converters_base.__file__))[0]

    for dir in os.listdir(conv_path):
        if dir in ("example", "template", "script_template"):
            continue

        qual_dir = os.path.join(conv_path, dir)
        if not os.path.isdir(qual_dir) or not os.path.isfile(os.path.join(qual_dir, "__init__.py")):
            continue

        # Load data about the converter and add it to the database
        db_conv: JsonMainDict = json.load(open(os.path.join(qual_dir, PLUGIN_DATAFILE)))
        conv_id: int = db_conv[db.DB_CONVERTER_KEY][db.DB_ID_KEY]
        d_conv_info: JsonDict = db_conv[db.DB_CONVERTER_KEY]

        # Rename keys as appropriate
        d_conv_info[db.DB_DESCRIPTION_KEY] = d_conv_info.pop(db.DB_DESC_KEY)
        d_conv_info[db.DB_FURTHER_INFO_KEY] = d_conv_info.pop(db.DB_INFO_KEY)

        l_db_converters.append(d_conv_info)

        # Determine possible conversions and add them all to the database
        s_in_formats = {*db_conv[db.DB_SUPPORTED_FORMATS_KEY]}.union({*db_conv[db.DB_IN_ONLY_FORMATS_KEY]})
        s_out_formats = {*db_conv[db.DB_SUPPORTED_FORMATS_KEY]}.union({*db_conv[db.DB_OUT_ONLY_FORMATS_KEY]})

        d_supported_conversions: dict[tuple[int, int], str] = {(x[db.DB_IN_ID_KEY], x[db.DB_IN_ID_KEY]):
                                                               x[db.DB_SUCCESS_KEY]
                                                               for x in db_conv[db.DB_SUPPORTED_CONVERSIONS_KEY]}
        s_unsupported_conversions: set[tuple[int, int]] = {(x[db.DB_IN_ID_KEY], x[db.DB_IN_ID_KEY])
                                                           for x in db_conv[db.DB_UNSUPPORTED_CONVERSIONS_KEY]}

        for in_id, out_id in product(s_in_formats, s_out_formats):
            # Skip conversions of any format to itself, any that are labeled as unsupported or supported (the latter
            # will be added in a separate loop to avoid duplicates)
            if (in_id == out_id or (in_id, out_id) in s_unsupported_conversions or
                    (in_id, out_id) in d_supported_conversions):
                continue
            l_db_converts_to.append({db.DB_CONV_ID_KEY: conv_id,
                                     db.DB_IN_ID_KEY: in_id,
                                     db.DB_OUT_ID_KEY: out_id,
                                     db.DB_SUCCESS_KEY: MSG_DOS_NOT_TESTED})

        # Now add any conversions which are explicitly labelled as supported
        for (in_id, out_id), dos in d_supported_conversions.items():
            l_db_converts_to.append({db.DB_CONV_ID_KEY: conv_id,
                                     db.DB_IN_ID_KEY: in_id,
                                     db.DB_OUT_ID_KEY: out_id,
                                     db.DB_SUCCESS_KEY: dos})

        # TODO: Add argument info here

    # Save the database
    json.dump(db_out, open(db_path.replace(".json", "2.json"), "w"), sort_keys=True, indent=4)


def main():
    """Standard entry-point function for this script.
    """

    args = parse_args()

    run_from_args(args)


if __name__ == "__main__":

    main()
