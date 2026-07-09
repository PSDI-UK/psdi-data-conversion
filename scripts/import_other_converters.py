#!/usr/bin/env python3

"""scripts/import_other_converters.py
=============

Created 2026-07-09 by Bryan Gillis.

Script to import database information on unsupported converters
"""

import json
import os
from argparse import ArgumentParser
from collections import OrderedDict
from itertools import product

from psdi_data_conversion import database as db
from psdi_data_conversion.converters import base as converters_base

s_exclude_convs = {"atomsk", "c2x", "example", "openbabel", "template", "script_template"}

# Sorting orders for dicts in the JSON output
L_CONVERTER_META_SORT_ORDER = [db.DB_NAME_KEY, db.DB_DESC_KEY, db.DB_INFO_KEY, db.DB_ID_KEY,
                               db.DB_URL_KEY, db.DB_KEY_PREFIX_KEY, db.DB_SUPPORT_AMBIG_EXT_KEY]
L_CONVERTER_SORT_ORDER = ["converter", "extra_formats", "supported_formats", "in_only_formats", "out_only_formats",
                          "supported_conversions", "unsupported_conversions", "in_flags", "in_options", "out_flags",
                          "out_options"]
L_CONVERSION_SORT_ORDER = ["in_id", "out_id", "degree_of_success"]

# Common types
JsonDict = dict[str, None | int | str | bool | dict | list]
JsonMainDict = dict[str, None | int | str | bool | JsonDict | list[JsonDict]]


def get_sorted_dict(d: dict, l_order: list | None = None):
    """Returns an ordered dict with a provided sorting order"""
    if l_order is None:
        return OrderedDict(sorted(d.items(), key=lambda item: item[0]))
    return OrderedDict(sorted(d.items(), key=lambda item: l_order.index(item[0])))


def sort_json_list(l_d: list[JsonDict], l_order: list | None = None):
    """Sorts a list of JSON dicts based on values of keys, using the provided list of descending-order importance of
    keys in sorting
    """
    def get_key(d: JsonDict):
        l_key = [d[x] for x in l_order if x in d]
        for i, key in enumerate(l_key):
            if isinstance(key, str):
                l_key[i] = (key.lower(), key)
        return tuple(l_key)
    l_d.sort(key=get_key)


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

    db_in: JsonMainDict = json.load(open(db.get_database_path()))

    # Get the directory containing converter plugins
    conv_root_path = os.path.split(os.path.realpath(converters_base.__file__))[0]

    for d_conv_in in db_in["converters"]:
        label = d_conv_in["name"].lower().replace(" ", "")
        if label in s_exclude_convs:
            continue

        conv_path = os.path.join(conv_root_path, label)
        os.makedirs(conv_path, exist_ok=True)

        conv_id: int = d_conv_in["id"]

        d_conv_out: JsonDict = {
            "converter": {
                "name": d_conv_in["name"],
                "desc": d_conv_in["description"],
                "id": d_conv_in["id"],
                "url": d_conv_in["url"],
                "info": "",
                "supports_ambiguous_extensions": None,
                "database_key_prefix": None
            },
            "extra_formats": [],
            "in_flags": [],
            "in_options": [],
            "out_flags": [],
            "out_options": []
        }

        s_in_formats: set[int] = set()
        s_out_formats: set[int] = set()
        s_supported_conversions: set[tuple[int, int]] = set()

        l_conversions: list[JsonDict] = []

        for d_conversion_in in db_in["converts_to"]:
            if d_conversion_in["converters_id"] != conv_id:
                continue
            in_id: int = d_conversion_in["in_id"]
            out_id: int = d_conversion_in["out_id"]
            dos: str = d_conversion_in["degree_of_success"]

            s_in_formats.add(in_id)
            s_out_formats.add(out_id)
            s_supported_conversions.add((in_id, out_id))

            if dos != "not tested":
                d_conversion_out: JsonDict = {
                    "in_id": in_id,
                    "out_id": out_id,
                    "degree_of_success": dos
                }
                d_conversion_out = get_sorted_dict(d_conversion_out, L_CONVERSION_SORT_ORDER)
                l_conversions.append(d_conversion_out)

        sort_json_list(l_conversions, L_CONVERSION_SORT_ORDER)

        s_supported_formats = s_in_formats.intersection(s_out_formats)
        s_in_only_formats = s_in_formats.difference(s_out_formats)
        s_out_only_formats = s_out_formats.difference(s_in_formats)

        l_supported_formats = list(s_supported_formats)
        l_supported_formats.sort()
        d_conv_out["supported_formats"] = l_supported_formats

        l_in_only_formats = list(s_in_only_formats)
        l_in_only_formats.sort()
        d_conv_out["in_only_formats"] = l_in_only_formats

        l_out_only_formats = list(s_out_only_formats)
        l_out_only_formats.sort()
        d_conv_out["out_only_formats"] = l_out_only_formats

        d_conv_out["supported_conversions"] = l_conversions

        l_unsupported_conversions: list[JsonDict] = []
        for in_id, out_id in product(s_in_formats, s_out_formats):
            if in_id == out_id or (in_id, out_id) in s_supported_conversions:
                continue
            d_conversion = {
                "in_id": in_id,
                "out_id": out_id
            }
            d_conversion = get_sorted_dict(d_conversion, L_CONVERSION_SORT_ORDER)
            l_unsupported_conversions.append(d_conversion)

        sort_json_list(l_unsupported_conversions, L_CONVERSION_SORT_ORDER)
        d_conv_out["unsupported_conversions"] = l_unsupported_conversions

        d_conv_out["converter"] = get_sorted_dict(d_conv_out["converter"], L_CONVERTER_META_SORT_ORDER)
        d_conv_out = get_sorted_dict(d_conv_out, L_CONVERTER_SORT_ORDER)

        db_out_path = os.path.join(conv_path, "data.json")
        json.dump(d_conv_out, open(db_out_path, "w"), indent=4)


def main():
    """Standard entry-point function for this script.
    """

    args = parse_args()

    run_from_args(args)


if __name__ == "__main__":

    main()
