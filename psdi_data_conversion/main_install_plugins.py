#!/usr/bin/env python3

"""psdi_data_conversion/main_install_plugins.py
=============

Created 2026-07-07 by Bryan Gillis.

Entry-point file for the script to install converter plugins.
"""

import json
import os
from argparse import ArgumentParser
from collections import OrderedDict
from copy import deepcopy
from itertools import product
from uuid import uuid4

from psdi_data_conversion import database as db
from psdi_data_conversion.converters import base as converters_base
from psdi_data_conversion.utils import (JsonDict, JsonMainDict, TextColors, confirm_editable_mode, get_wrapped_str,
                                        print_wrap)

# Constants
PLUGIN_DATAFILE = "data.json"
FORMATS_DATAFILE = "formats.json"

MSG_DOS_NOT_TESTED = "not tested"

L_CONVERTER_EXCLUDE_DIRS = ["example", "template", "script_template"]

# The maximum ID of an extra format which we assume is a manually-assigned ID rather than a UUID
THRESHOLD_FORMAT_ID = 9999

# Sorting orders for dicts in the JSON output
L_CONVERTER_SORT_ORDER = [db.DB_NAME_KEY, db.DB_DESCRIPTION_KEY, db.DB_FURTHER_INFO_KEY, db.DB_ID_KEY,
                          db.DB_URL_KEY, db.DB_KEY_PREFIX_KEY, db.DB_SUPPORT_AMBIG_EXT_KEY]
L_CONVERTS_TO_SORT_ORDER = [db.DB_CONV_ID_KEY, db.DB_IN_ID_KEY, db.DB_OUT_ID_KEY, db.DB_SUCCESS_KEY]
L_FORMATS_SORT_ORDER = [db.DB_FORMAT_EXT_KEY, db.DB_FORMAT_NOTE_KEY, db.DB_ID_KEY, db.DB_FORMAT_C2X_KEY,
                        db.DB_FORMAT_COMP_KEY, db.DB_FORMAT_CONN_KEY, db.DB_FORMAT_2D_KEY, db.DB_FORMAT_3D_KEY,
                        db.DB_FORMAT_CONFIRMED_NEW_KEY]
L_ARG_INFO_ORDER = [db.DB_FLAG_KEY, db.DB_BRIEF_KEY, db.DB_DESCRIPTION_KEY, db.DB_FURTHER_INFO_KEY, db.DB_ID_KEY]

REPLACEME_ARG_IN_OUT_ID = "REPLACEME_ARG_IN_OR_OUT_ID"
L_ARG_FORMATS_INFO_ORDER = [db.DB_FORMAT_ID_KEY, REPLACEME_ARG_IN_OUT_ID]


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

    parser.add_argument("--sort-only", action="store_true", help="If set, will not install plugins, and will only "
                        "apply sorting to the currently-installed database file.")

    args = parser.parse_args()

    return args


def get_format_info_str(format_info: JsonDict):
    """Returns a string which details information about a format"""
    return (f"{format_info[db.DB_FORMAT_EXT_KEY]} (ID: {format_info[db.DB_ID_KEY]}): "
            f"{format_info[db.DB_FORMAT_NOTE_KEY]}")


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


def sort_db(db_path: str):
    """Sort the existing database"""
    db_out: JsonMainDict = json.load(open(db_path))
    db_out = get_sorted_dict(db_out)

    sort_json_list(db_out[db.DB_FORMATS_KEY], L_FORMATS_SORT_ORDER)
    for i, d in enumerate(db_out[db.DB_FORMATS_KEY]):
        db_out[db.DB_FORMATS_KEY][i] = get_sorted_dict(d, L_FORMATS_SORT_ORDER)

    sort_json_list(db_out[db.DB_CONVERTERS_KEY], L_CONVERTER_SORT_ORDER)
    for i, d in enumerate(db_out[db.DB_CONVERTERS_KEY]):
        db_out[db.DB_CONVERTERS_KEY][i] = get_sorted_dict(d, L_CONVERTER_SORT_ORDER)

    sort_json_list(db_out[db.DB_CONVERTS_TO_KEY], L_CONVERTS_TO_SORT_ORDER)
    for i, d in enumerate(db_out[db.DB_CONVERTS_TO_KEY]):
        db_out[db.DB_CONVERTS_TO_KEY][i] = get_sorted_dict(d, L_CONVERTS_TO_SORT_ORDER)

    # Get the directory containing converter plugins
    conv_path = os.path.split(os.path.realpath(converters_base.__file__))[0]
    for dir in os.listdir(conv_path):
        if dir in ("example", "template", "script_template"):
            continue

        qual_dir = os.path.join(conv_path, dir)
        if not os.path.isdir(qual_dir) or not os.path.isfile(os.path.join(qual_dir, PLUGIN_DATAFILE)):
            continue

        # Load data about the converter and add it to the database
        db_conv: JsonMainDict = json.load(open(os.path.join(qual_dir, PLUGIN_DATAFILE)))
        d_conv_info: JsonDict = db_conv[db.DB_CONVERTER_KEY]
        conv_prefix: str = d_conv_info[db.DB_KEY_PREFIX_KEY]

        if not conv_prefix:
            continue
        for args_str, in_or_out in product(("flags", "options"), ("in", "out")):
            out_args_str = args_str
            if args_str == "options":
                out_args_str = "argflags"
            out_key = f"{conv_prefix}{out_args_str}_{in_or_out}"
            out_format_key = f"{conv_prefix}format_to_{out_args_str}_{in_or_out}"

            arg_in_out_id = f"{conv_prefix}{out_args_str}_{in_or_out}_id"

            l_arg_format_order = deepcopy(L_ARG_FORMATS_INFO_ORDER)
            l_arg_format_order = [x if x != REPLACEME_ARG_IN_OUT_ID else arg_in_out_id
                                  for x in l_arg_format_order]

            sort_json_list(db_out[out_key], L_ARG_INFO_ORDER)
            for i, d in enumerate(db_out[out_key]):
                db_out[out_key][i] = get_sorted_dict(d, L_ARG_INFO_ORDER)

            sort_json_list(db_out[out_format_key], l_arg_format_order)
            for i, d in enumerate(db_out[out_format_key]):
                db_out[out_format_key][i] = get_sorted_dict(d, l_arg_format_order)

    json.dump(db_out, open(db_path, "w"), indent=4)


def get_support_ambig_ext(db_conv: JsonMainDict, db_out: JsonMainDict):
    """Determine whether or not the converter supports any ambiguous extensions"""
    s_in_formats: set[int] = set(db_conv[db.DB_SUPPORTED_FORMATS_KEY] + db_conv[db.DB_IN_ONLY_FORMATS_KEY])
    s_out_formats: set[int] = set(db_conv[db.DB_SUPPORTED_FORMATS_KEY] + db_conv[db.DB_OUT_ONLY_FORMATS_KEY])

    l_formats: list[JsonDict] = db_out[db.DB_FORMATS_KEY]

    s_in_format_exts: set[str] = set()
    s_out_format_exts: set[str] = set()

    for d_format_info in l_formats:
        id: int = d_format_info[db.DB_ID_KEY]
        ext: str = d_format_info[db.DB_FORMAT_EXT_KEY]
        if id in s_in_formats:
            if ext in s_in_format_exts:
                return True
            s_in_format_exts.add(ext)
        if id in s_out_formats:
            if ext in s_out_format_exts:
                return True
            s_out_format_exts.add(ext)

    return False


def run_from_args(args):
    """Workhorse function to perform primary execution of this script, using the provided parsed arguments.

    Parameters
    ----------
    args : Namespace
        The parsed arguments for this script.
    """

    # Check if the package is installed in editable mode
    confirm_editable_mode()

    db_out: JsonMainDict = {}

    # Get the main database path and directory
    db_path = db.get_database_path()
    db_dir = os.path.split(db_path)[0]

    if args.sort_only:
        return sort_db(db_path)

    # Make a dict of converter paths and DBs we want to process
    conv_parent_path = os.path.split(os.path.realpath(converters_base.__file__))[0]
    d_conv_db: dict[str, JsonMainDict] = {}
    s_changed_conv_dbs: set[str] = set()
    for dir in os.listdir(conv_parent_path):
        if dir in L_CONVERTER_EXCLUDE_DIRS:
            continue

        qual_conv_path = os.path.join(conv_parent_path, dir)
        if not os.path.isdir(qual_conv_path) or not os.path.isfile(os.path.join(qual_conv_path, PLUGIN_DATAFILE)):
            continue

        d_conv_db[qual_conv_path] = json.load(open(os.path.join(qual_conv_path, PLUGIN_DATAFILE)))

    # Load the formats data and check and process new formats
    l_format_info: list[JsonDict] = json.load(open(os.path.join(db_dir, FORMATS_DATAFILE)))[db.DB_FORMATS_KEY]
    d_format_info_for_ext: dict[str, list[JsonDict]] = {}
    for format_info in l_format_info:
        ext = format_info[db.DB_FORMAT_EXT_KEY]
        if ext not in d_format_info_for_ext:
            d_format_info_for_ext[ext] = [format_info]
        else:
            d_format_info_for_ext[ext].append(format_info)
    format_info_updated = False

    questionable_formats_found = False
    first_questionable_format_found = True
    for qual_conv_path, db_conv in d_conv_db.items():

        d_format_id_changes: dict[int, int] = {}
        d_questionable_formats: dict[int, tuple[JsonDict, list[JsonDict]]] = {}
        l_extra_format_info: list[JsonDict] = db_conv[db.DB_EXTRA_FORMATS_KEY]

        for extra_format_info in l_extra_format_info:

            # Check if this format is questionable as to if it already exists in the database
            format_id: int = extra_format_info[db.DB_ID_KEY]
            ext: str = extra_format_info[db.DB_FORMAT_EXT_KEY]
            if (not args.force and format_id <= THRESHOLD_FORMAT_ID and ext in d_format_info_for_ext
                    and not extra_format_info.get(db.DB_FORMAT_CONFIRMED_NEW_KEY)):
                questionable_formats_found = True
                d_questionable_formats[format_id] = (extra_format_info, d_format_info_for_ext[ext])
                continue

            # This format is new, so generate a UUID for it if necessary
            if format_id <= THRESHOLD_FORMAT_ID:
                d_format_id_changes[format_id] = uuid4().int

        # If we found any questionable formats for this converter, report them now
        if d_questionable_formats:
            if first_questionable_format_found:
                first_questionable_format_found = False
                print_wrap(f"{TextColors.WARNING}!!! ALERT !!!{TextColors.ENDC}\n"
                           f"{TextColors.WARNING}-------------{TextColors.ENDC}\n"
                           "The following formats provided by the converter "
                           f"'{db_conv[db.DB_CONVERTER_KEY][db.DB_NAME_KEY]}' might already exist in the database. For "
                           "each, please check against the provided list of possible matches.\n")
                print(get_wrapped_str("- If it is indeed one of those, remove it from the list of extra formats in the "
                                      "converter database file", initial_indent="", subsequent_indent=" "*2) +
                      f" {TextColors.OKCYAN}{os.path.join(qual_conv_path, PLUGIN_DATAFILE)}{TextColors.ENDC}\n" +
                      get_wrapped_str("and update references to its ID in that file to instead use the ID of its "
                                      "entry in the database\n", initial_indent=" "*2, subsequent_indent=" "*2)
                      )
                print_wrap(f"- If it is not one of those, add a line '{TextColors.WARNING}\"confirmed_new\": "
                           f"true{TextColors.ENDC}' to its entry in the converter database file\n",
                           initial_indent="", subsequent_indent=" "*2)
                print_wrap("Once this is done for all formats listed here, rerun this script. Alternatively, if you "
                           "confirm that all listed formats are new, you can rerun the script with the "
                           f"'{TextColors.WARNING}-f/--force{TextColors.ENDC}' flag.\n\n---\n")
            else:
                print_wrap(f"\n\n------\n\nThe following formats provided by the converter "
                           f"'{db_conv[db.DB_CONVERTER_KEY][db.DB_NAME_KEY]}' might already exist in the database:"
                           "\n\n---\n")

            l_format_matches_strings: list[str] = []
            for extra_format_info, l_matching_format_info in d_questionable_formats.values():
                format_matches_string = ("Extra format:\n- " +
                                         get_wrapped_str(get_format_info_str(extra_format_info),
                                                         initial_indent="", subsequent_indent=" "*2) +
                                         "\n\nPotential matches:\n")
                format_matches_string += "\n".join(["- " + get_wrapped_str(get_format_info_str(format_info),
                                                                           initial_indent="", subsequent_indent=" "*2)
                                                   for format_info in l_matching_format_info])
                l_format_matches_strings.append(format_matches_string)
            print_wrap("\n\n---\n\n".join(l_format_matches_strings))

            continue

        # If we get here, all formats look good. Now update the IDs of any formats to UUIDs wherever they appear, if
        # necessary

        if len(l_extra_format_info) > 0:
            format_info_updated = True

        for extra_format_info in l_extra_format_info:
            current_id: int = extra_format_info[db.DB_ID_KEY]
            if current_id in d_format_id_changes:
                extra_format_info[db.DB_ID_KEY] = d_format_id_changes[current_id]

        for l_format_ids in (db_conv[db.DB_SUPPORTED_FORMATS_KEY],
                             db_conv[db.DB_IN_ONLY_FORMATS_KEY],
                             db_conv[db.DB_OUT_ONLY_FORMATS_KEY]):
            for i, current_id in enumerate(l_format_ids):
                if current_id in d_format_id_changes:
                    l_format_ids[i] = d_format_id_changes[current_id]

        for l_conversions in (db_conv[db.DB_SUPPORTED_CONVERSIONS_KEY],
                              db_conv[db.DB_UNSUPPORTED_CONVERSIONS_KEY],):
            for d_conversion_info in l_conversions:
                current_in_id: int = d_conversion_info[db.DB_IN_ID_KEY]
                if current_in_id in d_format_id_changes:
                    d_conversion_info[db.DB_IN_ID_KEY] = d_format_id_changes[current_in_id]
                current_out_id: int = d_conversion_info[db.DB_OUT_ID_KEY]
                if current_out_id in d_format_id_changes:
                    d_conversion_info[db.DB_OUT_ID_KEY] = d_format_id_changes[current_out_id]

        for l_args in (db_conv[db.DB_IN_FLAGS_KEY],
                       db_conv[db.DB_OUT_FLAGS_KEY],
                       db_conv[db.DB_IN_OPTIONS_KEY],
                       db_conv[db.DB_OUT_OPTIONS_KEY]):
            for d_arg_info in l_args:
                l_format_ids: list[int] = d_arg_info[db.DB_FORMAT_ID_LIST_KEY]
                for i, current_id in enumerate(l_format_ids):
                    if current_id in d_format_id_changes:
                        l_format_ids[i] = d_format_id_changes[current_id]

        # Mark if this converter DB has changed from any format IDs changing
        if d_format_id_changes:
            s_changed_conv_dbs.add(qual_conv_path)

        # Add all extra formats to the database now
        for format_info in l_extra_format_info:
            assert format_info[db.DB_ID_KEY] > THRESHOLD_FORMAT_ID
            if db.DB_FORMAT_CONFIRMED_NEW_KEY in format_info:
                del format_info[db.DB_FORMAT_CONFIRMED_NEW_KEY]
            l_format_info.append(format_info)
        db_conv[db.DB_EXTRA_FORMATS_KEY] = []

    # If we found any questionable formats, end execution here to let the user deal with them appropriately
    if questionable_formats_found:
        return

    # Sort the format info and add it to the output database
    l_format_info = [get_sorted_dict(x, L_FORMATS_SORT_ORDER) for x in l_format_info]
    sort_json_list(l_format_info, L_FORMATS_SORT_ORDER)
    db_out[db.DB_FORMATS_KEY] = l_format_info

    # Update the main format info file if any changes have been made to it
    if format_info_updated:
        json.dump({db.DB_FORMATS_KEY: l_format_info},
                  open(os.path.join(db_dir, FORMATS_DATAFILE), "w"), indent=4)

    # Create initial entries for the output dict
    l_db_converters: list[JsonDict] = []
    db_out[db.DB_CONVERTERS_KEY] = l_db_converters
    l_db_converts_to: list[JsonDict] = []
    db_out[db.DB_CONVERTS_TO_KEY] = l_db_converts_to

    for qual_conv_path, db_conv in d_conv_db.items():

        # Load data about the converter and add it to the database
        d_conv_info: JsonDict = db_conv[db.DB_CONVERTER_KEY]
        conv_prefix: str | None = d_conv_info[db.DB_KEY_PREFIX_KEY]

        # Generate any missing data
        if d_conv_info[db.DB_ID_KEY] is None:
            s_changed_conv_dbs.add(qual_conv_path)
            d_conv_info[db.DB_ID_KEY] = uuid4().int

        if d_conv_info[db.DB_SUPPORT_AMBIG_EXT_KEY] is None:
            s_changed_conv_dbs.add(qual_conv_path)
            d_conv_info[db.DB_SUPPORT_AMBIG_EXT_KEY] = get_support_ambig_ext(db_conv, db_out)

        if qual_conv_path in s_changed_conv_dbs:
            json.dump(db_conv, open(os.path.join(qual_conv_path, PLUGIN_DATAFILE), "w"), indent=4)

        conv_id: int = d_conv_info[db.DB_ID_KEY]

        # Rename keys as appropriate
        d_conv_info[db.DB_DESCRIPTION_KEY] = d_conv_info.pop(db.DB_DESC_KEY)
        d_conv_info[db.DB_FURTHER_INFO_KEY] = d_conv_info.pop(db.DB_INFO_KEY)

        l_db_converters.append(get_sorted_dict(d_conv_info, L_CONVERTER_SORT_ORDER))

        # Determine possible conversions and add them all to the database
        s_in_formats = {*db_conv[db.DB_SUPPORTED_FORMATS_KEY]}.union({*db_conv[db.DB_IN_ONLY_FORMATS_KEY]})
        s_out_formats = {*db_conv[db.DB_SUPPORTED_FORMATS_KEY]}.union({*db_conv[db.DB_OUT_ONLY_FORMATS_KEY]})

        d_supported_conversions: dict[tuple[int, int], str] = {(x[db.DB_IN_ID_KEY], x[db.DB_OUT_ID_KEY]):
                                                               x[db.DB_SUCCESS_KEY]
                                                               for x in db_conv[db.DB_SUPPORTED_CONVERSIONS_KEY]}
        s_unsupported_conversions: set[tuple[int, int]] = {(x[db.DB_IN_ID_KEY], x[db.DB_OUT_ID_KEY])
                                                           for x in db_conv[db.DB_UNSUPPORTED_CONVERSIONS_KEY]}

        for in_id, out_id in product(s_in_formats, s_out_formats):
            # Skip conversions of any format to itself, any that are labeled as unsupported or supported (the latter
            # will be added in a separate loop to avoid duplicates)
            if (in_id == out_id or (in_id, out_id) in s_unsupported_conversions or
                    (in_id, out_id) in d_supported_conversions):
                continue
            d_conv_to = {db.DB_CONV_ID_KEY: conv_id,
                         db.DB_IN_ID_KEY: in_id,
                         db.DB_OUT_ID_KEY: out_id,
                         db.DB_SUCCESS_KEY: MSG_DOS_NOT_TESTED}
            l_db_converts_to.append(get_sorted_dict(d_conv_to, L_CONVERTS_TO_SORT_ORDER))

        # Now add any conversions which are explicitly labelled as supported
        for (in_id, out_id), dos in d_supported_conversions.items():
            d_conv_to = {db.DB_CONV_ID_KEY: conv_id,
                         db.DB_IN_ID_KEY: in_id,
                         db.DB_OUT_ID_KEY: out_id,
                         db.DB_SUCCESS_KEY: dos}
            l_db_converts_to.append(get_sorted_dict(d_conv_to, L_CONVERTS_TO_SORT_ORDER))

        # Add the info on input/output arguments and which formats support them, if applicable
        if not conv_prefix:
            continue
        for args_str, in_or_out in product(("flags", "options"), ("in", "out")):
            out_args_str = args_str
            if args_str == "options":
                out_args_str = "argflags"
            in_key = f"{in_or_out}_{args_str}"
            out_key = f"{conv_prefix}{out_args_str}_{in_or_out}"
            out_format_key = f"{conv_prefix}format_to_{out_args_str}_{in_or_out}"

            l_arg_info: list[JsonDict] = []
            l_arg_formats: list[JsonDict] = []
            arg_in_out_id = f"{conv_prefix}{out_args_str}_{in_or_out}_id"

            l_arg_format_order = deepcopy(L_ARG_FORMATS_INFO_ORDER)
            l_arg_format_order = [x if x != REPLACEME_ARG_IN_OUT_ID else arg_in_out_id
                                  for x in l_arg_format_order]

            for d_in_arg_info in db_conv[in_key]:
                d_out_arg_info: JsonDict = {
                    db.DB_DESCRIPTION_KEY: d_in_arg_info[db.DB_DESCRIPTION_KEY],
                    db.DB_FLAG_KEY: d_in_arg_info[db.DB_FLAG_KEY],
                    db.DB_FURTHER_INFO_KEY: d_in_arg_info[db.DB_FURTHER_INFO_KEY],
                    db.DB_ID_KEY: d_in_arg_info[db.DB_ID_KEY],
                }
                if db.DB_BRIEF_KEY in d_in_arg_info:
                    d_out_arg_info[db.DB_BRIEF_KEY] = d_in_arg_info[db.DB_BRIEF_KEY]
                d_out_arg_info = get_sorted_dict(d_out_arg_info, L_ARG_INFO_ORDER)
                l_arg_info.append(d_out_arg_info)

                for format_id in d_in_arg_info[db.DB_FORMAT_ID_LIST_KEY]:
                    d_out_arg_format: JsonDict = {
                        db.DB_FORMAT_ID_KEY: format_id,
                        arg_in_out_id: d_in_arg_info[db.DB_ID_KEY]
                    }
                    d_out_arg_format = get_sorted_dict(d_out_arg_format, l_arg_format_order)
                    l_arg_formats.append(d_out_arg_format)

            sort_json_list(l_arg_info, L_ARG_INFO_ORDER)
            db_out[out_key] = l_arg_info

            sort_json_list(l_arg_formats, l_arg_format_order)
            db_out[out_format_key] = l_arg_formats

    # Sort the top-level dict and lists, then save the database
    db_out = get_sorted_dict(db_out)
    sort_json_list(l_db_converters, L_CONVERTER_SORT_ORDER)
    sort_json_list(l_db_converts_to, L_CONVERTS_TO_SORT_ORDER)
    json.dump(db_out, open(db_path, "w"), indent=4)


def main():
    """Standard entry-point function for this script.
    """

    args = parse_args()

    run_from_args(args)


if __name__ == "__main__":

    main()
