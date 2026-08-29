#!/usr/bin/env python3

"""scripts/replace_aliases.py
=============

Created 2026-08-27 by Bryan Gillis.
"""

import json
import logging
from argparse import ArgumentParser
from pathlib import Path
from typing import Iterable

from psdi_data_conversion.utils import JsonMainDict

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


def get_cleaned_list(input_list: list[int], l_aliases: Iterable[tuple[int, int]]):
    s = set(input_list)
    for id, alias_for in l_aliases:
        if id in s:
            s.add(alias_for)
            s.remove(id)
    output_list = list(s)
    output_list.sort()
    return output_list


def run_from_args(args):
    """Workhorse function to perform primary execution of this script, using the provided parsed arguments.

    Parameters
    ----------
    args : Namespace
        The parsed arguments for this script.
    """

    logger.debug("# Entering function `run_from_args`")

    l_formats: list[dict] = json.load(open("psdi_data_conversion/static/data/formats.json"))['formats']
    d_aliases = {x['id']: x['alias_of'] for x in l_formats if x.get('alias_of')}

    plugins_path = Path("psdi_data_conversion/converters").resolve()

    for conv_path in plugins_path.iterdir():
        if not conv_path.is_dir():
            return
        data_path = conv_path / "data.json"
        d_conv: JsonMainDict = json.load(open(data_path))

        for key in "supported_formats", "in_only_formats", "out_only_formats":
            d_conv[key] = get_cleaned_list(d_conv[key], d_aliases.items())
        for key in "in_flags", "out_flags", "in_options", "out_options":
            for d_arg in d_conv[key]:
                d_arg["format_ids"] = get_cleaned_list(d_arg["format_ids"], d_aliases.items())

        l_supported_conversions = d_conv["supported_conversions"]
        l_supp_conv = [(x["in_id"], x["out_id"], x["degree_of_success"]) for x in l_supported_conversions]
        s_supp_conv = set(l_supp_conv)
        for in_id, out_id, dos in l_supp_conv:
            in_tuple = in_id, out_id, dos
            if in_id in d_aliases:
                in_id = d_aliases[in_id]
            if out_id in d_aliases:
                out_id = d_aliases[out_id]
            out_tuple = in_id, out_id, dos
            if out_tuple != in_tuple:
                s_supp_conv.remove(in_tuple)
                s_supp_conv.add(out_tuple)
        l_supp_conv = list(s_supp_conv)
        l_supp_conv.sort()
        d_conv["supported_conversions"] = [{"in_id": x[0], "out_id": x[1], "degree_of_success": x[2]}
                                           for x in l_supp_conv]

        l_unsupported_conversions = d_conv["unsupported_conversions"]
        l_unsupp_conv = [(x["in_id"], x["out_id"]) for x in l_unsupported_conversions]
        s_unsupp_conv = set(l_unsupp_conv)
        for in_id, out_id in l_unsupp_conv:
            in_tuple = in_id, out_id
            if in_id in d_aliases:
                in_id = d_aliases[in_id]
            if out_id in d_aliases:
                out_id = d_aliases[out_id]
            out_tuple = in_id, out_id
            if out_tuple != in_tuple:
                s_unsupp_conv.remove(in_tuple)
                s_unsupp_conv.add(out_tuple)
        l_unsupp_conv = list(s_unsupp_conv)
        l_unsupp_conv.sort()
        d_conv["unsupported_conversions"] = [{"in_id": x[0], "out_id": x[1]}
                                             for x in l_unsupp_conv]

        json.dump(d_conv, open(data_path, "w"), indent=4)

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
