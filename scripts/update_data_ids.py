#!/usr/bin/env python3

"""scripts/update_data_ids.py
=============

Created 2026-06-19 by Bryan Gillis.

Script to convert converter IDs in the database to new UUIDs
"""

import json
import logging
from argparse import ArgumentParser
from copy import deepcopy
from itertools import product
from typing import Any
from uuid import uuid4

from psdi_data_conversion import database as db

L_DB_CONVERTER_PREFIXES = ["ob"]
L_DB_ARG_SUFFIXES = [db.DB_IN_FLAGS_FORMATS_KEY_BASE, db.DB_OUT_FLAGS_FORMATS_KEY_BASE,
                     db.DB_IN_OPTIONS_FORMATS_KEY_BASE, db.DB_OUT_OPTIONS_FORMATS_KEY_BASE]

logger = logging.getLogger(__name__)


def get_argument_parser():
    """Get an argument parser for this script.

    Returns
    -------
    parser : ArgumentParser
        An argument parser set up with the allowed command-line arguments for this script.
    """

    parser = ArgumentParser()

    parser.add_argument("out", type=str, help="The output modified database file")
    parser.add_argument("id_changes", type=str, help="Saved dict of ID changes")

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

    # Load the database
    qual_db_path = db.get_database_path()
    d_db_in: dict[str, Any] = json.load(open(qual_db_path))
    d_db_out = deepcopy(d_db_in)

    # Start by changing all format IDs to UUIDs, and keep a dict of the changes
    d_id_changes: dict[int, int] = {}

    for d_converter_info_in, d_converter_info_in in zip(d_db_in[db.DB_CONVERTERS_KEY], d_db_out[db.DB_CONVERTERS_KEY]):
        id_in: int = d_converter_info_in[db.DB_ID_KEY]
        id_out = uuid4().int
        d_converter_info_in[db.DB_ID_KEY] = id_out
        d_id_changes[id_in] = id_out

    # Now update the IDs wherever else they appear into the new IDs
    # for prefix, suffix in product(L_DB_CONVERTER_PREFIXES, L_DB_ARG_SUFFIXES):
    #     key = prefix + suffix
    #     for d_arg_info_in, d_arg_info_out in zip(d_db_in[key], d_db_out[key]):
    #         d_arg_info_out[db.DB_FORMAT_ID_KEY] = d_id_changes[d_arg_info_in[db.DB_FORMAT_ID_KEY]]
    for d_conversion_info_in, d_conversion_info_out in zip(d_db_in[db.DB_CONVERTS_TO_KEY],
                                                           d_db_out[db.DB_CONVERTS_TO_KEY]):
        d_conversion_info_out[db.DB_CONV_ID_KEY] = d_id_changes[d_conversion_info_in[db.DB_CONV_ID_KEY]]

    # Write the updated database file and ID changes dict
    json.dump(d_db_out, open(args.out, "w"))
    json.dump(d_id_changes, open(args.id_changes, "w"))

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
