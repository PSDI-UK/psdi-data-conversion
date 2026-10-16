#!/usr/bin/env python3

"""scripts/list_format_support.py
=============

Created 2026-10-06 by Bryan Gillis.

Script to sort formats by which converters support them
"""

import logging
from argparse import ArgumentParser

from psdi_data_conversion.converter import L_REGISTERED_CONVERTERS
from psdi_data_conversion.database import ConverterInfo, FormatInfo, get_converter_info, get_possible_formats
from psdi_data_conversion.utils import print_wrap

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


def run_from_args(args):
    """Workhorse function to perform primary execution of this script, using the provided parsed arguments.

    Parameters
    ----------
    args : Namespace
        The parsed arguments for this script.
    """

    logger.debug("# Entering function `run_from_args`")

    l_conv = [get_converter_info(x) for x in L_REGISTERED_CONVERTERS]

    d_conv_s_in_fmts: dict[ConverterInfo, set[FormatInfo]] = {}
    d_conv_s_out_fmts: dict[ConverterInfo, set[FormatInfo]] = {}

    # Get sets of all formats supported as input/output by each converter
    for conv in l_conv:
        l_in_fmts, l_out_fmts = get_possible_formats(conv)
        d_conv_s_in_fmts[conv] = set(l_in_fmts)
        d_conv_s_out_fmts[conv] = set(l_out_fmts)

    # Keep track of how many converters support each format
    s_all_in_fmts: set[FormatInfo] = set()
    s_all_out_fmts: set[FormatInfo] = set()
    d_in_fmt_support_count: dict[FormatInfo, int] = {x: 0 for x in s_all_in_fmts}
    d_out_fmt_support_count: dict[FormatInfo, int] = {x: 0 for x in s_all_out_fmts}

    for conv in l_conv:
        s_uniq_in_fmts = d_conv_s_in_fmts[conv]
        s_uniq_out_fmts = d_conv_s_out_fmts[conv]

        s_all_in_fmts = s_all_in_fmts.union(d_conv_s_in_fmts[conv])
        s_all_out_fmts = s_all_out_fmts.union(d_conv_s_out_fmts[conv])

        d_in_fmt_support_count.update({x: 0 for x in s_all_in_fmts if x not in d_in_fmt_support_count})
        d_out_fmt_support_count.update({x: 0 for x in s_all_out_fmts if x not in d_out_fmt_support_count})

        for fmt in d_conv_s_in_fmts[conv]:
            d_in_fmt_support_count[fmt] += 1
        for fmt in d_conv_s_out_fmts[conv]:
            d_out_fmt_support_count[fmt] += 1

        l_other_conv = [x for x in l_conv if x != conv]
        for other_conv in l_other_conv:
            s_uniq_in_fmts = s_uniq_in_fmts.difference(d_conv_s_in_fmts[other_conv])
            s_uniq_out_fmts = s_uniq_out_fmts.difference(d_conv_s_out_fmts[other_conv])

        # List formats supported by only this converter
        for (in_or_out, s_fmts) in (("Input", s_uniq_in_fmts), ("Output", s_uniq_out_fmts)):
            print_wrap(f"{in_or_out} formats only supported by {conv.format_word()}:", newline=True)
            l_fmts = list(s_fmts)
            l_fmts.sort(key=lambda x: x.disambiguated_name)
            if not l_fmts:
                print_wrap("- (None)", newline=True)
                continue
            for fmt in l_fmts:
                print_wrap(f"- {fmt.format_oneline()}")
            print("")

    # List formats supported by at least one converter
    s_in_multiple: set[FormatInfo] = set()
    s_out_multiple: set[FormatInfo] = set()
    for (in_or_out, d_support_count) in (("Input", d_in_fmt_support_count), ("Output", d_out_fmt_support_count)):
        print_wrap(f"{in_or_out} formats supported by more than one converter:", newline=True)
        l_fmts = [key for key, val in d_support_count.items() if val > 1]
        if in_or_out == "Input":
            s_in_multiple = s_in_multiple.union(set(l_fmts))
        else:
            s_out_multiple = s_out_multiple.union(set(l_fmts))
        l_fmts.sort(key=lambda x: x.disambiguated_name)
        if not l_fmts:
            print_wrap("- (None)", newline=True)
            continue
        for fmt in l_fmts:
            print_wrap(f"- {fmt.format_oneline()}")
        print("")

    print_wrap("Formats supported as both input and output by more than one converter:", newline=True)
    l_inout_multiple = list(s_in_multiple.intersection(s_out_multiple))
    l_inout_multiple.sort(key=lambda x: x.disambiguated_name)
    if not l_inout_multiple:
        print_wrap("- (None)", newline=True)
    else:
        for fmt in l_inout_multiple:
            print_wrap(f"- {fmt.format_oneline()}")
        print("")

    # List formats supported by all converters
    s_in_all: set[FormatInfo] = set()
    s_out_all: set[FormatInfo] = set()
    for (in_or_out, d_support_count) in (("Input", d_in_fmt_support_count), ("Output", d_out_fmt_support_count)):
        print_wrap(f"{in_or_out} formats supported by all converters:", newline=True)
        l_fmts = [key for key, val in d_support_count.items() if val == len(l_conv)]
        if in_or_out == "Input":
            s_in_all = s_in_all.union(set(l_fmts))
        else:
            s_out_all = s_out_all.union(set(l_fmts))
        l_fmts.sort(key=lambda x: x.disambiguated_name)
        if not l_fmts:
            print_wrap("- (None)", newline=True)
            continue
        for fmt in l_fmts:
            print_wrap(f"- {fmt.format_oneline()}")
        print("")

    print_wrap("Formats supported as both input and output by all converters:", newline=True)
    l_all_multiple = list(s_in_all.intersection(s_out_all))
    l_all_multiple.sort(key=lambda x: x.disambiguated_name)
    if not l_all_multiple:
        print_wrap("- (None)", newline=True)
    else:
        for fmt in l_all_multiple:
            print_wrap(f"- {fmt.format_oneline()}")
        print("")

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
