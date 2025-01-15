"""@file psdi-data-conversion/tests/cli_test.py

Created 2025-01-15 by Bryan Gillis.

Tests of the command-line interface
"""

import shlex
import sys
from unittest.mock import patch

from psdi_data_conversion.main import parse_args


def get_parsed_args(s):
    """Performs argument parsing on a string which represents what the arguments would be after the function call
    """
    s = "-c " + s
    with patch.object(sys, 'argv', shlex.split(s)):
        return parse_args


def test_input_validity():
    """Unit tests to ensure that the CLI properly checks for valid input
    """
    # TODO
    pass


def test_input_processing():
    """Unit tests to ensure that the CLI properly processes input arguments to determine values that are needed but
    weren't provided
    """
    # TODO
    pass
