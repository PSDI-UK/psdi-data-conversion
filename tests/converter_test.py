"""
# converter_test.py

Unit tests of the converter class. This module uses the common test specifications defined in
psdi_data_conversion/testing/conversion_test_specs.py so that a common set of conversion tests is performed through
the Python library (this module), the command-line application, and the GUI.
"""

import logging
import math
import os
import pytest

from psdi_data_conversion import constants as const
from psdi_data_conversion.converter import L_REGISTERED_CONVERTERS
from psdi_data_conversion.converters.openbabel import OBFileConverter
from psdi_data_conversion.testing.utils import run_test_conversion_with_library
from psdi_data_conversion.testing import conversion_test_specs as specs


@pytest.fixture(autouse=True)
def setup_test() -> None:
    """Reset global aspects before a test, so that different tests won't interfere with each other"""

    # Remove the global log file if one exists
    try:
        os.remove(const.GLOBAL_LOG_FILENAME)
    except FileNotFoundError:
        pass

    # Clear any existing loggers so new ones will be created fresh
    logging.Logger.manager.loggerDict.clear()


def test_default():
    """Test that the default converter is registered.
    """
    assert const.CONVERTER_DEFAULT in L_REGISTERED_CONVERTERS


def test_basic_conversions():
    """Run a basic set of conversions with various converters and file formats which we expect to succeed without
    issue.
    """
    run_test_conversion_with_library(specs.basic_tests)


def test_archive_convert():
    """Run a test of converting an archive of files
    """
    run_test_conversion_with_library(specs.archive_tests)


def test_archive_wrong_format():
    """Run a test that converting an archive but specifying the wrong input format will produce a warning
    """
    run_test_conversion_with_library(specs.archive_wrong_format_test)


def test_log_mode():
    """Test that the various log modes result in the expected log files being created
    """
    run_test_conversion_with_library(specs.log_mode_tests)


def test_stdout():
    """Test that the output is sent to stdout when requested
    """
    run_test_conversion_with_library(specs.stdout_test)


def test_quiet():
    """Test that quiet mode suppresses all output
    """
    run_test_conversion_with_library(specs.quiet_test)


def test_open_babel_warning():
    """Run a test that expected warnings from Open Babel are captured in the log
    """
    run_test_conversion_with_library(specs.open_babel_warning_test)


def test_exceed_output_file_size():
    """Run a test of the converter to ensure it reports an error properly if the output file size is too large
    """
    run_test_conversion_with_library(specs.max_size_test)


def test_invalid_converter():
    """Run a test of the converter to ensure it reports an error properly if an invalid converter is requested
    """
    run_test_conversion_with_library(specs.invalid_converter_test)


def test_quality_note():
    """Run a test of the converter on an `.xyz` to `.inchi` conversion which we expect to have warnings about data
    loss and extrapolation
    """
    run_test_conversion_with_library(specs.quality_note_test)


def test_cleanup():
    """Test that input files are deleted if requested
    """
    run_test_conversion_with_library(specs.cleanup_input_test)


def test_failed_conversion():
    """Run a test of the converter on a conversion we expect to fail
    """
    run_test_conversion_with_library(specs.failed_conversion_test)


def test_format_args():
    """Run a test that format args are processed correctly
    """
    run_test_conversion_with_library(specs.format_args_test)


def test_coord_gen():
    """Run a test that coordinate generation args are processed correctly
    """
    run_test_conversion_with_library(specs.coord_gen_test)


def test_envvars():
    """Test that setting appropriate envvars will set them for a file converter
    """

    test_file_size = 1234
    os.environ[const.MAX_FILESIZE_EV] = str(test_file_size)

    converter = OBFileConverter(filename="1NE6.mmcif",
                                to_format="pdb",
                                use_envvars=True,)
    assert math.isclose(converter.max_file_size, test_file_size*const.MEGABYTE)

    # And also check it isn't applied if we don't ask it to use envvars
    converter_no_ev = OBFileConverter(filename="1NE6.mmcif",
                                      to_format="pdb",
                                      use_envvars=False,)
    assert not math.isclose(converter_no_ev.max_file_size, test_file_size*const.MEGABYTE)
