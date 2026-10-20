"""@file psdi-data-conversion/tests/logging_test.py

Created 2024-12-09 by Bryan Gillis.

Tests of functions relating to logging
"""

import os
import re
import time
import logging

from psdi_data_conversion import log_utility
from psdi_data_conversion import constants as const


def test_date_time():
    """Tests of getting the current date and time
    """

    # Use a regex match to test that the date and time are in the right format
    date_re = re.compile(const.DATE_RE_RAW)
    time_re = re.compile(const.TIME_RE_RAW)
    datetime_re = re.compile(const.DATETIME_RE_RAW)

    date_str_1 = log_utility.get_date()
    time_str_1 = log_utility.get_time()
    datetime_str_1 = log_utility.get_date_time()

    assert date_re.match(date_str_1)
    assert time_re.match(time_str_1)
    assert datetime_re.match(datetime_str_1)

    # Test that the time changes after a second, and is still in the correct format
    time.sleep(1.2)
    time_str_2 = log_utility.get_time()
    datetime_str_2 = log_utility.get_date_time()

    assert time_re.match(time_str_2)
    assert datetime_re.match(datetime_str_2)

    assert time_str_2 != time_str_1
    assert datetime_str_2 != datetime_str_1


def test_setup_logger(tmp_path):
    """Tests of `log_utility.setUpDataConversionLogger`
    """
    # Get a logger to test with
    logger = log_utility.set_up_data_conversion_logger("test")

    # Test getting a second logger with the same name does not return the same as the first
    same_name_logger = log_utility.set_up_data_conversion_logger("test")
    assert same_name_logger is not logger

    # Test getting a logger with a different name returns a different logger
    diff_name_logger = log_utility.set_up_data_conversion_logger("not.test")
    assert diff_name_logger is not logger

    # Test that a logger without a name provided will also differ
    no_name_logger = log_utility.set_up_data_conversion_logger()
    assert no_name_logger is not logger

    # Test that the filenames are as expected
    test_log_filename = os.path.join(tmp_path, "log.txt")
    test_log_level = logging.WARN
    test_error_filename = os.path.join(tmp_path, "err.txt")
    test_error_level = logging.CRITICAL
    fn_logger = log_utility.set_up_data_conversion_logger("fn-test",
                                                          local_log_file=test_log_filename,
                                                          local_logger_level=test_log_level,
                                                          extra_loggers=[(test_error_filename,
                                                                          test_error_level,
                                                                          False, "a")])

    # Search through the logger's handlers to get all files it logs to and at what levels
    l_files_and_levels = []
    for handler in fn_logger.handlers:
        if isinstance(handler, logging.FileHandler):
            l_files_and_levels.append((handler.baseFilename, handler.level))
    assert (os.path.abspath(const.GLOBAL_LOG_FILENAME),
            const.GLOBAL_LOGGER_LEVEL) in l_files_and_levels
    assert (os.path.abspath(test_log_filename), test_log_level) in l_files_and_levels
    assert (os.path.abspath(test_error_filename), test_error_level) in l_files_and_levels


def test_logging(tmp_path):
    """Test that logging works as expected
    """

    test_filename = os.path.join(tmp_path, "log.txt")

    # Delete any existing error logs
    if os.path.isfile(const.GLOBAL_LOG_FILENAME):
        os.remove(const.GLOBAL_LOG_FILENAME)
    if os.path.isfile(test_filename):
        os.remove(test_filename)

    logger_name = "log_utility-test"

    # Create a logger to work with
    logger = log_utility.set_up_data_conversion_logger(logger_name, test_filename)
    logger.setLevel(logging.INFO)

    # Try logging a few messages at different levels
    debug_msg = "FINDME_DEBUG"
    info_msg = "FINDME_INFO"
    error_msg = "FINDME_ERROR"
    logger.debug(debug_msg)
    logger.info(info_msg)
    logger.error(error_msg)

    # Open the files and check that only the expected messages are present - by default, the global log will log
    # at ERROR level and above, and the local log will log at INFO level and above

    with open(const.GLOBAL_LOG_FILENAME, "r") as fi:
        global_log_content = fi.read()

        assert debug_msg not in global_log_content
        assert info_msg not in global_log_content
        assert error_msg in global_log_content

    with open(test_filename, "r") as fi:
        local_log_content = fi.read()

        assert debug_msg not in local_log_content
        assert info_msg in local_log_content
        assert error_msg in local_log_content
