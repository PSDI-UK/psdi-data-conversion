"""@file psdi-data-conversion/tests/logging_test.py

Created 2024-12-09 by Bryan Gillis.

Tests of functions relating to logging
"""

import os
import re
import time
import logging

from psdi_data_conversion import log_utility


def test_date_time():
    """Tests of getting the current date and time
    """

    # Use a regex match to test that the date and time are in the right format
    date_re_raw = r"\d{4}-[0-1]\d-[0-3]\d"
    time_re_raw = r"[0-2]\d:[0-5]\d:[0-5]\d"
    datetime_re_raw = f"{date_re_raw} {time_re_raw}"

    date_re = re.compile(date_re_raw)
    time_re = re.compile(time_re_raw)
    datetime_re = re.compile(datetime_re_raw)

    date_str_1 = log_utility.get_date()
    time_str_1 = log_utility.get_time()
    datetime_str_1 = log_utility.get_date_time()

    assert date_re.match(date_str_1)
    assert time_re.match(time_str_1)
    assert datetime_re.match(datetime_str_1)

    # Test that the time changes after a second, and is still in the correct format
    time.sleep(1)
    time_str_2 = log_utility.get_time()
    datetime_str_2 = log_utility.get_date_time()

    assert time_re.match(time_str_2)
    assert datetime_re.match(datetime_str_2)

    assert time_str_2 != time_str_1
    assert datetime_str_2 != datetime_str_1


def test_get_logger(tmp_path):
    """Tests of `log_utility.getDataConversionLogger`
    """
    # Get a logger to test with
    logger = log_utility.getDataConversionLogger("test")

    # Test getting a second logger with the same name returns the same as the first
    same_name_logger = log_utility.getDataConversionLogger("test")
    assert same_name_logger is logger

    # Test getting a logger with a different name returns a different logger
    diff_name_logger = log_utility.getDataConversionLogger("not.test")
    assert diff_name_logger is not logger

    # Test that a logger without a name provided will also differ
    no_name_logger = log_utility.getDataConversionLogger()
    assert no_name_logger is not logger

    # Test that the filenames are as expected
    test_filename = os.path.join(tmp_path, "log.txt")
    test_level = logging.CRITICAL
    fn_logger = log_utility.getDataConversionLogger("fn-test", local_log_file=test_filename,
                                                    local_logger_level=test_level)

    # Search through the logger's handlers to get all files it logs to and at what levels
    l_files_and_levels = []
    for handler in fn_logger.handlers:
        if isinstance(handler, logging.FileHandler):
            l_files_and_levels.append((handler.baseFilename, handler.level))
    assert (os.path.abspath(log_utility.GLOBAL_LOG_FILENAME), log_utility.GLOBAL_LOGGER_LEVEL) in l_files_and_levels
    assert (os.path.abspath(test_filename), test_level) in l_files_and_levels


def test_logging(tmp_path):
    """Test that logging works as expected
    """

    test_filename = os.path.join(tmp_path, "log.txt")

    # Delete any existing error logs
    if os.path.isfile(log_utility.GLOBAL_LOG_FILENAME):
        os.remove(log_utility.GLOBAL_LOG_FILENAME)
    if os.path.isfile(test_filename):
        os.remove(test_filename)

    logger_name = "log_utility-test"

    # Create a logger to work with
    logger = log_utility.getDataConversionLogger(logger_name, test_filename)
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

    with open(log_utility.GLOBAL_LOG_FILENAME, "r") as fi:
        global_log_content = fi.read()

        assert debug_msg not in global_log_content
        assert info_msg not in global_log_content
        assert error_msg in global_log_content

    with open(test_filename, "r") as fi:
        local_log_content = fi.read()

        assert debug_msg not in local_log_content
        assert info_msg in local_log_content
        assert error_msg in local_log_content
