"""@file psdi-data-conversion/tests/logging_test.py

Created 2024-12-09 by Bryan Gillis.

Tests of functions relating to logging
"""

import os
import re
import time

from psdi_data_conversion import logging


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

    date_str_1 = logging.get_date()
    time_str_1 = logging.get_time()
    datetime_str_1 = logging.get_date_time()

    assert date_re.match(date_str_1)
    assert time_re.match(time_str_1)
    assert datetime_re.match(datetime_str_1)

    # Test that the time changes after a second, and is still in the correct format
    time.sleep(1)
    time_str_2 = logging.get_time()
    datetime_str_2 = logging.get_date_time()

    assert time_re.match(time_str_2)
    assert datetime_re.match(datetime_str_2)

    assert time_str_2 != time_str_1
    assert datetime_str_2 != datetime_str_1


def test_get_logger():
    """Tests of `logging.getLogger`
    """
    # Get a logger to test with
    logger = logging.getLogger("test")

    # Test getting a second logger with the same name returns the same as the first
    same_name_logger = logging.getLogger("test")
    assert same_name_logger is logger

    # Test getting a logger with a different name returns a different logger
    diff_name_logger = logging.getLogger("not.test")
    assert diff_name_logger is not logger

    # Test that the filenames are as expected
    assert logger.getGlobalFilename() == os.path.abspath(logging.GLOBAL_ERROR_LOG)
    assert logger.getLocalFilename() is None

    test_filename = "./static/downloads/local_error_log.txt"
    os.makedirs("./static/downloads/", exist_ok=True)
    fn_logger = logging.getLogger("fn.test", test_filename)
    assert fn_logger.getGlobalFilename() == os.path.abspath(logging.GLOBAL_ERROR_LOG)
    assert fn_logger.getLocalFilename() == os.path.abspath(test_filename)


def test_logging(caplog):
    """Test that logging works as expected
    """

    test_filename = "./static/downloads/local_error_log.txt"
    os.makedirs("./static/downloads/", exist_ok=True)

    # Delete any existing error logs
    if os.path.isfile(logging.GLOBAL_ERROR_LOG):
        os.remove(logging.GLOBAL_ERROR_LOG)
    if os.path.isfile(test_filename):
        os.remove(test_filename)

    # Create a logger to work with
    logger = logging.getLogger("logging-test", test_filename)

    # Try logging a few messages at different levels
    debug_msg = "FINDME_DEBUG"
    info_msg = "FINDME_INFO"
    warning_msg = "FINDME_WARNING"
    logger.debug(debug_msg)
    logger.info(info_msg)
    logger.warning(warning_msg)

    # Open the files and check that only the expected messages are present - by default, the global log will log
    # at WARNING level and above, and the local log will log at INFO level and above

    log_content = caplog.text

    assert not re.compile("test:" + r"\S*\s*\S*" + debug_msg).search(log_content)
    assert not re.compile("test:" + r"\S*\s*\S*" + info_msg).search(log_content)
    assert re.compile("test:" + r"\S*\s*\S*" + warning_msg).search(log_content)

    assert not re.compile("test-local:" + r"\S*\s*\S*" + debug_msg).search(log_content)
    assert re.compile("test-local:" + r"\S*\s*\S*" + info_msg).search(log_content)
    assert re.compile("test-local:" + r"\S*\s*\S*" + warning_msg).search(log_content)
