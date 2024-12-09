"""@file psdi-data-conversion/tests/logging_test.py

Created 2024-12-09 by Bryan Gillis.

Tests of functions relating to logging
"""

import os
import re
import time

from psdi_data_conversion import logging
from app import get_date, get_date_time, get_time


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

    date_str_1 = get_date()
    time_str_1 = get_time()
    datetime_str_1 = get_date_time()

    assert date_re.match(date_str_1)
    assert time_re.match(time_str_1)
    assert datetime_re.match(datetime_str_1)

    # Test that the time changes after a second, and is still in the correct format
    time.sleep(1)
    time_str_2 = get_time()
    datetime_str_2 = get_date_time()

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
    fn_logger = logging.getLogger("fn.test", test_filename)
    assert fn_logger.getGlobalFilename() == os.path.abspath(logging.GLOBAL_ERROR_LOG)
    assert fn_logger.getLocalFilename() == os.path.abspath(test_filename)
