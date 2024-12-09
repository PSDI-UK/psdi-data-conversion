"""@file psdi-data-conversion/tests/logging_test.py

Created 2024-12-09 by Bryan Gillis.

Tests of functions relating to logging
"""

import re
import time

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
