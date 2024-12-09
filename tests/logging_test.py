"""@file psdi-data-conversion/tests/logging_test.py

Created 2024-12-09 by Bryan Gillis.

Tests of functions relating to logging
"""

import re
import time

from app import get_date_time


def test_date_time():
    """Tests of getting the current date and time
    """

    # Use a regex match to test that the date is in the right format
    ex_re = re.compile(r"\d{4}-[0-1]\d-[0-3]\d")

    dt_str_1 = get_date_time()
    assert ex_re.match(dt_str_1)

    # Test that the time changes after a second, and is still in the correct format
    time.sleep(1)
    dt_str_2 = get_date_time()
    assert ex_re.match(dt_str_2)
    assert dt_str_2 != dt_str_1
