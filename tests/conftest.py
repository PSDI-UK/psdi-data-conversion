"""@file tests/conftest.py

Created 2025-01-17 by Bryan Gillis.

Common fixtures used by unit tests
"""

import os
import pytest

TEST_DATA_LOC = os.path.abspath("./test_data")


@pytest.fixture
def test_data_loc():
    return TEST_DATA_LOC
