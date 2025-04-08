"""
# constants.py

Constants related to unit testing
"""

import os

# Locations relative to the root directory of the project - to ensure the files are found, tests should chdir to this
# directory before searching for files

TEST_DATA_LOC_IN_PROJECT = "./test_data"

INPUT_TEST_DATA_LOC_IN_PROJECT = TEST_DATA_LOC_IN_PROJECT
OUTPUT_TEST_DATA_LOC_IN_PROJECT = os.path.join(TEST_DATA_LOC_IN_PROJECT, "output")
