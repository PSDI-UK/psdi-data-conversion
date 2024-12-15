"""@file psdi-data-conversion/tests/converter_test.py

Created 2024-12-15 by Bryan Gillis.

Unit tests of the converter class
"""

import os
import pytest
from app import FILE_KEY, FILE_TO_UPLOAD_KEY

TEST_DATA_LOC = "./test_data"


class MockFileStorage:
    """Mock version of the `FileStorage` class which provides the needed functionality in a way convenient for unit
    tests.
    """
    filename: str | None = None
    source_filename: str | None = None

    def __init__(self, source_filename):
        self.source_filename = source_filename
        self.filename = os.path.split(self.source_filename)[1]

    def save(self, dest_filename):
        """To speed things up, symlink the file instead of creating a copy
        """

        # Silently make sure the destination directory exists
        os.makedirs(os.path.split(dest_filename)[0], exist_ok=True)

        os.symlink(self.source_filename, dest_filename)


def get_mock_files(source_filename):
    """Convenience function for unit test to get a mock `files` dict to pass as an argument to initializing a converter
    """
    mock_file_storage = MockFileStorage(source_filename)
    return {FILE_KEY: mock_file_storage,
            FILE_TO_UPLOAD_KEY: mock_file_storage,
            }


@pytest.fixture()
def base_mock_form():
    return {'to': None,
            'from': None,
            'converter': None,
            'from_flags': None,
            'to_flags': None,
            'from_arg_flags': None,
            'to_arg_flags': None,
            'from_args': None,
            'to_args': None,
            'coordinates': None,
            'coordOption': None,
            'from_full': None,
            'to_full': None,
            'success': None, }
