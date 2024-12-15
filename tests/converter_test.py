"""@file psdi-data-conversion/tests/converter_test.py

Created 2024-12-15 by Bryan Gillis.

Unit tests of the converter class
"""

import os
import pytest

from app import FILE_KEY, FILE_TO_UPLOAD_KEY
from psdi_data_conversion.converter import CONVERTER_OB, STATUS_CODE_SIZE, FileConverter

TEST_DATA_LOC = os.path.abspath("./test_data")


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
    """A fixture providing a default `form` object which can be used to instantiate a converter
    """
    return {'token': '1041c0a661d118d5f28e7c6830375dd0',
            'converter': CONVERTER_OB,
            'from': 'mmcif',
            'to': 'pdb',
            'from_full': 'mmcif: Macromolecular Crystallographic Info',
            'to_full': 'pdb: Protein Data Bank',
            'success': 'good',
            'from_flags': '',
            'to_flags': '',
            'from_arg_flags': '',
            'from_args': '',
            'to_arg_flags': '',
            'to_args': '',
            'coordinates': 'neither',
            'coordOption': 'medium',
            'upload_file': 'true'}


@pytest.fixture()
def tmp_upload_path(tmp_path):
    """Fixture providing a temporary directory for upload files
    """
    upload_path = os.path.join(tmp_path, "uploads")
    os.makedirs(upload_path, exist_ok=True)
    return upload_path


@pytest.fixture()
def tmp_download_path(tmp_path):
    """Fixture providing a temporary directory for download files
    """
    download_path = os.path.join(tmp_path, "downloads")
    os.makedirs(download_path, exist_ok=True)
    return download_path


p_status_code = [0]


def abort_raise(status_code):
    """Callback for aborting during a test, which re-raises any raised exceptions
    """
    try:
        raise
    except RuntimeError as e:
        if "No active exception to reraise" not in str(e):
            raise
        p_status_code[0] = status_code


def test_mmcif_to_pdb(base_mock_form, tmp_upload_path, tmp_download_path):
    """Run a test of the converter on a straightforward `.mmcif` to `.pdb` conversion
    """

    source_filename = os.path.join(TEST_DATA_LOC, "1NE6.mmcif")

    # Reset the stored status code to zero - will be changed if something goes wrong
    p_status_code[0] = 0

    test_converter = FileConverter(files=get_mock_files(source_filename),
                                   form=base_mock_form,
                                   file_to_convert=FILE_TO_UPLOAD_KEY,
                                   abort_callback=abort_raise,
                                   upload_dir=tmp_upload_path,
                                   download_dir=tmp_download_path)

    test_converter.run()

    assert p_status_code[0] == 0, f"Converter aborted with status code {p_status_code[0]}"


def test_exceed_output_file_size(base_mock_form, tmp_upload_path, tmp_download_path):
    """Run a test of the converter to ensure it reports an error properly if the output file size is too large
    """

    source_filename = os.path.join(TEST_DATA_LOC, "1NE6.mmcif")

    # Reset the stored status code to zero - we want to make sure it gets set to 413 correctly
    p_status_code[0] = 0

    test_converter = FileConverter(files=get_mock_files(source_filename),
                                   form=base_mock_form,
                                   file_to_convert=FILE_TO_UPLOAD_KEY,
                                   abort_callback=abort_raise,
                                   upload_dir=tmp_upload_path,
                                   download_dir=tmp_download_path,
                                   max_file_size=0)

    test_converter.run()

    assert p_status_code[0] == STATUS_CODE_SIZE, f"Converter did not abort with status code {STATUS_CODE_SIZE}"

    test_converter.run()
