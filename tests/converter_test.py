"""@file psdi-data-conversion/tests/converter_test.py

Created 2024-12-15 by Bryan Gillis.

Unit tests of the converter class
"""

import os
import re
import pytest

from app import FILE_KEY, FILE_TO_UPLOAD_KEY
from psdi_data_conversion.log_utility import GLOBAL_LOG_FILENAME
from psdi_data_conversion.converter import (CONVERTER_OB, LOCAL_LOG_EXT, OUTPUT_LOG_EXT, STATUS_CODE_SIZE,
                                            FileConverter, FileConverterAbortException)

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


def test_mmcif_to_pdb(base_mock_form, tmp_upload_path, tmp_download_path):
    """Run a test of the converter on a straightforward `.mmcif` to `.pdb` conversion
    """

    # Remove the global log file if one exists
    try:
        os.remove(GLOBAL_LOG_FILENAME)
    except FileNotFoundError:
        pass

    # Save some variables from input we'll be using throughout this test
    source_filename = os.path.join(TEST_DATA_LOC, "1NE6.mmcif")
    files = get_mock_files(source_filename)
    filename = files[FILE_TO_UPLOAD_KEY].filename
    filename_base = os.path.splitext(filename)[0]
    to_format = base_mock_form["to"]

    test_converter = FileConverter(files=files,
                                   form=base_mock_form,
                                   file_to_convert=FILE_TO_UPLOAD_KEY,
                                   upload_dir=tmp_upload_path,
                                   download_dir=tmp_download_path)

    # Check that the input file now exists where we expect it to
    ex_input_filename = os.path.join(tmp_upload_path, filename)
    assert os.path.exists(ex_input_filename)

    test_converter.run()

    # Check that the input file has been properly deleted from the uploads directory
    ex_input_filename = os.path.join(tmp_upload_path, filename)
    assert not os.path.exists(ex_input_filename)

    # Check that the expected output file is found in the downloads directory
    ex_output_filename = os.path.join(tmp_download_path, f"{filename_base}.{to_format}")
    assert os.path.isfile(ex_output_filename)

    # Check that the logs are as we expect

    # Check that the global log file exists and is empty
    assert os.path.isfile(GLOBAL_LOG_FILENAME)
    global_log_text = open(GLOBAL_LOG_FILENAME).read()
    assert len(global_log_text) == 0

    # Check that the local log and output logs exist and contain expected information
    local_log_filename = os.path.join(tmp_download_path,
                                      f"{filename}-{filename_base}.{to_format}.{LOCAL_LOG_EXT}")
    assert os.path.isfile(local_log_filename)
    output_log_filename = os.path.join(tmp_download_path,
                                       f"{filename_base}.{OUTPUT_LOG_EXT}")
    assert os.path.isfile(local_log_filename)

    for filename in (local_log_filename, output_log_filename):
        log_text = open(filename).read()
        assert re.compile(r"File name:\s+"+filename_base).search(log_text)
        assert "Open Babel Warning" in log_text
        assert "Failed to kekulize aromatic bonds" in log_text
        assert log_text


def test_exceed_output_file_size(base_mock_form, tmp_upload_path, tmp_download_path):
    """Run a test of the converter to ensure it reports an error properly if the output file size is too large
    """

    # Remove the global log file if one exists
    try:
        os.remove(GLOBAL_LOG_FILENAME)
    except FileNotFoundError:
        pass

    source_filename = os.path.join(TEST_DATA_LOC, "1NE6.mmcif")
    files = get_mock_files(source_filename)

    test_converter = FileConverter(files=files,
                                   form=base_mock_form,
                                   file_to_convert=FILE_TO_UPLOAD_KEY,
                                   upload_dir=tmp_upload_path,
                                   download_dir=tmp_download_path,
                                   max_file_size=0)

    with pytest.raises(FileConverterAbortException) as esc_info:
        test_converter.run()
    assert esc_info.value.status_code == STATUS_CODE_SIZE

    # Check that the input file has been properly deleted from the uploads directory
    ex_input_filename = os.path.join(tmp_upload_path, files[FILE_TO_UPLOAD_KEY].filename)
    assert not os.path.exists(ex_input_filename)

    # Check that the expected output file is not found in the downloads directory
    ex_output_filename_base = os.path.splitext(files[FILE_TO_UPLOAD_KEY].filename)[0]
    ex_output_ext = base_mock_form["to"]
    ex_output_filename = os.path.join(tmp_download_path, f"{ex_output_filename_base}.{ex_output_ext}")
    assert not os.path.exists(ex_output_filename)
