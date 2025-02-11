"""@file tests/file_io_test.py

Created 2025-02-11 by Bryan Gillis.

Unit tests of functions in the `file_io` module
"""

import os

from psdi_data_conversion.file_io import is_archive, is_supported_archive, unpack_zip_or_tar
from psdi_data_conversion import constants as const

# Archive files prepared in the test data directory to be used for testing
ARCHIVE_FILENAME_TAR = "caffeine-smi.tar"
ARCHIVE_FILENAME_TARGZ = "caffeine-smi.tar.gz"
ARCHIVE_FILENAME_ZIP = "caffeine-smi.zip"


def test_is_archive():
    """Test the functions which check whether a file is an archive or not
    """

    # Check supported types return as expected
    for ext in const.L_SUPPORTED_ARCHIVE_EXTENSIONS:
        filename = f"file{ext}"
        assert is_archive(filename), ext
        assert is_supported_archive(filename), ext

    # Check unsupported types return as expected
    for ext in const.L_UNSUPPORTED_ARCHIVE_EXTENSIONS:
        filename = f"file{ext}"
        assert is_archive(filename), ext
        assert not is_supported_archive(filename), ext

    # Check non-archive type returns as expected
    assert not is_archive("foo.txt")
    assert not is_supported_archive("foo.txt")


def test_unpack_archive(test_data_loc, tmp_path_factory):
    """Test that archives can successfully be unpacked
    """

    for archive_filename in (ARCHIVE_FILENAME_TAR, ARCHIVE_FILENAME_TARGZ, ARCHIVE_FILENAME_ZIP):

        # Get a temporary directory to work with for each different archive file we test
        extract_dir = tmp_path_factory.mktemp("unpack-test")

        # Check that we can extract each archive without error
        qual_archive_filename = os.path.join(test_data_loc, archive_filename)

        l_filenames = unpack_zip_or_tar(qual_archive_filename, extract_dir=extract_dir)

        # Check that files were successfully extracted
        assert l_filenames, archive_filename
        for filename in l_filenames:
            assert os.path.isfile(filename), (archive_filename, filename)
