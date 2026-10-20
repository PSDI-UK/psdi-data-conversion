"""@file tests/file_io_test.py

Created 2025-02-11 by Bryan Gillis.

Unit tests of functions in the `file_io` module
"""

from itertools import product
import os

import pytest

from psdi_data_conversion.file_io import (is_archive, is_supported_archive, pack_zip_or_tar, split_archive_ext,
                                          unpack_zip_or_tar)
from psdi_data_conversion import constants as const
from psdi_data_conversion.testing.utils import get_input_test_data_loc

# Archive files prepared in the test data directory to be used for testing
ARCHIVE_FILENAME_TAR = "caffeine-smi.tar"
ARCHIVE_FILENAME_TARGZ = "caffeine-smi.tar.gz"
ARCHIVE_FILENAME_ZIP = "caffeine-smi.zip"


def test_is_archive():
    """Test the functions which check whether a file is an archive or not
    """

    # Check supported types return as expected
    for ext in const.D_SUPPORTED_ARCHIVE_FORMATS:
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


def test_split_archive_ext():
    """Test the function to split a filename with a possible archive extension
    """

    l_bases = ["foo", "foo.txt", "footar", "foo-1.0.3"]
    l_exts = [".zip", ".tar", ".tar.gz", ".gz"]

    for base, ext in product(l_bases, l_exts):
        filename = base+ext
        assert split_archive_ext(filename) == (base, ext)


def test_archive(tmp_path_factory):
    """Test that archives can successfully be unpacked
    """

    input_test_data_loc = get_input_test_data_loc()

    for archive_filename in (ARCHIVE_FILENAME_TAR, ARCHIVE_FILENAME_TARGZ, ARCHIVE_FILENAME_ZIP):

        # Get a temporary directory to work with for each different archive file we test
        extract_dir = tmp_path_factory.mktemp("unpack-test")

        # Check that we can extract each archive without error
        qual_archive_filename = os.path.join(input_test_data_loc, archive_filename)

        l_qualified_filenames = unpack_zip_or_tar(qual_archive_filename, extract_dir=extract_dir)

        # Check that files were successfully extracted
        assert l_qualified_filenames, archive_filename
        for filename in l_qualified_filenames:
            assert os.path.isfile(filename), (archive_filename, filename)

        # Test that we can repack the files into a new archive in various ways
        new_archive_filename = f"new-{archive_filename}"
        pack_dir = tmp_path_factory.mktemp("pack-test")
        qual_new_archive_filename = os.path.join(pack_dir, new_archive_filename)

        # Test using the fully-qualified filenames we were given when the files were unpacked
        qual_new_archive_filename_returned = pack_zip_or_tar(qual_new_archive_filename,
                                                             l_filenames=l_qualified_filenames)
        assert os.path.isfile(qual_new_archive_filename_returned)

        # Check that input files were not deleted
        for filename in l_qualified_filenames:
            assert os.path.isfile(filename), (archive_filename, filename)

        # Now try using relative filenames, just the root name for the archive, and clean up afterwards
        l_filenames = [os.path.split(x)[1] for x in l_qualified_filenames]
        os.remove(qual_new_archive_filename)

        new_archive_filename_base, ext = os.path.splitext(qual_new_archive_filename)
        if new_archive_filename_base.endswith(".tar"):
            new_archive_filename_base, pre_ext = os.path.splitext(new_archive_filename_base)
            ext = pre_ext + ext

        pack_zip_or_tar(new_archive_filename_base,
                        l_filenames=l_filenames,
                        source_dir=extract_dir,
                        archive_format=ext,
                        cleanup=True)
        assert os.path.isfile(qual_new_archive_filename)

        # Check that input files were deleted
        for filename in l_qualified_filenames:
            assert not os.path.exists(filename), (archive_filename, filename)

    # Check that we get expected exceptions for unsupported archives
    with pytest.raises(ValueError, match="unsupported"):
        unpack_zip_or_tar("foo.rar")
    with pytest.raises(ValueError, match="valid"):
        unpack_zip_or_tar("foo.txt")
