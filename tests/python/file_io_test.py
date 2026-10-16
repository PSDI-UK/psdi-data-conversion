"""@file tests/file_io_test.py

Created 2025-02-11 by Bryan Gillis.

Unit tests of functions in the `file_io` module
"""

import os
from itertools import product

import pytest

from psdi_data_conversion import constants as const
from psdi_data_conversion.file_io import (is_archive, is_supported_archive, pack_zip_or_tar, split_archive_ext,
                                          unpack_zip_or_tar)
from psdi_data_conversion.testing.utils import get_input_test_data_loc

# Archive files prepared in the test data directory to be used for testing
ARCHIVE_FILENAME_TAR = "caffeine-smi.tar"
ARCHIVE_FILENAME_TARGZ = "caffeine-smi.tar.gz"
ARCHIVE_FILENAME_ZIP = "caffeine-smi.zip"


@pytest.mark.parametrize("ext", const.D_SUPPORTED_ARCHIVE_FORMATS)
def test_is_archive_supported_types(ext, subtests):
    """Test that supported archive types are correctly reported as such"""
    filename = f"file{ext}"
    with subtests.test("correctly determined to be an archive"):
        assert is_archive(filename)
    with subtests.test("correctly determined to be supported"):
        assert is_supported_archive(filename)


@pytest.mark.parametrize("ext", const.L_UNSUPPORTED_ARCHIVE_EXTENSIONS)
def test_is_archive_unsupported_types(ext, subtests):
    """Test that unsupported archive types are correctly reported as such"""
    filename = f"file{ext}"
    with subtests.test("correctly determined to be an archive"):
        assert is_archive(filename)
    with subtests.test("correctly determined to be not supported"):
        assert not is_supported_archive(filename)


@pytest.mark.parametrize("ext", (".txt", ".zip.bak"))
def test_is_archive_invalid_types(ext, subtests):
    """Test that invalid archive types are correctly reported as such"""
    filename = f"file{ext}"
    with subtests.test("correctly determined to not be an archive"):
        assert not is_archive(filename)
    with subtests.test("correctly determined to be not supported"):
        assert not is_supported_archive(filename)


l_bases = ["foo", "foo.txt", "footar", "foo-1.0.3"]
l_exts = [".zip", ".tar", ".tar.gz", ".gz"]


@pytest.mark.parametrize("base, ext", list(product(l_bases, l_exts)))
def test_split_archive_ext(base, ext):
    """Test the function to split a filename with a possible archive extension"""

    filename = base+ext
    assert split_archive_ext(filename) == (base, ext)


@pytest.mark.parametrize("archive_filename", (ARCHIVE_FILENAME_TAR, ARCHIVE_FILENAME_TARGZ, ARCHIVE_FILENAME_ZIP))
def test_archive(archive_filename, tmp_path_factory, subtests):
    """Test that archives can successfully be unpacked"""

    # Get a temporary directory to work with
    extract_dir = tmp_path_factory.mktemp("unpack-test")

    qual_archive_filename = os.path.join(get_input_test_data_loc(), archive_filename)

    # Check that we can extract each archive without error
    l_qualified_filenames = unpack_zip_or_tar(qual_archive_filename, extract_dir=extract_dir)

    # Check that files were successfully extracted
    with subtests.test("list of files extracted is provided"):
        assert l_qualified_filenames, archive_filename
    for filename in l_qualified_filenames:
        with subtests.test("extracted file exists", filename=filename):
            assert os.path.isfile(filename)

    # Test that we can repack the files into a new archive in various ways
    new_archive_filename = f"new-{archive_filename}"
    pack_dir = tmp_path_factory.mktemp("pack-test")
    qual_new_archive_filename = os.path.join(pack_dir, new_archive_filename)

    # Test using the fully-qualified filenames we were given when the files were unpacked
    qual_new_archive_filename_returned = pack_zip_or_tar(qual_new_archive_filename,
                                                         l_filenames=l_qualified_filenames)
    with subtests.test("packed file exists"):
        assert os.path.isfile(qual_new_archive_filename_returned)

    # Check that input files were not deleted
    for filename in l_qualified_filenames:
        with subtests.test("input file not deleted", filename=filename):
            assert os.path.isfile(filename)

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
    with subtests.test("packed file exists when using relative filenames"):
        assert os.path.isfile(qual_new_archive_filename)

    # Check that input files were deleted
    for filename in l_qualified_filenames:
        with subtests.test("input file deleted", filename=filename):
            assert not os.path.exists(filename), (archive_filename, filename)


def test_unpack_unsupported():
    """Test that we can't unpack an unsupported archive"""
    with pytest.raises(ValueError, match="unsupported"):
        unpack_zip_or_tar("foo.rar")


def test_unpack_invalid():
    """Test that we can't unpack an invalid archive"""
    with pytest.raises(ValueError, match="valid"):
        unpack_zip_or_tar("foo.txt")
