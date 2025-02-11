"""@file psdi_data_conversion/file_io.py

Created 2025-02-11 by Bryan Gillis.

Functions and classes related to general filesystem input/output
"""

import glob
import os
from shutil import unpack_archive

from psdi_data_conversion import constants as const


def unpack_zip_or_tar(archive_filename: str,
                      extract_dir: str = ".") -> list[str]:
    """Unpack a zip or tar archive into a temporary directory and return a list of the extracted files

    Parameters
    ----------
    archive_filename : str
        Filename of the archive to unpack, either relative or fully-qualified
    extract_dir : str
        The directory to extract the contents of the archive to. By default, the current working directory will
        be used

    Returns
    -------
    list[str]
        List of the fully-qualified paths to the extracted files. This is determined by checking the directory contents
        before and after extraction, so is NOT thread-safe, unless it is otherwise ensured e.g. by using a unique
        temporary directory for each thread
    """

    qual_archive_filename = os.path.realpath(archive_filename)

    # Determine if the file is of a known (un)supported archive type, and if it is, whether it's a zip or tar, and
    # set up arguments appropriately to ensure security
    unpack_kwargs: dict[str, str] = {}

    if any([qual_archive_filename.endswith(x) for x in const.L_UNSUPPORTED_ARCHIVE_EXTENSIONS]):

        raise ValueError(f"The archive file '{qual_archive_filename}' is of an unsupported archive type")

    elif any([qual_archive_filename.endswith(x) for x in const.L_ZIP_EXTENSIONS]):

        # Zip types don't support the "filter" kwarg, but use similar security measures by default. This may prompt
        # a warning, which can be ignored
        pass

    elif any([qual_archive_filename.endswith(x) for x in const.L_TAR_EXTENSIONS]):

        # Tar types need to set up the "filter" argument to ensure no files are unpacked outside the base directory
        unpack_kwargs["filter"] = "data"

    else:

        raise ValueError(f"The archive file '{qual_archive_filename}' is not recognised as a valid archive type")

    # To determine the names of extracted files, we call `os.listdir` before and after unpacking and look for the new
    # elements

    s_dir_before = set(os.listdir(extract_dir))
    unpack_archive(qual_archive_filename, extract_dir=extract_dir, **unpack_kwargs)
    s_dir_after = set(os.listdir(extract_dir))

    # Get the new files, and in case they're in a directory, use glob to get their contents
    s_new_files = s_dir_after.difference(s_dir_before)
    l_new_globs = [glob.glob(x) if os.path.isfile(x)
                   else glob.glob(os.path.join(x, "**"))
                   for x in s_new_files]

    # This gives us a list of globs (individual files are set up as globs for consistency), so we unpack to a single
    # list with nested list comprehension
    l_new_files = [x for l_glob_files in l_new_globs for x in l_glob_files]

    return l_new_files
