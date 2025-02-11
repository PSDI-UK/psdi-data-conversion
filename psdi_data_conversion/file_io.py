"""@file psdi_data_conversion/file_io.py

Created 2025-02-11 by Bryan Gillis.

Functions and classes related to general filesystem input/output
"""

import os
from shutil import unpack_archive

from psdi_data_conversion import constants as const


def unpack_zip_or_tar(archive_filename: str) -> list[str]:
    """Unpack a zip or tar archive into a temporary directory and return a list of the extracted files

    Parameters
    ----------
    archive_filename : str
        Filename of the archive to unpack, either relative or fully-qualified

    Returns
    -------
    list[str]
        List of the fully-qualified paths to the extracted files
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

    unpack_archive(qual_archive_filename, **unpack_kwargs)
