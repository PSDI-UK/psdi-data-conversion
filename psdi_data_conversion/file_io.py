"""@file psdi_data_conversion/file_io.py

Created 2025-02-11 by Bryan Gillis.

Functions and classes related to general filesystem input/output
"""

import glob
import os
from shutil import copyfile, make_archive, unpack_archive
from tempfile import TemporaryDirectory

from psdi_data_conversion import constants as const


def is_archive(filename: str) -> bool:
    """Uses a file's extension to check if it's an archive or not
    """
    return any([filename.endswith(x) for x in const.L_ALL_ARCHIVE_EXTENSIONS])


def is_supported_archive(filename: str) -> bool:
    """Uses a file's extension to check if it's an archive of a supported type or not
    """
    return any([filename.endswith(x) for x in const.D_SUPPORTED_ARCHIVE_FORMATS])


def split_archive_ext(filename: str) -> tuple[str, str]:
    """Splits a file into a base and an extension, with handling for compound .tar.* extensions
    """
    base, ext = os.path.splitext(filename)
    if base.endswith(const.TAR_EXTENSION):
        base, pre_ext = os.path.splitext(base)
        ext = pre_ext+ext
    return base, ext


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

    elif any([qual_archive_filename.endswith(x) for x in const.D_ZIP_FORMATS]):

        # Zip types don't support the "filter" kwarg, but use similar security measures by default. This may prompt
        # a warning, which can be ignored
        pass

    elif any([qual_archive_filename.endswith(x) for x in const.D_TAR_FORMATS]):

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
    l_qual_new_files = [os.path.join(extract_dir, x) for x in s_new_files]
    l_new_globs = [glob.glob(x) if os.path.isfile(x)
                   else glob.glob(os.path.join(x, "**"))
                   for x in l_qual_new_files]

    # This gives us a list of globs (individual files are set up as globs for consistency), so we unpack to a single
    # list with nested list comprehension
    l_new_files = [x for l_glob_files in l_new_globs for x in l_glob_files]

    # Sort the file list for consistency in output
    l_new_files.sort(key=lambda s: s.lower())

    return l_new_files


def pack_zip_or_tar(archive_filename: str,
                    l_filenames: list[str],
                    archive_format: str | None = None,
                    source_dir: str = ".",
                    cleanup=False) -> str:
    """_summary_

    Parameters
    ----------
    archive_filename : str
        The desired name of the output archive to create, either fully-qualified or relative to the current directory
    l_filenames : list[str]
        List of files to be archived, either fully-qualified or relative to `source_dir`. If provided fully-qualified,
        they will be placed in the root directory of the archive
    source_dir : str, optional
        Path to directory containing the files to be archived (default current directory). If filenames are provided
        fully-qualified, this is ignored
    cleanup : bool, optional
        If True, source files will be deleted after the archive is successfully created

    Returns
    -------
    str
        The name of the created archive file

    Raises
    ------
    ValueError
        If `archive_filename` is not of a valid format
    FileNotFoundError
        If one of the listed files does not exist
    """

    if not archive_format and not is_supported_archive(archive_filename):
        raise ValueError(f"Desired archive filename '{archive_filename}' is not of a supported type. Supported types "
                         f"are: {const.D_SUPPORTED_ARCHIVE_FORMATS.keys()}")

    # It's supported, so determine the specific format, and provide it and the base of the filename in the forms that
    # `make_archive` wants
    if archive_format is None:
        for _ext, _format in const.D_SUPPORTED_ARCHIVE_FORMATS.items():
            if archive_filename.endswith(_ext):
                archive_format = _format
                archive_root_filename = split_archive_ext(archive_filename)[0]
                break
        else:
            raise AssertionError("Invalid execution path entered - filename wasn't found with a valid archive "
                                 "extension, but it did pass the `is_supported_archive` check")
    else:
        archive_root_filename = archive_filename

        # Check that the provided archive format is valid, and add the appropriate extension to the filename
        archive_extension: str | None = None
        for _ext, _format in const.D_SUPPORTED_ARCHIVE_FORMATS.items():
            if archive_format == _ext:
                # Extension was provided instead of the format; we can work with that
                archive_extension = archive_format
                archive_format = _format
                break
            elif archive_format == _format:
                archive_extension = _ext
                break
        if archive_extension is None:
            raise ValueError(f"Invalid archive format '{archive_format}'. Valid formats are: "
                             f"{const.D_SUPPORTED_ARCHIVE_FORMATS.keys()}")

        # Check if the root filename already contained the extension so we don't add it again, and strip it from the
        # root
        if archive_root_filename.endswith(archive_extension):
            archive_filename = archive_root_filename
            archive_root_filename = archive_root_filename[:-len(archive_extension)]
        else:
            archive_filename = archive_root_filename+archive_extension

    with TemporaryDirectory() as root_dir:

        # Copy all files from the source dir to the root dir, which is what will be packed

        l_files_to_cleanup: list[str] = []

        for filename in l_filenames:

            # Check if the filename is fully-qualified, and copy it from wherever it's found
            if os.path.isfile(filename):
                copyfile(filename, os.path.join(root_dir, os.path.basename(filename)))
                l_files_to_cleanup.append(filename)
                continue

            qualified_filename = os.path.join(source_dir, filename)
            if os.path.isfile(qualified_filename):
                copyfile(qualified_filename, os.path.join(root_dir, os.path.basename(filename)))
                l_files_to_cleanup.append(qualified_filename)
            else:
                raise FileNotFoundError(f"File '{filename}' could not be found, either fully-qualified or relative to "
                                        f"{source_dir}")

        make_archive(archive_root_filename,
                     format=archive_format,
                     root_dir=root_dir)

    if cleanup:
        for filename in l_files_to_cleanup:
            try:
                os.remove(filename)
            except Exception:
                pass

    # Return the name of the created file
    return archive_filename
