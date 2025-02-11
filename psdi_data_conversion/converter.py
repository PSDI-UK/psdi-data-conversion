"""@file psdi-data-conversion/psdi_data_conversion/converter.py

Created 2024-12-10 by Bryan Gillis.

Class and functions to perform file conversion
"""

import os
import importlib
import sys
from tempfile import TemporaryDirectory
import traceback
from typing import NamedTuple
from psdi_data_conversion import constants as const
from psdi_data_conversion.converters import base

import glob

from psdi_data_conversion.file_io import is_archive, is_supported_archive, unpack_zip_or_tar

# Find all modules for specific converters
l_converter_modules = glob.glob(os.path.dirname(base.__file__) + "/*.py")

try:

    class NameAndClass(NamedTuple):
        name: str
        converter_class: type[base.FileConverter]

    def get_converter_name_and_class(module_path: str) -> NameAndClass | None:

        module_name = os.path.splitext(os.path.basename(module_path))[0]

        # Skip the base module and the package __init__
        if module_name in ("base", "__init__"):
            return None

        package_name = "psdi_data_conversion.converters"
        module = importlib.import_module(f".{module_name}", package=package_name)

        # Check that the module defines a converter
        if not hasattr(module, "converter") or not issubclass(module.converter, base.FileConverter):
            print(f"ERROR: Module `{module_name}` in package `{package_name}` fails to define a converter to the "
                  "variable `converter` which is a subclass of `FileConverter`.", file=sys.stderr)
            return None

        converter_class = module.converter
        name = converter_class.name

        return NameAndClass(name, converter_class)

    # Get a list of all converter names and FileConverter subclasses
    l_converter_names_and_classes = [get_converter_name_and_class(module_name) for
                                     module_name in l_converter_modules]
    # Remove the None entry from the list, which corresponds to the 'base' module
    l_converter_names_and_classes = [x for x in l_converter_names_and_classes if x is not None]

    # Make constant dict and list of registered converters
    D_REGISTERED_CONVERTERS: dict[str, type[base.FileConverter]] = dict(l_converter_names_and_classes)
    L_REGISTERED_CONVERTERS: list[str] = [name for name in D_REGISTERED_CONVERTERS.keys()]

except Exception:
    print(f"ERROR: Failed to register converters. Exception was: \n{traceback.format_exc()}", file=sys.stderr)
    D_REGISTERED_CONVERTERS: dict[str, type[base.FileConverter]] = {}
    L_REGISTERED_CONVERTERS: list[str] = []


def get_converter(*args, name=const.CONVERTER_DEFAULT, **converter_kwargs) -> base.FileConverter:
    """Get a FileConverter of the proper subclass for the requested converter type

    Parameters
    ----------
    filename : str
        The filename of the input file to be converted, either relative to current directory or fully-qualified
    to_format : str
        The desired format to convert to, as the file extension (e.g. "cif")
    from_format : str | None
        The format to convert from, as the file extension (e.g. "pdb"). If None is provided (default), will be
        determined from the extension of `filename`
    name : str
        The desired converter type, by default 'Open Babel'
    data : dict[str | Any] | None
        A dict of any other data needed by a converter or for extra logging information, default empty dict
    abort_callback : Callable[[int], None]
        Function to be called if the conversion hits an error and must be aborted, default `abort_raise`, which
        raises an appropriate exception
    use_envvars : bool
        If set to True, environment variables will be checked for any that set options for this class and used,
        default False
    upload_dir : str
        The location of input files relative to the current directory
    download_dir : str
        The location of output files relative to the current directory
    max_file_size : float
        The maximum allowed file size for input/output files, in MB, default 1 MB. If 0, will be unlimited
    log_file : str | None
        If provided, all logging will go to a single file or stream. Otherwise, logs will be split up among multiple
        files for server-style logging.
    log_mode : str
        How logs should be stores. Allowed values are:
        - 'full' - Multi-file logging, only recommended when running as a public web app
        - 'simple' - Logs saved to one file
        - 'stdout' - Output logs and errors only to stdout
        - 'none' - Output only errors to stdout
    log_level : int | None
        The level to log output at. If None (default), the level will depend on the chosen `log_mode`:
        - 'full' or 'simple': INFO
        - 'stdout' - INFO to stdout, no logging to file
        - 'none' - ERROR to stdout, no logging to file
    refresh_local_log : bool
        If True, the local log generated from this run will be overwritten. If False it will be appended to. Default
        True
    delete_input : bool
        Whether or not to delete input files after conversion, default False

    Returns
    -------
    FileConverter
        A subclassed FileConverter for the desired converter type

    Raises
    ------
    FileConverterInputException
        If the converter isn't recognized or there's some other issue with the input
    """
    if name not in L_REGISTERED_CONVERTERS:
        raise base.FileConverterInputException(f"Converter {name} not recognized. Allowed converters are: " +
                                               f"{L_REGISTERED_CONVERTERS}")
    converter_class = D_REGISTERED_CONVERTERS[name]

    return converter_class(*args, **converter_kwargs)


def run_converter(filename,
                  *args,
                  from_format: str | None = None,
                  archive_output=True,
                  **converter_kwargs) -> str:
    """Shortcut to create and run a FileConverter in one step

    Parameters
    ----------
    filename : str
        Either the filename of the input file to be converted or of an archive file containing files to be converted
        (zip and tar supported), either relative to current directory or fully-qualified
    to_format : str
        The desired format to convert to, as the file extension (e.g. "cif")
    from_format : str | None
        The format to convert from, as the file extension (e.g. "pdb"). If None is provided (default), will be
        determined from the extension of `filename` if it's a simple file, or the contained files if `filename` is an
        archive file
    name : str
        The desired converter type, by default 'Open Babel'
    data : dict[str | Any] | None
        A dict of any other data needed by a converter or for extra logging information, default empty dict
    abort_callback : Callable[[int], None]
        Function to be called if the conversion hits an error and must be aborted, default `abort_raise`, which
        raises an appropriate exception
    use_envvars : bool
        If set to True, environment variables will be checked for any that set options for this class and used,
        default False
    upload_dir : str
        The location of input files relative to the current directory
    download_dir : str
        The location of output files relative to the current directory
    archive_output : bool
        If True (default) and the input file is an archive (i.e. zip or tar file), the converted files will be archived
        into a file of the same format, their logs will be combined into a single log, and the converted files and
        individual logs will be deleted
    max_file_size : float
        The maximum allowed file size for input/output files, in MB, default 1 MB. If 0, will be unlimited
    log_file : str | None
        If provided, all logging will go to a single file or stream. Otherwise, logs will be split up among multiple
        files for server-style logging.
    log_mode : str
        How logs should be stores. Allowed values are:
        - 'full' - Multi-file logging, only recommended when running as a public web app
        - 'simple' - Logs saved to one file
        - 'stdout' - Output logs and errors only to stdout
        - 'none' - Output only errors to stdout
    log_level : int | None
        The level to log output at. If None (default), the level will depend on the chosen `log_mode`:
        - 'full' or 'simple': INFO
        - 'stdout' - INFO to stdout, no logging to file
        - 'none' - ERROR to stdout, no logging to file
    refresh_local_log : bool
        If True, the local log generated from this run will be overwritten. If False it will be appended to. Default
        True
    delete_input : bool
        Whether or not to delete input files after conversion, default False

    Returns
    -------
    str
        A message briefly stating the conversion being performed

    Raises
    ------
    FileConverterInputException
        If the converter isn't recognized or there's some other issue with the input
    FileConverterAbortException
        If something goes wrong during the conversion process
    """

    # Check if the filename is for an archive file, and handle appropriately

    if not is_archive(filename):
        # Not an archive, so just get and run the converter straightforwardly
        return get_converter(*args, filename=filename, from_format=from_format, **converter_kwargs).run()

    if not is_supported_archive(filename):
        raise base.FileConverterInputException(f"{filename} is an unsupported archive type. Supported types are: "
                                               f"{const.L_SUPPORTED_ARCHIVE_EXTENSIONS}")

    # If we get here, the filename is of a supported archive type. Make a temporary directory to extract its contents
    # to, then run the converter on each file extracted
    out_str = ""
    with TemporaryDirectory() as extract_dir:
        l_filenames = unpack_zip_or_tar(filename, extract_dir=extract_dir)
        for extracted_filename in l_filenames:
            out_str += get_converter(*args, filename=extracted_filename,
                                     from_format=from_format, **converter_kwargs).run()
    return out_str
