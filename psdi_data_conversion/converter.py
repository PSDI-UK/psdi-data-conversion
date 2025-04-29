"""@file psdi-data-conversion/psdi_data_conversion/converter.py

Created 2024-12-10 by Bryan Gillis.

Class and functions to perform file conversion
"""

from dataclasses import dataclass, field
import json
import glob
import importlib
import os
import sys
import traceback
from typing import Any, Callable, NamedTuple
from multiprocessing import Lock
from psdi_data_conversion import log_utility
from collections.abc import Callable
from dataclasses import dataclass, field
from tempfile import TemporaryDirectory
from typing import Any, NamedTuple

from psdi_data_conversion import constants as const
from psdi_data_conversion.converters import base
from psdi_data_conversion.converters.openbabel import CONVERTER_OB
from psdi_data_conversion.file_io import (is_archive, is_supported_archive, pack_zip_or_tar, split_archive_ext,
                                          unpack_zip_or_tar)
from psdi_data_conversion.utils import regularize_name

# A lock to prevent multiple threads from logging at the same time
logLock = Lock()

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

        # To make querying case/space-insensitive, we store all names in lowercase with spaces stripped
        name = converter_class.name.lower().replace(" ", "")

        return NameAndClass(name, converter_class)

    # Get a list of all converter names and FileConverter subclasses
    l_converter_names_and_classes = [get_converter_name_and_class(module_name) for
                                     module_name in l_converter_modules]
    # Remove the None entry from the list, which corresponds to the 'base' module
    l_converter_names_and_classes = [x for x in l_converter_names_and_classes if x is not None]

    # Make constant dict and list of supported converters
    D_SUPPORTED_CONVERTERS: dict[str, type[base.FileConverter]] = dict(l_converter_names_and_classes)
    L_SUPPORTED_CONVERTERS: list[str] = [name for name in D_SUPPORTED_CONVERTERS.keys()]

    # Make constant dict and list of registered converters
    D_REGISTERED_CONVERTERS = {converter_name: converter_class for converter_name, converter_class in
                               D_SUPPORTED_CONVERTERS.items() if converter_class.can_be_registered()}
    L_REGISTERED_CONVERTERS: list[str] = [name for name in D_REGISTERED_CONVERTERS.keys()]

    # Make dicts of flags, options, and args (combined flags and options) for each converter
    _d_converter_flags, _d_converter_options, _d_converter_args = {}, {}, {}
    for name, converter_class in D_SUPPORTED_CONVERTERS.items():
        l_flags = converter_class.allowed_flags if converter_class.allowed_flags else ()
        l_options = converter_class.allowed_options if converter_class.allowed_options else ()
        _d_converter_flags[name] = l_flags
        _d_converter_options[name] = l_options
        _d_converter_args[name] = (*l_flags, *l_options)
    D_CONVERTER_FLAGS: dict[str, tuple[tuple[str, dict[str, Any], Callable]]] = _d_converter_flags
    D_CONVERTER_OPTIONS: dict[str, tuple[tuple[str, dict[str, Any], Callable]]] = _d_converter_options
    D_CONVERTER_ARGS: dict[str, tuple[tuple[str, dict[str, Any], Callable]]] = _d_converter_args

except Exception:
    print(f"ERROR: Failed to register converters. Exception was: \n{traceback.format_exc()}", file=sys.stderr)
    D_SUPPORTED_CONVERTERS: dict[str, type[base.FileConverter]] = {}
    L_SUPPORTED_CONVERTERS: list[str] = []
    D_REGISTERED_CONVERTERS: dict[str, type[base.FileConverter]] = {}
    L_REGISTERED_CONVERTERS: list[str] = []
    D_CONVERTER_FLAGS = {}
    D_CONVERTER_OPTIONS = {}
    D_CONVERTER_ARGS = {}


def get_supported_converter_class(name: str):
    """Get the appropriate converter class matching the provided name from the dict of supported converters

    Parameters
    ----------
    name : str
        Converter name (case- and space-insensitive)

    Returns
    -------
    type[base.FileConverter]
    """
    return D_SUPPORTED_CONVERTERS[regularize_name(name)]


def get_registered_converter_class(name: str):
    """Get the appropriate converter class matching the provided name from the dict of supported converters

    Parameters
    ----------
    name : str
        Converter name (case- and space-insensitive)

    Returns
    -------
    type[base.FileConverter]
    """
    return D_REGISTERED_CONVERTERS[regularize_name(name)]


def converter_is_supported(name: str):
    """Checks if a converter is supported in principle by this project

    Parameters
    ----------
    name : str
        Converter name (case- and space-insensitive)

    Returns
    -------
    bool
    """
    return regularize_name(name) in L_SUPPORTED_CONVERTERS


def converter_is_registered(name: str):
    """Checks if a converter is registered (usable)

    Parameters
    ----------
    name : str
        Converter name (case- and space-insensitive)

    Returns
    -------
    bool
    """
    return regularize_name(name) in L_REGISTERED_CONVERTERS


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
        A dict of any other data needed by a converter or for extra logging information, default empty dict. See the
        docstring of each converter for supported keys and values that can be passed to `data` here
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
        The maximum allowed file size for input/output files, in MB, default 1 MB for Open Babel, unlimited for other
        converters. If 0, will be unlimited. If an archive of files is provided, this will apply to the total of all
        files contained in it
    no_check : bool
        If False (default), will check at setup whether or not a conversion between the desired file formats is
        supported with the specified converter
    log_file : str | None
        If provided, all logging will go to a single file or stream. Otherwise, logs will be split up among multiple
        files for server-style logging.
    log_mode : str
        How logs should be stored. Allowed values are:
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
    name = regularize_name(name)
    if name not in L_REGISTERED_CONVERTERS:
        raise base.FileConverterInputException(const.ERR_CONVERTER_NOT_RECOGNISED.format(name) +
                                               f"{L_REGISTERED_CONVERTERS}")
    converter_class = get_registered_converter_class(name)

    return converter_class(*args, **converter_kwargs)


@dataclass
class FileConversionRunResult:
    """An object of this class will be output by the `run_converter` function on success to provide key info on
    the files created
    """
    # Lists of results from each individual conversion
    l_output_filenames: list[str] = field(default_factory=list)
    l_log_filenames: list[str] = field(default_factory=list)
    l_in_size: list[int] = field(default_factory=list)
    l_out_size: list[int] = field(default_factory=list)
    status_code: int = 0

    # If only one conversion was performed, these variables will hold the results for that conversion. Otherwise they
    # will point to summary files / hold the combined size
    output_filename: str | None = None
    log_filename: str | None = None
    in_size: int = field(init=False)
    out_size: int = field(init=False)

    def __post_init__(self):
        """Calculate appropriate values where possible - in_size and out_size are the sum of individual sizes, and if
        only one run was performed, we can set the output and log filenames to the filenames from that one run
        """
        if self.output_filename is None and len(self.l_output_filenames) == 1:
            self.output_filename = self.l_output_filenames[0]
        if self.log_filename is None and len(self.l_log_filenames) == 1:
            self.log_filename = self.l_log_filenames[0]

        self.in_size = sum(self.l_in_size)
        self.out_size = sum(self.l_out_size)


def check_from_format(filename: str,
                      from_format: str | int,
                      strict=False) -> bool:
    """Check that the filename for an input file ends with the expected extension

    Parameters
    ----------
    filename : str
        The filename
    from_format : str | int
        The expected format (extension)
    strict : bool, optional
        If True, will raise an exception on failure. Otherwise will print a warning and return False

    Returns
    -------
    bool
        Whether the file ends with the expected extension or not

    Raises
    ------
    base.FileConverterInputException
        If `strict` is True and the the file does not end with the expected exception
    """

    # Get the name of the format
    if isinstance(from_format, str):
        from_format_name = from_format
    else:
        from psdi_data_conversion.database import get_format_info
        from_format_name = get_format_info(from_format).name

    # Silently make sure `from_format` starts with a dot
    if not from_format_name.startswith("."):
        from_format_name = f".{from_format}"

    if filename.endswith(from_format_name):
        return True

    msg = const.ERR_WRONG_EXTENSIONS.format(file=os.path.basename(filename), ext=from_format_name)

    if strict:
        raise base.FileConverterInputException(msg)

    print(f"WARNING: {msg}", file=sys.stderr)

    return False


def run_converter(filename: str,
                  to_format: str,
                  *args,
                  from_format: str | None = None,
                  download_dir=const.DEFAULT_DOWNLOAD_DIR,
                  max_file_size=None,
                  log_file: str | None = None,
                  log_mode=const.LOG_SIMPLE,
                  strict=False,
                  archive_output=True,
                  **converter_kwargs) -> FileConversionRunResult:
    """Shortcut to create and run a FileConverter in one step

    Parameters
    ----------
    filename : str
        Either the filename of the input file to be converted or of an archive file containing files to be converted
        (zip and tar supported), either relative to current directory or fully-qualified. If an archive is provided,
        the contents will be converted and then packed into an archive of the same type
    to_format : str
        The desired format to convert to, as the file extension (e.g. "cif")
    from_format : str | None
        The format to convert from, as the file extension (e.g. "pdb"). If None is provided (default), will be
        determined from the extension of `filename` if it's a simple file, or the contained files if `filename` is an
        archive file
    name : str
        The desired converter type, by default 'Open Babel'
    data : dict[str | Any] | None
        A dict of any other data needed by a converter or for extra logging information, default empty dict. See the
        docstring of each converter for supported keys and values that can be passed to `data` here
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
    strict : bool
        If True and `from_format` is not None, will fail if any input file has the wrong extension (including files
        within archives, but not the archives themselves). Otherwise, will only print a warning in this case
    archive_output : bool
        If True (default) and the input file is an archive (i.e. zip or tar file), the converted files will be archived
        into a file of the same format, their logs will be combined into a single log, and the converted files and
        individual logs will be deleted
    max_file_size : float
        The maximum allowed file size for input/output files, in MB, default 1 MB for Open Babel, unlimited for other
        converters. If 0, will be unlimited. If an archive of files is provided, this will apply to the total of all
        files contained in it
    no_check : bool
        If False (default), will check at setup whether or not a conversion between the desired file formats is
        supported with the specified converter
    log_file : str | None
        If provided, all logging will go to a single file or stream. Otherwise, logs will be split up among multiple
        files for server-style logging.
    log_mode : str
        How logs should be stored. Allowed values are:
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
    FileConversionRunResult
        An object containing the filenames of output files and logs created, and input/output file sizes

    Raises
    ------
    FileConverterInputException
        If the converter isn't recognized or there's some other issue with the input
    FileConverterAbortException
        If something goes wrong during the conversion process
    """

    # Set the maximum file size based on which converter is being used if using the default
    if max_file_size is None:
        max_file_size = const.DEFAULT_MAX_FILE_SIZE_OB if name == CONVERTER_OB else const.DEFAULT_MAX_FILE_SIZE

    # Set the log file if it was unset - note that in server logging mode, this value won't be used within the
    # converter class, so it needs to be set up here to match what will be set up there
    if log_file is None:
        base_filename = os.path.basename(split_archive_ext(filename)[0])
        log_file = os.path.join(download_dir, base_filename + const.OUTPUT_LOG_EXT)

    # Check if the filename is for an archive file, and handle appropriately

    l_run_output: list[base.FileConversionResult] = []

    file_is_archive = is_archive(filename)

    # Status code for the overall success of the process
    status_code = 0

    if not file_is_archive:
        # Not an archive, so just get and run the converter straightforwardly
        if from_format is not None:
            check_from_format(filename, from_format, strict=strict)
        l_run_output.append(get_converter(filename,
                            to_format,
                            *args,
                            from_format=from_format,
                            download_dir=download_dir,
                            max_file_size=max_file_size,
                            log_file=log_file,
                            log_mode=log_mode,
                            **converter_kwargs).run())

    elif not is_supported_archive(filename):
        raise base.FileConverterInputException(f"{filename} is an unsupported archive type. Supported types are: "
                                               f"{const.D_SUPPORTED_ARCHIVE_FORMATS}")

    else:
        # The filename is of a supported archive type. Make a temporary directory to extract its contents
        # to, then run the converter on each file extracted
        with TemporaryDirectory() as extract_dir:
            l_filenames = unpack_zip_or_tar(filename, extract_dir=extract_dir)

            # Check for no files in archive
            if len(l_filenames) == 0:
                raise base.FileConverterInputException(const.ERR_EMPTY_ARCHIVE)

            # First check for files of invalid type, to avoid converting if one will cause a failure
            if from_format is not None:
                for extracted_filename in l_filenames:
                    check_from_format(extracted_filename, from_format, strict=strict)

            # Keep track of the file size budget
            remaining_file_size = max_file_size

            for extracted_filename in l_filenames:
                # Make a filename for the log for this particular conversion
                individual_log_file = os.path.join(extract_dir,
                                                   os.path.basename(os.path.splitext(extracted_filename)[0]) +
                                                   const.OUTPUT_LOG_EXT)

                # If the log mode is "full", set it to "full-force" for the individual runs to force use of the log file
                # name we set up for it
                individual_log_mode = log_mode if log_mode != const.LOG_FULL else const.LOG_FULL_FORCE

                try:
                    individual_run_output = get_converter(extracted_filename,
                                                          to_format,
                                                          *args,
                                                          from_format=from_format,
                                                          download_dir=download_dir,
                                                          log_file=individual_log_file,
                                                          log_mode=individual_log_mode,
                                                          max_file_size=remaining_file_size,
                                                          **converter_kwargs).run()
                except base.FileConverterAbortException as e:
                    # If the run fails, create a run output object to indicate that
                    individual_run_output = base.FileConversionResult(log_filename=individual_log_file,
                                                                      status_code=e.status_code)
                    status_code = max((status_code, e.status_code))
                    # If we specifically have a failure due to the size being exceeded, stop here, since no further
                    # runs are allowed
                    if isinstance(e, base.FileConverterSizeException):
                        l_run_output.append(individual_run_output)
                        break

                l_run_output.append(individual_run_output)

                # Reduce the file size limit by how much was used here
                remaining_file_size -= max((individual_run_output.in_size, individual_run_output.out_size))

            # Combine the output logs into a single log
            with open(log_file, "w") as fo:
                for individual_run_output in l_run_output:
                    individual_log_filename = individual_run_output.log_filename
                    if not os.path.exists(individual_log_filename):
                        raise base.FileConverterException(f"Expected log file '{individual_log_filename}' cannot be "
                                                          "found")
                    fo.write(open(individual_log_filename, "r").read() + "\n")
                    os.remove(individual_log_filename)

    # Combine the possibly-multiple FileConversionResults objects into a single FileConversionRunResult
    run_output = FileConversionRunResult(*zip(*[(x.output_filename,
                                                 x.log_filename,
                                                 x.in_size,
                                                 x.out_size) for x in l_run_output]),
                                         log_filename=log_file,
                                         status_code=status_code)

    if file_is_archive and archive_output:
        # If we get here, the file is an archive and we want to archive the output

        # Prune any unsuccessful runs from the list of output files
        l_successful_files = [x for x in run_output.l_output_filenames if x is not None]

        if len(l_successful_files) > 0:

            # Determine the directory for the output from the output filenames
            downloads_dir = os.path.split(l_successful_files[0])[0]

            # Create new names for the archive file and log file
            filename_base, ext = split_archive_ext(os.path.basename(filename))
            run_output.output_filename = os.path.join(downloads_dir, f"{filename_base}-{to_format}{ext}")

            # Pack the output files into an archive, cleaning them up afterwards
            pack_zip_or_tar(run_output.output_filename,
                            l_successful_files,
                            cleanup=True)

        # If the run was ultimately unsuccessful, raise an exception now, referencing the output log and including
        # error lines in it
        if status_code:
            msg = const.ERR_CONVERSION_FAILED.format(run_output.log_filename)
            l_output_log_lines = open(run_output.log_filename, "r").read().splitlines()
            l_error_lines = [line for line in l_output_log_lines if "ERROR" in line]
            msg += "\n".join(l_error_lines)
            if status_code == const.STATUS_CODE_SIZE:
                exception_class = base.FileConverterSizeException
            else:
                exception_class = base.FileConverterAbortException
            raise exception_class(status_code, msg)

    # Log conversion information if in service mode
    service_mode_ev = os.environ.get(const.SERVICE_MODE_EV)
    service_mode = (service_mode_ev is not None) and (service_mode_ev.lower() == "true")
    if service_mode:
        try:
            l_index = filename.rfind('/') + 1
            r_index = len(filename)
            in_filename = filename[l_index:r_index]

            l_index = run_output.output_filename.rfind('/') + 1
            r_index = len(run_output.output_filename)

            input_size = set_size_units(run_output.in_size)
            output_size = set_size_units(run_output.out_size)

            if status_code:
                outcome = "failed"
                fail_reason = l_error_lines
            else:
                outcome = "succeeded"
                fail_reason = ""

            entry = {
                "datetime": log_utility.get_date_time(),
                "input_format": converter_kwargs['data']['from_full'],
                "output_format": converter_kwargs['data']['to_full'],
                "input_filename": in_filename,
                "output_filename": run_output.output_filename[l_index:r_index],
                "input_size": input_size,
                "output_size": output_size }

            for key in [ "converter", "coordinates", "coordOption", "from_flags",
                "to_flags", "from_arg_flags", "to_arg_flags" ]:
                if key in converter_kwargs['data'] and converter_kwargs['data'][key] != "" and not \
                    ((key == "coordinates" or key == "coordOption") and converter_kwargs['data']['coordinates'] == "neither") :
                    entry[key] = converter_kwargs['data'][key]

            entry["outcome"] = outcome

            if fail_reason != "":
                entry["fail_reason"] = fail_reason

            logLock.acquire()
            sys.__stdout__.write(f"{json.dumps(entry) + '\n'}")
            logLock.release()
        except Exception:
            sys.__stdout__.write({"datetime": log_utility.get_date_time(),
                                  "logging_error": "An error occurred during logging of conversion information."})

    return run_output

def set_size_units(size):
    if size >= 1024:
        return str('%.3f' % (size / 1024)) + ' kB'
    elif size >= const.MEGABYTE:
        return str(size / const.MEGABYTE) + ' MB'
    else:
        return str(size) + ' B'
