"""@file psdi_data_conversion/converters/base.py

Created 2025-01-23 by Bryan Gillis.

Base class and information for file format converters
"""


import abc
import logging
import os
import subprocess
import sys
import traceback
from collections.abc import Callable
from copy import deepcopy
from dataclasses import dataclass
from typing import Any

from psdi_data_conversion import constants as const
from psdi_data_conversion import log_utility
from psdi_data_conversion.dist import bin_exists, get_bin_path, get_dist
from psdi_data_conversion.security import SAFE_STRING_RE, string_is_safe

try:
    # werkzeug is installed in the optional dependency Flask. It's only used here to recognize an exception type,
    # and if Flask isn't installed, that exception will never be raised, so we can just replace it with None and later
    # not try to catch it if werkzeug isn't found
    from werkzeug.exceptions import HTTPException
except ImportError:
    HTTPException = None


class FileConverterException(RuntimeError):
    """Exception class to represent any runtime error encountered by this package.
    """

    def __init__(self,
                 *args,
                 logged: bool = False,
                 help: bool = False,
                 msg_preformatted: bool = False):
        super().__init__(*args)
        self.logged = logged
        self.help = help
        self.msg_preformatted = msg_preformatted


class FileConverterAbortException(FileConverterException):
    """Class representing an exception triggered by a call to abort a file conversion
    """

    def __init__(self,
                 status_code: int,
                 *args,
                 **kwargs):
        super().__init__(*args, **kwargs)
        self.status_code = status_code


class FileConverterSizeException(FileConverterAbortException):
    """Class representing an exception triggered by the maximum size being exceeded
    """

    def __init__(self,
                 *args,
                 in_size: int | None = None,
                 out_size: int | None = None,
                 max_file_size: int | None = None,
                 **kwargs):
        super().__init__(*args, **kwargs)
        self.in_size = in_size
        self.out_size = out_size
        self.max_file_size = max_file_size


class FileConverterInputException(FileConverterException):
    """Exception class to represent errors encountered with input parameters for the data conversion script.
    """
    pass


if HTTPException is not None:
    l_abort_exceptions = (HTTPException, FileConverterAbortException)
else:
    l_abort_exceptions = (FileConverterAbortException,)


@dataclass
class FileConversionResult:
    """An object of this class will be output by the file converter's `run` function on success to provide key info on
    the files created
    """
    output_filename: str | None = None
    log_filename: str | None = None
    in_size: int = 0
    out_size: int = 0
    status_code: int = 0


def abort_raise(status_code: int,
                *args,
                e: Exception | None = None,
                **kwargs):
    """Callback for aborting during a file conversion, which passes relevant information to an exception of the
    appropriate type
    """
    if e:
        raise e
    elif status_code == const.STATUS_CODE_SIZE:
        exception_class = FileConverterSizeException
    else:
        exception_class = FileConverterAbortException
    raise exception_class(status_code, *args, **kwargs)


class FileConverter:
    """Class to handle conversion of files from one type to another
    """

    # Class variables and methods which must/can be overridden by subclasses
    # ----------------------------------------------------------------------

    name: str | None = None
    """Name of the converter - must be overridden in each subclass to name each converter uniquely"""

    info: str | None = None
    """General info about the converter - can be overridden in a subclass to add information about a converter which
    isn't covered in its database entry, such as notes on its support."""

    allowed_flags: tuple[tuple[str, dict, Callable], ...] | None = None
    """List of flags allowed for the converter (flags are arguments that are set by being present, and don't require a
    value specified - e.g. "-v" to enable verbose mode) - should be overridden with a tuple of tuples containing the
    flag names, a dict of kwargs to be passed to the argument parser's `add_argument` method, and callable function to
    get a dict of needed info for them. If the converter does not accept any flags, an empty tuple should be supplied
    (e.g `allowed_flags = ()`), as `None` will be interpreted as this value not having been overridden"""

    allowed_options: tuple[tuple[str, dict, Callable], ...] | None = None
    """List of options allowed for the converter (options are arguments that take one or more values, e.g. "-o out.txt")
    - should be overridden with a tuple of tuples containing the option names, a dict of kwargs to be passed to the
    argument parser's `add_argument` method, and callable function to get a dict of needed info for them.
    As with flags, an empty tuple should be provided if the converter does not accept any options"""

    database_key_prefix: str | None = None
    """The prefix used in the database for keys related to this converter"""

    supports_ambiguous_extensions: bool = False
    """Whether or not this converter supports formats which share the same extension. This is used to enforce stricter
    but less user-friendly requirements on format specification"""

    @abc.abstractmethod
    def _convert(self):
        """Run the conversion with the desired converter. This must be implemented for each converter class.
        """
        pass

    @classmethod
    def can_be_registered(cls) -> bool:
        """If the converter class may not be able to be registered (for instance, it relies on a binary which isn't
        supported on all platforms), this method should be overridden to perform necessary checks to indicate if it
        can be registered or not.
        """
        return True

    # If the converter supports flags specific to the input file format, set the below to True for the subclass so help
    # text will be properly displayed notifying the user that they can request this by providing an input format (and
    # similar for the other similar class variables below)
    has_in_format_flags_or_options = False
    has_out_format_flags_or_options = False

    @staticmethod
    def get_in_format_flags(in_format: str) -> tuple[tuple[str, str], ...]:
        """Gets flags which are applicable for a specific input file format, returned as a tuple of (flag, description).
        This should be overridden for each converter class if it uses any format-specific input flags.
        """
        return ()

    @staticmethod
    def get_out_format_flags(in_format: str) -> tuple[tuple[str, str], ...]:
        """Gets flags which are applicable for a specific output file format, returned as a tuple of (flag,
        description). This should be overridden for each converter class if it uses any format-specific output flags.
        """
        return ()

    @staticmethod
    def get_in_format_options(in_format: str) -> tuple[tuple[str, str], ...]:
        """Gets options which are applicable for a specific input file format, returned as a tuple of (option,
        description). This should be overridden for each converter class if it uses any format-specific input options.
        """
        return ()

    @staticmethod
    def get_out_format_options(in_format: str) -> tuple[tuple[str, str], ...]:
        """Gets options which are applicable for a specific output file format, returned as a tuple of (option,
        description). This should be overridden for each converter class if it uses any format-specific output options.
        """
        return ()

    # Base class functionality
    # ------------------------

    def __init__(self,
                 filename: str,
                 to_format: str,
                 from_format: str | None = None,
                 data: dict[str, Any] | None = None,
                 abort_callback: Callable[[int], None] = abort_raise,
                 use_envvars=False,
                 upload_dir=const.DEFAULT_UPLOAD_DIR,
                 download_dir=const.DEFAULT_DOWNLOAD_DIR,
                 max_file_size=None,
                 no_check=False,
                 log_file: str | None = None,
                 log_mode=const.LOG_FULL,
                 log_level: int | None = None,
                 refresh_local_log: bool = True,
                 delete_input=False):
        """Initialize the object, storing needed data and setting up loggers.

        Parameters
        ----------
        filename : str
            The filename of the input file to be converted, either relative to current directory or fully-qualified
        to_format : str
            The desired format to convert to, as the file extension (e.g. "cif")
        from_format : str | None
            The format to convert from, as the file extension (e.g. "pdb"). If None is provided (default), will be
            determined from the extension of `filename`
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
            The maximum allowed file size for input/output files, in MB. If 0, will be unlimited. Default 0 (unlimited)
        no_check : bool
            If False (default), will check at setup whether or not a conversion between the desired file formats is
            supported with the specified converter
        log_file : str | None
            If provided, all logging will go to a single file or stream. Otherwise, logs will be split up among multiple
            files for server-style logging.
        log_mode : str
            How logs should be stores. Allowed values are:
            - 'full' - Multi-file logging, only recommended when running as a public web app
            - 'full-force' - Multi-file logging, only recommended when running as a public web app, with the log file
                name forced to be used for the output log
            - 'simple' - Logs saved to one file
            - 'stdout' - Output logs and errors only to stdout
            - 'none' - Output only errors to stdout
        log_level : int | None
            The level to log output at. If None (default), the level will depend on the chosen `log_mode`:
            - 'full', 'full-force', or 'simple': INFO
            - 'stdout' - INFO to stdout, no logging to file
            - 'none' - ERROR to stdout, no logging to file
        refresh_local_log : bool
            If True, the local log generated from this run will be overwritten. If False it will be appended to. Default
            True
        delete_input : bool
            Whether or not to delete input files after conversion, default False
        """

        # Wrap the initialisation in a try block, calling the abort_callback function if anything goes wrong
        self.abort_callback = abort_callback

        try:

            if max_file_size is None:
                from psdi_data_conversion.converters.openbabel import CONVERTER_OB
                if self.name == CONVERTER_OB:
                    self.max_file_size = const.DEFAULT_MAX_FILE_SIZE_OB
                else:
                    self.max_file_size = const.DEFAULT_MAX_FILE_SIZE
            else:
                self.max_file_size = max_file_size*const.MEGABYTE

            # Set values from envvars if desired
            if use_envvars:
                # Get the maximum allowed size from the envvar for it
                from psdi_data_conversion.converters.openbabel import CONVERTER_OB
                if self.name == CONVERTER_OB:
                    ev_max_file_size = os.environ.get(const.MAX_FILESIZE_OB_EV)
                else:
                    ev_max_file_size = os.environ.get(const.MAX_FILESIZE_EV)

                if ev_max_file_size is not None:
                    self.max_file_size = float(ev_max_file_size)*const.MEGABYTE

            # Set member variables directly from input
            self.in_filename = filename
            self.to_format = to_format
            self.upload_dir = upload_dir
            self.download_dir = download_dir
            self.log_file = log_file
            self.log_mode = log_mode
            self.log_level = log_level
            self.refresh_local_log = refresh_local_log
            self.delete_input = delete_input

            # Use an empty dict for data if None was provided
            if data is None:
                self.data = {}
            else:
                self.data = dict(deepcopy(data))

            # Get from_format from the input file extension if not supplied
            if from_format is None:
                self.from_format = os.path.splitext(self.in_filename)[1]
            else:
                self.from_format = from_format

            # Convert in and out formats to FormatInfo, and raise an exception if one is ambiguous
            from psdi_data_conversion.database import disambiguate_formats
            (self.from_format_info,
             self.to_format_info) = disambiguate_formats(self.name, self.from_format, self.to_format)

            # Set placeholders for member variables which will be set when conversion is run
            self.in_size: int | None = None
            self.out_size: int | None = None
            self.out: str | None = None
            self.err: str | None = None
            self.quality: str | None = None

            # Create directory 'uploads' if not extant.
            if not os.path.exists(self.upload_dir):
                os.makedirs(self.upload_dir, exist_ok=True)

            # Create directory 'downloads' if not extant.
            if not os.path.exists(self.download_dir):
                os.makedirs(self.download_dir, exist_ok=True)

            self.local_filename = os.path.split(self.in_filename)[1]
            self.filename_base = os.path.splitext(self.local_filename)[0]
            self.out_filename = f"{self.download_dir}/{self.filename_base}.{self.to_format_info.name}"

            # Set up files to log to
            self._setup_loggers()

            # Check that the requested conversion is valid and warn of any issues unless suppressed
            if not no_check:
                from psdi_data_conversion.database import get_conversion_quality
                qual = get_conversion_quality(self.name,
                                              self.from_format_info.id,
                                              self.to_format_info.id)
                if not qual:
                    raise FileConverterInputException(f"Conversion from {self.from_format_info.name} to "
                                                      f"{self.to_format_info.name} "
                                                      f"with {self.name} is not supported.", help=True)
                if qual.details:
                    msg = (":\nPotential data loss or extrapolation issues with the conversion from "
                           f"{self.from_format_info.name} to {self.to_format_info.name}:\n")
                    for detail_line in qual.details.split("\n"):
                        msg += f"- {detail_line}\n"
                    self.logger.warning(msg)

            self.logger.debug("Finished FileConverter initialisation")

        except Exception as e:
            # Don't catch a deliberate abort; let it pass through
            if isinstance(e, l_abort_exceptions):
                if not hasattr(e, "logged") or e.logged is False:
                    self.logger.error(f"Unexpected exception raised while running the converter, of type '{type(e)}' "
                                      f"with message: {str(e)}")
                    if e:
                        e.logged = True
                raise
            # Try to run the standard abort method. There's a good chance this will fail though depending on what went
            # wrong when during init, so we fallback to printing the exception to stderr
            try:
                if not (isinstance(e, FileConverterException) and e.help):
                    self.logger.error(f"Exception triggering an abort was raised while initializing the converter. "
                                      f"Exception was type '{type(e)}', with message: {str(e)}")
                    if e:
                        e.logged = True
                self._abort(message="The application encountered an error while initializing the converter:\n" +
                            traceback.format_exc(), e=e)
            except Exception as ee:
                if isinstance(ee, l_abort_exceptions) or (isinstance(ee, FileConverterException) and ee.help):
                    # Don't catch a deliberate abort or help exception; let it pass through
                    raise
                message = ("ERROR: The application encounted an error during initialization of the converter and "
                           "could not cleanly log the error due to incomplete init: " + traceback.format_exc())
                print(message, file=sys.stderr)
                try:
                    self.abort_callback(const.STATUS_CODE_GENERAL, message, e=e)
                except TypeError:
                    self.abort_callback(const.STATUS_CODE_GENERAL)

    def _setup_loggers(self):
        """Run at init to set up loggers for this object.
        """

        # Determine level to log at based on quiet status
        if self.log_level:
            self._local_logger_level = self.log_level
            self._stdout_output_level = self.log_level
        else:
            if self.log_mode == const.LOG_NONE:
                self._local_logger_level = None
                self._stdout_output_level = logging.ERROR
            elif self.log_mode == const.LOG_STDOUT:
                self._local_logger_level = None
                self._stdout_output_level = logging.INFO
            elif self.log_mode in (const.LOG_FULL, const.LOG_FULL_FORCE, const.LOG_SIMPLE):
                self._local_logger_level = const.DEFAULT_LOCAL_LOGGER_LEVEL
                self._stdout_output_level = logging.ERROR
            else:
                raise FileConverterInputException(f"ERROR: Unrecognised logging option: {self.log_mode}. Allowed "
                                                  f"options are: {const.L_ALLOWED_LOG_MODES}")
        if self.log_mode in (const.LOG_FULL, const.LOG_FULL_FORCE):
            return self._setup_server_loggers()

        self.output_log = self.log_file

        write_mode = "w" if self.refresh_local_log else "a"
        self.logger = log_utility.set_up_data_conversion_logger(local_log_file=self.log_file,
                                                                local_logger_level=self._local_logger_level,
                                                                stdout_output_level=self._stdout_output_level,
                                                                suppress_global_handler=True,
                                                                mode=write_mode)

        self.logger.debug(f"Set up logging in log mode '{self.log_mode}'")
        if self.log_level:
            self.logger.debug(f"Logging level set to {self.log_level}")
        else:
            self.logger.debug(f"Logging level left to defaults. Using {self._local_logger_level} for local logger "
                              f"and {self._stdout_output_level} for stdout output")

    def _setup_server_loggers(self):
        """Run at init to set up loggers for this object in server-style execution
        """
        # For server mode, we need a specific log name, so set that up unless the mode is set to force the use of
        # the input log file
        if self.log_mode == const.LOG_FULL_FORCE:
            self.output_log = self.log_file
        else:
            self.output_log = os.path.join(self.download_dir, f"{self.filename_base}{const.OUTPUT_LOG_EXT}")

        # If any previous log exists, delete it
        if os.path.exists(self.output_log):
            os.remove(self.output_log)

        write_mode = "w" if self.refresh_local_log else "a"
        # Set up loggers - one for general-purpose log_utility, and one just for what we want to output to the user
        self.logger = log_utility.set_up_data_conversion_logger(local_log_file=self.output_log,
                                                                local_logger_level=self._local_logger_level,
                                                                stdout_output_level=self._stdout_output_level,
                                                                local_logger_raw_output=False,
                                                                mode=write_mode)

        self.logger.debug(f"Set up server-style logging, with user logging at level {self._local_logger_level}")

    def run(self):
        """Run the file conversion
        """

        try:
            self.logger.debug("Checking input file size")
            self._check_input_file_size_and_status()

            self.logger.debug("Starting file conversion")
            self._convert()

            self.logger.debug("Finished file conversion; performing cleanup tasks")
            self._finish_convert()
        except Exception as e:
            # Don't catch a deliberate abort; let it pass through
            if isinstance(e, l_abort_exceptions):
                # Log the error if it hasn't yet been logged
                if not hasattr(e, "logged") or e.logged is False:
                    self.logger.error(f"Unexpected exception raised while running the converter, of type '{type(e)}' "
                                      f"with message: {str(e)}")
                    e.logged = True
                raise
            if not (isinstance(e, FileConverterException) and e.help):
                self.logger.error(f"Exception triggering an abort was raised while running the converter. Exception "
                                  f"was type '{type(e)}', with message: {str(e)}")
                if e:
                    e.logged = True
            self._abort(message="The application encountered an error while running the converter:\n" +
                        traceback.format_exc(), e=e)

        return FileConversionResult(output_filename=self.out_filename,
                                    log_filename=self.output_log,
                                    in_size=self.in_size,
                                    out_size=self.out_size)

    def _abort(self,
               status_code: int = const.STATUS_CODE_GENERAL,
               message: str | None = None,
               e: Exception | None = None,
               **kwargs):
        """Abort the conversion, reporting the desired message to the user at the top of the output

        Parameters
        ----------
        status_code : int
            The HTTP status code to exit with. Default is 422: Unprocessable Content
        message : str | None
            If provided, this message will be logged in the user output log at the top of the file and will appear in
            any raised exception if possible. This should typically explain the reason the process failed
        e : Exception | None
            The caught exception which triggered this abort, if any
        **kwargs : Any
            Any additional keyword arguments are passed to the `self.abort_callback` function if it accepts them

        """

        def try_debug_log(msg, *args, **kwargs):
            try:
                self.logger.debug(msg, *args, **kwargs)
            except AttributeError:
                pass

        def error_log(msg, *args, **kwargs):
            try:
                self.logger.error(msg, *args, **kwargs)
            except AttributeError:
                print(msg, file=sys.stderr)

        # Remove the input and output files if they exist
        if self.delete_input:
            self.logger.debug(f"Cleaning up input file {self.in_filename}")
            try:
                os.remove(self.in_filename)
            except FileNotFoundError:
                pass

        try:
            os.remove(self.out_filename)
        except (FileNotFoundError, AttributeError):
            try_debug_log("Application aborting; no output file found to clean up")
        else:
            try_debug_log(f"Application aborting, so cleaning up output file {self.out_filename}")

        # If we have a Help exception, override the message with its message
        if isinstance(e, FileConverterException) and e.help:
            try_debug_log("Help exception triggered, so only using its message for output")
            message = str(e)

        if message:
            # If we're adding a message in server mode, read in any prior logs, clear the log, write the message, then
            # write the prior logs
            if self.log_file is None:
                try_debug_log("Adding abort message to the top of the output log so it will be the first thing "
                              "read by the user")
                prior_output_log = open(self.output_log, "r").read()
                os.remove(self.output_log)
                with open(self.output_log, "w") as fo:
                    fo.write(message + "\n")
                    fo.write(prior_output_log)

            # Note this message in the error logger as well
            if not (isinstance(e, FileConverterException) and e.help):
                error_log(message)
                if e:
                    e.logged = True

        # Call the abort callback function now. We first try passing information to the callback function
        try:
            self.abort_callback(status_code, message, e=e, **kwargs)
        except TypeError:
            # The callback function doesn't support arguments, so we instead call the callback, catch any exception it
            # raises, monkey-patch on the extra info, and reraise it
            try:
                self.abort_callback(status_code)
            except Exception as ee:
                ee.status_code = status_code
                ee.message = message
                ee.e = e
                ee.abort_kwargs = kwargs
                raise ee

    def _abort_from_err(self):
        """Call an abort after a call to the converter has completed, but it's returned an error. Create a message for
        the logger including this error and other relevant information.
        """
        self.logger.error(self._create_message_start() +
                          self._create_message() +
                          self.out + '\n' +
                          self.err)
        self._abort(message=self.err, logged=True)

    def _create_message(self) -> str:
        """Create a log of options passed to the converter - this method should be overloaded to log any information
        unique to a specific converter.
        """

        self.logger.debug("Default _create_message method called - not outputting any additional information specific "
                          "to this converter")

        return ""

    def _create_message_start(self) -> str:
        """Create beginning of message for log files

        Returns
        -------
        str
            The beginning of a message for log files, containing generic information about what was trying to be done
        """
        # We want the entries to all line up, so we need a dummy line at the top to force a newline break - anything
        # empty or whitespace will be stripped by the logger, so we use a lone colon, which looks least obtrusive
        return (":\n"
                f"File name:         {self.filename_base}\n"
                f"From:              {self.from_format_info.name} ({self.from_format_info.note})\n"
                f"To:                {self.to_format} ({self.to_format_info.note})\n"
                f"Converter:         {self.name}\n")

    def _log_success(self):
        """Write conversion information to server-side file, ready for downloading to user
        """

        message = (self._create_message_start()+self._create_message() +
                   'Quality:           ' + self.quality + '\n'
                   'Success:           Assuming that the data provided was of the correct format, the conversion\n'
                   '                   was successful (to the best of our knowledge) subject to any warnings below.\n' +
                   self.out + '\n' + self.err).strip() + '\n'

        self.logger.info(message)

    def _check_input_file_size_and_status(self):
        """Get input file size and status, checking that the file isn't too large
        """

        try:
            self.in_size = os.path.getsize(os.path.realpath(self.in_filename))
        except FileNotFoundError:
            # Something went wrong and the output file doesn't exist
            err_message = f"Expected output file {self.in_filename} does not exist."
            self.logger.error(err_message)
            self.err += f"ERROR: {err_message}\n"
            self._abort_from_err()

        # Check that the input file doesn't exceed the maximum allowed size
        if self.max_file_size > 0 and self.in_size > self.max_file_size:

            self._abort(const.STATUS_CODE_SIZE,
                        f"ERROR converting {os.path.basename(self.in_filename)} to " +
                        os.path.basename(self.out_filename) + ": "
                        f"Input file exceeds maximum size.\nInput file size is "
                        f"{self.in_size/const.MEGABYTE:.2f} MB; maximum input file size is "
                        f"{self.max_file_size/const.MEGABYTE:.2f} MB.\n",
                        max_file_size=self.max_file_size,
                        in_size=self.in_size,
                        out_size=None)

    def _check_output_file_size_and_status(self):
        """Get output file size and status, checking that the file isn't too large
        """

        try:
            self.out_size = os.path.getsize(os.path.realpath(self.out_filename))
        except FileNotFoundError:
            # Something went wrong and the output file doesn't exist
            err_message = f"Expected output file {self.out_filename} does not exist."
            self.logger.error(err_message)
            self.err += f"ERROR: {err_message}\n"
            self._abort_from_err()

        # Check that the output file doesn't exceed the maximum allowed size
        if self.max_file_size > 0 and self.out_size > self.max_file_size:

            self._abort(const.STATUS_CODE_SIZE,
                        f"ERROR converting {os.path.basename(self.in_filename)} to " +
                        os.path.basename(self.out_filename) + ": "
                        f"Output file exceeds maximum size.\nInput file size is "
                        f"{self.in_size/const.MEGABYTE:.2f} MB; Output file size is {self.out_size/const.MEGABYTE:.2f} "
                        f"MB; maximum output file size is {self.max_file_size/const.MEGABYTE:.2f} MB.\n",
                        max_file_size=self.max_file_size,
                        in_size=self.in_size,
                        out_size=self.out_size)

        self.logger.debug(f"Output file found to have size {self.out_size/const.MEGABYTE:.2f} MB")

    def get_quality(self) -> str:
        """Query the JSON file to obtain conversion quality
        """
        from psdi_data_conversion.database import get_conversion_quality

        conversion_quality = get_conversion_quality(converter_name=self.name,
                                                    in_format=self.from_format_info.id,
                                                    out_format=self.to_format_info.id)
        if not conversion_quality:
            return "unknown"
        return conversion_quality.qual_str

    def _finish_convert(self):
        """Run final common steps to clean up a conversion and log success or abort due to an error
        """

        self._check_output_file_size_and_status()

        if self.delete_input:
            os.remove(self.in_filename)
        if "success" in self.data:
            self.quality = self.data["success"]
        else:
            self.quality = self.get_quality()

        self._log_success()


class ScriptFileConverter(FileConverter):
    """File Converter specialized to run a shell script to call the converter
    """

    script: str | None = None
    """The name of the script to run this converter, relative to the ``psdi_data_conversion/scripts`` directory"""

    required_bin: str | None = None
    """The name of the binary called by the script, relative to the ``psdi_data_conversion/bin/$DIST`` directory,
    where `DIST` is 'linux', 'windows', and/or 'mac', depending on the user's platform. The code will check
    that a binary by this name exists for the user's distribution, and will only register this converter if one is
    found.
    """

    @classmethod
    def can_be_registered(cls) -> bool:
        """If a binary is required for this script, check that it exists for the user's OS/distribution. If one isn't
        required (`cls.required_bin` is None), also return True
        """
        if cls.required_bin is None:
            return True
        return bin_exists(cls.required_bin)

    def _convert(self):

        self.logger.debug(f"Performing conversion with ScriptFileConverter using script '{self.script}'")

        env = {"DIST": get_dist()}
        if self.required_bin is not None:
            env["BIN_PATH"] = get_bin_path(self.required_bin)

        process = subprocess.run(['sh', f'psdi_data_conversion/scripts/{self.script}', *self._get_script_args()],
                                 env=env, capture_output=True, text=True)

        self.out = process.stdout
        self.err = process.stderr

        if process.returncode != 0:
            self.logger.error(f"Conversion process completed with non-zero returncode {process.returncode}; aborting")
            self._abort_from_err()
        else:
            self.logger.debug("Conversion process completed successfully")

    def _get_script_args(self):
        """Get the list of arguments which will be passed to the script"""

        from_flags = self.data.get("from_flags", "")
        to_flags = self.data.get("from_flags", "")
        from_options = self.data.get("from_options", "")
        to_options = self.data.get("from_options", "")

        # Check that all user-provided input passes security checks
        for user_args in [from_flags, to_flags, from_options, to_options]:
            if not string_is_safe(user_args):
                raise FileConverterInputException(f"Provided argument '{user_args}' does not pass security check - it "
                                                  f"must match the regex {SAFE_STRING_RE.pattern}.", help=True)

        return ['--' + self.to_format_info.name, self.in_filename, self.out_filename, from_flags, to_flags,
                from_options, to_options]
