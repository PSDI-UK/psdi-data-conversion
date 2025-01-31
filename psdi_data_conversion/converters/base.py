"""@file psdi_data_conversion/converters/base.py

Created 2025-01-23 by Bryan Gillis.

Base class and information for file format converters
"""


from copy import deepcopy
import json
import logging
from collections.abc import Callable
import os
import subprocess
import abc

import traceback
from typing import Any

from psdi_data_conversion import constants as const, log_utility

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
    pass


class FileConverterAbortException(FileConverterException):
    """Class representing an exception triggered by a call to abort a file conversion
    """

    def __init__(self, status_code, *args):
        super().__init__(*args)
        self.status_code = status_code


class FileConverterInputException(FileConverterException):
    """Exception class to represent errors encountered with input parameters for the data conversion script.
    """
    pass


if HTTPException is not None:
    l_abort_exceptions = (HTTPException, FileConverterAbortException)
else:
    l_abort_exceptions = (FileConverterAbortException,)


def abort_raise(status_code):
    """Callback for aborting during a file conversion, which re-raises any raised exceptions
    """
    raise FileConverterAbortException(status_code)


class FileConverter:
    """Class to handle conversion of files from one type to another
    """

    # Name of the converter - must be overridden in each subclass to name each converter uniquely
    name: str | None = None

    # General info about the converter - should be overridden in each subclass to describe the converter
    info: str | None = None

    # List of flags allowed for the converter (flags are arguments that are set by being present, and don't require a
    # value specified - e.g. "-v" to enable verbose mode) - should be overridden with a tuple of tuples containing the
    # flag names and help texts for them
    allowed_flags: tuple[tuple[str, str], ...] | None = None

    # List of options allowed for the converter (options are arguments that take one or more values, e.g. "-o out.txt")
    # - should be overridden with a tuple of tuples containing the option names and help texts for them
    allowed_options: tuple[tuple[str, str], ...] | None = None

    def __init__(self,
                 filename: str,
                 to_format: str,
                 from_format: str | None = None,
                 data: dict[str, Any] | None = None,
                 abort_callback: Callable[[int], None] = abort_raise,
                 use_envvars=False,
                 upload_dir=const.DEFAULT_UPLOAD_DIR,
                 download_dir=const.DEFAULT_DOWNLOAD_DIR,
                 max_file_size=const.DEFAULT_MAX_FILE_SIZE,
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
            The maximum allowed file size for input/output files, in MB. If 0, will be unlimited. Default 0 (unlimited)
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
        """

        # Set member variables directly from input
        self.in_filename = filename
        self.to_format = to_format
        self.abort_callback = abort_callback
        self.upload_dir = upload_dir
        self.download_dir = download_dir
        self.max_file_size = max_file_size*const.MEGABYTE
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

        # Remove any leading periods from to/from_format
        if self.to_format.startswith("."):
            self.to_format = self.to_format[1:]
        if self.from_format.startswith("."):
            self.from_format = self.from_format[1:]

        # Set placeholders for member variables which will be set when conversion is run
        self.in_size: int | None = None
        self.out_size: int | None = None
        self.out: str | None = None
        self.err: str | None = None
        self.quality: str | None = None

        # Set values from envvars if desired
        if use_envvars:
            # Get the maximum allowed size from the envvar for it
            ev_max_file_size = os.environ.get(const.MAX_FILESIZE_ENVVAR)
            if ev_max_file_size is not None:
                self.max_file_size = float(ev_max_file_size)*const.MEGABYTE

        # Create directory 'uploads' if not extant.
        if not os.path.exists(self.upload_dir):
            os.makedirs(self.upload_dir, exist_ok=True)

        # Create directory 'downloads' if not extant.
        if not os.path.exists(self.download_dir):
            os.makedirs(self.download_dir, exist_ok=True)

        self.local_filename = os.path.split(self.in_filename)[1]
        self.filename_base = os.path.splitext(self.local_filename)[0]
        self.out_filename = f"{self.download_dir}/{self.filename_base}.{self.to_format}"

        # Set up files to log to
        self._setup_loggers()

        self.logger.debug("Finished FileConverter initialisation")

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
            elif self.log_mode == const.LOG_SIMPLE or self.log_mode == const.LOG_FULL:
                self._local_logger_level = const.DEFAULT_LOCAL_LOGGER_LEVEL
                self._stdout_output_level = logging.ERROR
            else:
                raise FileConverterInputException(f"ERROR: Unrecognised logging option: {self.log_mode}. Allowed "
                                                  f"options are: {const.L_ALLOWED_LOG_MODES}")
        if self.log_mode == const.LOG_FULL:
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
        self.output_log = os.path.join(self.download_dir, f"{self.filename_base}{const.OUTPUT_LOG_EXT}")

        # If any previous log exists, delete it
        if os.path.exists(self.output_log):
            os.remove(self.output_log)

        write_mode = "w" if self.refresh_local_log else "a"
        # Set up loggers - one for general-purpose log_utility, and one just for what we want to output to the user
        self.logger = log_utility.set_up_data_conversion_logger(local_log_file=self.output_log,
                                                                local_logger_level=self._local_logger_level,
                                                                local_logger_raw_output=False,
                                                                mode=write_mode)

        self.logger.debug(f"Set up server-style logging, with user logging at level {self._local_logger_level}")

    def run(self):
        """Run the file conversion
        """

        try:
            self.logger.debug("Starting file conversion")
            self._convert()

            self.logger.debug("Finished file conversion; performing cleanup tasks")
            self._finish_convert()
        except Exception as e:
            if isinstance(e, l_abort_exceptions):
                # Don't catch a deliberate abort; let it pass through
                self.logger.error(f"Unexpected exception raised of type '{type(e)}' with message: {str(e)}")
                raise
            self.logger.error(f"Exception triggering an abort was raised. Exception was type '{type(e)}', with "
                              f"message: {str(e)}")
            self._abort(message=f"The application encountered an error:\n{traceback.format_exc()}")

        return ('\nConverting from ' + self.filename_base + '.' + self.from_format + ' to ' + self.filename_base +
                '.' + self.to_format + '\n')

    def _abort(self, status_code=const.STATUS_CODE_GENERAL, message=None):
        """Abort the conversion, reporting the desired message to the user at the top of the output

        Parameters
        ----------
        status_code : int
            The HTTP status code to exit with. Default is 422: Unprocessable Content
        message : str | None
            If provided, this message will be logged in the user output log at the top of the file. This should
            typically explain the reason the process failed
        """

        # Remove the input and output files if they exist
        if self.delete_input:
            self.logger.debug(f"Cleaning up input file {self.in_filename}")
            try:
                os.remove(self.in_filename)
            except FileNotFoundError:
                pass
        try:
            os.remove(self.out_filename)
        except FileNotFoundError:
            self.logger.debug("Application aborting; no output file found to clean up")
        else:
            self.logger.debug(f"Application aborting, so cleaning up output file {self.out_filename}")

        if message:
            # If we're adding a message in server mode, read in any prior logs, clear the log, write the message, then
            # write the prior logs
            if self.log_file is None:
                self.logger.debug("Adding abort message to the top of the output log so it will be the first thing "
                                  "read by the user")
                prior_output_log = open(self.output_log, "r").read()
                os.remove(self.output_log)
                with open(self.output_log, "w") as fo:
                    fo.write(message + "\n")
                    fo.write(prior_output_log)

            # Note this message in the dev logger as well
            self.logger.error(message)

        self.abort_callback(status_code)

    def _abort_from_err(self):
        """Call an abort after a call to the converter has completed, but it's returned an error. Create a message for
        the logger including this error and other relevant information.
        """
        self.logger.error(self._create_message_start() +
                          self._create_message() +
                          self.out + '\n' +
                          self.err)
        self._abort(message=self.err)

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
                f"From:              {self.from_format}\n"
                f"To:                {self.to_format}\n"
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

    def _check_file_size_and_status(self):
        """Get file sizes, checking that output file isn't too large

        Returns
        -------
        in_size : int
            Size of input file in bytes
        out_size : int
            Size of output file in bytes
        """
        in_size = os.path.getsize(os.path.realpath(self.in_filename))
        try:
            out_size = os.path.getsize(os.path.realpath(self.out_filename))
        except FileNotFoundError:
            # Something went wrong and the output file doesn't exist
            err_message = f"Expected output file {self.out_filename} does not exist."
            self.logger.error(err_message)
            self.err += f"ERROR: {err_message}\n"
            self._abort_from_err()

        # Check that the output file doesn't exceed the maximum allowed size
        if self.max_file_size > 0 and out_size > self.max_file_size:

            self._abort(const.STATUS_CODE_SIZE,
                        f"ERROR converting {os.path.basename(self.in_filename)} to " +
                        os.path.basename(self.out_filename) + ": "
                        f"Output file exceeds maximum size.\nInput file size is "
                        f"{in_size/const.MEGABYTE:.2f} MB; Output file size is {out_size/const.MEGABYTE:.2f} "
                        f"MB; maximum output file size is {self.max_file_size/const.MEGABYTE:.2f} MB.\n")
        self.logger.debug(f"Output file found to have size {out_size/const.MEGABYTE:.2f} MB")

        return in_size, out_size

    def get_quality(self):
        """Query the JSON file to obtain conversion quality
        """

        try:

            # Load JSON file.
            with open("static/data/data.json") as datafile:
                data = json.load(datafile)

            from_format = [d for d in data["formats"] if d["extension"] == self.from_format]
            to_format = [d for d in data["formats"] if d["extension"] == self.to_format]
            open_babel = [d for d in data["converters"] if d["name"] == "Open Babel"]

            open_babel_id = open_babel[0]["id"]
            from_id = from_format[0]["id"]
            to_id = to_format[0]["id"]

            converts_to = [d for d in data["converts_to"] if
                           d["converters_id"] == open_babel_id and d["in_id"] == from_id and d["out_id"] == to_id]

            return converts_to[0]["degree_of_success"]

        except Exception as e:

            self.logger.warning(f"Unable to determine conversion quality from {self.from_format} to {self.to_format}. "
                                f"Received exception of type '{type(e)}' and message: {str(e)}")

            return "unknown"

    def _finish_convert(self):
        """Run final common steps to clean up a conversion and log success or abort due to an error
        """

        self.in_size, self.out_size = self._check_file_size_and_status()

        if self.delete_input:
            os.remove(self.in_filename)
        if "from_full" in self.data:
            self.from_format = self.data["from_full"]
        if "to_full" in self.data:
            self.to_format = self.data["from_full"]
        if "success" in self.data:
            self.quality = self.data["success"]
        else:
            self.quality = self.get_quality()

        self._log_success()

    @ abc.abstractmethod
    def _convert(self):
        """Run the conversion with the desired converter
        """
        pass


class ScriptFileConverter(FileConverter):
    """File Converter specialized to run a shell script to call the converter
    """

    script: str | None = None

    def _convert(self):

        self.logger.debug(f"Performing conversion with ScriptFileConverter using script '{self.script}'")

        if "from_flags" not in self.data:
            self.data["from_flags"] = ""
        if "to_flags" not in self.data:
            self.data["to_flags"] = ""

        process = subprocess.run(['sh', f'psdi_data_conversion/scripts/{self.script}',
                                  self.in_filename, self.out_filename, self.data["from_flags"], self.data["to_flags"]],
                                 capture_output=True, text=True)

        self.out = process.stdout
        self.err = process.stderr

        if process.returncode != 0:
            self.logger.error(f"Conversion process completed with non-zero returncode {process.returncode}; aborting")
            self._abort_from_err()
        else:
            self.logger.debug("Conversion process completed successfully")
