"""@file psdi_data_conversion/converters/base.py

Created 2025-01-23 by Bryan Gillis.

Base class and information for file format converters
"""


import json
import logging
from collections.abc import Callable
import os
import subprocess
import sys
import abc

import traceback

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

    name: str | None = None

    def __init__(self,
                 filename: str,
                 form: dict[str, str],
                 abort_callback: Callable[[int], None] = abort_raise,
                 use_envvars=False,
                 upload_dir=const.DEFAULT_UPLOAD_DIR,
                 download_dir=const.DEFAULT_DOWNLOAD_DIR,
                 max_file_size=const.DEFAULT_MAX_FILE_SIZE,
                 log_file: str | None = None,
                 log_mode=const.LOG_FULL,
                 delete_input=False):
        """Initialize the object, storing needed data and setting up loggers.

        Parameters
        ----------
        filename : str
            The filename of the input file to be converted, either relative to current directory or fully-qualified
        form : ImmutableMultiDict[str, str]
            The form dict provided by Flask at `request.form`
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
        delete_input : bool
            Whether or not to delete input files after conversion, default False
        """

        # Set member variables directly from input
        self.in_filename = filename
        self.form = form
        self.abort_callback = abort_callback
        self.upload_dir = upload_dir
        self.download_dir = download_dir
        self.max_file_size = max_file_size*const.MEGABYTE
        self.log_file = log_file
        self.log_mode = log_mode
        self.delete_input = delete_input

        # Set member variables from dict values in input
        self.from_format = self.form['from']
        self.to_format = self.form['to']

        # Set placeholders for member variables which will be set when conversion is run
        self.in_size: int | None = None
        self.out_size: int | None = None
        self.out: str | None = None
        self.err: str | None = None
        self.quality: str | None = None

        # Set placeholders for member variables used only by OB conversion
        self.from_flags: str | None = None
        self.to_flags: str | None = None
        self.read_flags_args: list[str] | None = None
        self.write_flags_args: list[str] | None = None
        self.calc_type: str | None = None
        self.option: str | None = None

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

    def _setup_loggers(self):
        """Run at init to set up loggers for this object.
        """

        # Determine level to log at based on quiet status
        if self.log_mode == const.LOG_NONE:
            self._local_logger_level = None
            self._stdout_output_level = logging.ERROR
        elif self.log_mode == const.LOG_STDOUT:
            self._local_logger_level = None
            self._stdout_output_level = logging.INFO
        elif self.log_mode == const.LOG_SIMPLE or self.log_mode == const.LOG_FULL:
            self._local_logger_level = const.DEFAULT_LOCAL_LOGGER_LEVEL
            self._stdout_output_level = logging.ERROR
            if self.log_mode == const.LOG_FULL:
                return self._setup_server_loggers()
        else:
            raise FileConverterInputException(f"ERROR: Unrecognised logging option: {self.log_mode}. Allowed options "
                                              f"are: {const.L_ALLOWED_LOG_MODES}", file=sys.stderr)

        self.output_log = self.log_file

        self.logger = log_utility.set_up_data_conversion_logger(local_log_file=self.log_file,
                                                                local_logger_level=self._local_logger_level,
                                                                stdout_output_level=self._stdout_output_level,
                                                                suppress_global_handler=True)
        self.output_logger = self.logger

    def _setup_server_loggers(self):
        """Run at init to set up loggers for this object in server-style execution
        """
        local_log_base = f"{self.download_dir}/{self.local_filename}-{self.filename_base}.{self.to_format}"
        local_log = f"{local_log_base}{const.LOCAL_LOG_EXT}"
        self.output_log = f"{self.download_dir}/{self.filename_base}{const.OUTPUT_LOG_EXT}"

        # If any previous local logs exist, delete them
        if os.path.exists(local_log):
            os.remove(local_log)
        if os.path.exists(self.output_log):
            os.remove(self.output_log)

        # Set up loggers - one for general-purpose log_utility, and one just for what we want to output to the user
        self.logger = log_utility.set_up_data_conversion_logger(local_log_file=local_log,
                                                                local_logger_level=self._local_logger_level,
                                                                stdout_output_level=self._stdout_output_level)
        self.output_logger = log_utility.set_up_data_conversion_logger(name="output",
                                                                       local_log_file=self.output_log,
                                                                       local_logger_raw_output=True,
                                                                       extra_loggers=[(local_log,
                                                                                       self._local_logger_level,
                                                                                       False)])

    def run(self):
        """Run the file conversion
        """

        try:
            self._convert()
            self._finish_convert()
        except Exception as e:
            if isinstance(e, l_abort_exceptions):
                # Don't catch a deliberate abort; let it pass through
                raise
            self._abort(message=f"The application encountered an unexpected error:\n{traceback.format_exc()}")

        self._append_to_log_file("conversions")

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
            try:
                os.remove(self.in_filename)
            except FileNotFoundError:
                pass
        try:
            os.remove(self.out_filename)
        except FileNotFoundError:
            pass

        if message:
            # If we're adding a message in server mode, read in any prior logs, clear the log, write the message, then
            # write the prior logs
            if self.log_file is None:
                prior_output_log = open(self.output_log, "r").read()
                os.remove(self.output_log)
                with open(self.output_log, "w") as fo:
                    fo.write(message + "\n")
                    fo.write(prior_output_log)

            # Note this message in the dev logger as well
            self.logger.error(message)

        self.abort_callback(status_code)

    def _abort_from_err(self):
        """Write conversion error information to server-side log file and abort the conversion
        """
        self.output_logger.info(self._create_message())
        self._abort(message=self.err)

    def _create_message(self):

        message = ''

        if self.calc_type == 'neither':
            message = 'Coord. gen.:       none\n'
        elif self.calc_type:
            message += 'Coord. gen.:       ' + self.calc_type + '\n'

        if self.option:
            message += 'Coord. option:     ' + self.option + '\n'

        if self.from_flags == '':
            message += 'Read options:      none\n'
        elif self.from_flags:
            message += 'Read options:      ' + self.from_flags + '\n'

        if self.to_flags == '':
            message += 'Write options:     none\n'
        elif self.to_flags:
            message += 'Write options:     ' + self.to_flags + '\n'

        if self.read_flags_args is None:
            pass
        elif len(self.read_flags_args) == 0:
            message += 'Read opts + args:  none\n'
        else:
            heading_added = False

            for pair in self.read_flags_args:
                if not heading_added:
                    message += 'Read opts + args:  ' + pair + '\n'
                    heading_added = True
                else:
                    message += '                   ' + pair + '\n'

        if self.write_flags_args is None:
            pass
        elif len(self.write_flags_args) == 0:
            message += 'Write opts + args: none\n'
        else:
            heading_added = False

            for pair in self.write_flags_args:
                if not heading_added:
                    message += 'Write opts + args: ' + pair + '\n'
                    heading_added = True
                else:
                    message += '                   ' + pair + '\n'

        return self._create_message_start() + message

    def _create_message_start(self):
        """Create beginning of message for log files

        Returns
        -------
        str
            The beginning of a message for log files, containing generic information about what was trying to be done
        """
        return ('\n'
                'File name:         ' + self.filename_base + '\n'
                'From:              ' + self.from_format + '\n'
                'To:                ' + self.to_format + '\n'
                'Converter:         ' + self.name + '\n')

    def _log_success(self):
        """Write conversion information to server-side file, ready for downloading to user
        """

        message = (self._create_message() +
                   'Quality:           ' + self.quality + '\n'
                   'Success:           Assuming that the data provided was of the correct format, the conversion\n'
                   '                   was successful (to the best of our knowledge) subject to any warnings below.\n' +
                   self.out + '\n' + self.err).strip() + '\n'

        self.output_logger.info(message)

    def _append_to_log_file(self, log_name):
        """Append data to a log file

        Parameters
        ----------
        log_name : _type_
            _description_
        """

        data = {
            "datetime": log_utility.get_date_time(),
            "fromFormat": self.from_format,
            "toFormat": self.to_format,
            "converter": self.name,
            "fname": self.filename_base,
            "inSize": self.in_size,
            "outSize": self.out_size,
            "quality": self.quality,
        }

        if (os.environ.get('ENABLE_DCS_LOG') is not None):
            print(json.dumps(data))

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
            self.err += f"ERROR: Expected output file {self.out_filename} does not exist.\n"
            self._abort_from_err()

        # Check that the output file doesn't exceed the maximum allowed size
        if self.max_file_size > 0 and out_size > self.max_file_size:

            self._abort(const.STATUS_CODE_SIZE,
                        f"ERROR converting {os.path.basename(self.in_filename)} to " +
                        os.path.basename(self.out_filename) + ": "
                        f"Output file exceeds maximum size.\nInput file size is "
                        f"{in_size/const.MEGABYTE:.2f} MB; Output file size is {out_size/const.MEGABYTE:.2f} "
                        f"MB; maximum output file size is {self.max_file_size/const.MEGABYTE:.2f} MB.\n")

        return in_size, out_size

    @staticmethod
    def get_quality(from_ext, to_ext):
        """Query the JSON file to obtain conversion quality
        """

        try:

            # Load JSON file.
            with open("static/data/data.json") as datafile:
                data = json.load(datafile)

            from_format = [d for d in data["formats"] if d["extension"] == from_ext]
            to_format = [d for d in data["formats"] if d["extension"] == to_ext]
            open_babel = [d for d in data["converters"] if d["name"] == "Open Babel"]

            open_babel_id = open_babel[0]["id"]
            from_id = from_format[0]["id"]
            to_id = to_format[0]["id"]

            converts_to = [d for d in data["converts_to"] if
                           d["converters_id"] == open_babel_id and d["in_id"] == from_id and d["out_id"] == to_id]

            return converts_to[0]["degree_of_success"]

        except Exception:

            return "unknown"

    def _finish_convert(self):
        """Run final common steps to clean up a conversion and log success or abort due to an error
        """

        self.in_size, self.out_size = self._check_file_size_and_status()

        if self.delete_input:
            os.remove(self.in_filename)
        if "from_full" in self.form:
            self.from_format = self.form["from_full"]
        if "to_full" in self.form:
            self.to_format = self.form["from_full"]
        if "success" in self.form:
            self.quality = self.form["success"]
        else:
            self.quality = self.get_quality(self.from_format, self.to_format)

        self._log_success()

    @abc.abstractmethod
    def _convert(self):
        """Run the conversion with the desired converter
        """
        pass


class ScriptFileConverter(FileConverter):
    """File Converter specialized to run a shell script to call the converter
    """

    script: str | None = None

    def _convert(self):

        if not self.from_flags:
            self.from_flags = ""
        if not self.to_flags:
            self.to_flags = ""

        process = subprocess.run(['sh', f'psdi_data_conversion/scripts/{self.script}',
                                 self.in_filename, self.out_filename, self.from_flags, self.to_flags],
                                 capture_output=True, text=True)

        self.out = process.stdout
        self.err = process.stderr

        if process.returncode != 0:
            self._abort_from_err()
