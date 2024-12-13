"""@file psdi-data-conversion/psdi_data_conversion/converter.py

Created 2024-12-10 by Bryan Gillis.

Class and functions to perform file conversion
"""

import json
import logging
import os
import py.io
import subprocess
from openbabel import openbabel
from flask import abort
from werkzeug.datastructures import FileStorage

from psdi_data_conversion import log_utility

# Maximum output file size in bytes
MEGABYTE = 1024*1024
MAX_FILE_SIZE = 1*MEGABYTE

# Create directory 'uploads' if not extant.
UPLOAD_DIR = './static/uploads'
if not os.path.exists(UPLOAD_DIR):
    os.makedirs(UPLOAD_DIR, exist_ok=True)

# Create directory 'downloads' if not extant.
DOWNLOAD_DIR = './static/downloads'
if not os.path.exists(DOWNLOAD_DIR):
    os.makedirs(DOWNLOAD_DIR, exist_ok=True)


def check_file_size(in_filename, out_filename):
    """Get file sizes, checking that output file isn't too large

    Parameters
    ----------
    in_filename : str
        Fully-qualified name of input file
    out_filename : str
        Fully-qualified name of output file

    Returns
    -------
    in_size : int
        Size of input file in bytes
    out_size : int
        Size of output file in bytes
    """
    in_size = os.path.getsize(in_filename)
    out_size = os.path.getsize(out_filename)

    # Check that the output file doesn't exceed the maximum allowed size
    if out_size > MAX_FILE_SIZE:
        log_utility.getDataConversionLogger("output").error(
            f"ERROR converting {os.path.basename(in_filename)} to {os.path.basename(out_filename)}: "
            f"Output file exceeds maximum size.\nInput file size is "
            f"{in_size/MEGABYTE:.2f} MB; Output file size is {out_size/MEGABYTE:.2f} "
            f"MB; maximum output file size is {MAX_FILE_SIZE/MEGABYTE:.2f} MB.\n")

        # Delete output and input files
        os.remove(in_filename)
        os.remove(out_filename)

        abort(405)   # return http status code 405

    return in_size, out_size


class FileConverter:
    """Class to handle conversion of files from one type to another
    """

    def __init__(self,
                 files: dict[str, FileStorage],
                 form: dict[str, str],
                 file_to_convert: str):
        """Initialize the object, storing needed data and setting up loggers.

        Parameters
        ----------
        files : ImmutableMultiDict[str, FileStorage]
            The file dict provided by Flask at `request.files`
        form : ImmutableMultiDict[str, str]
            The form dict provided by Flask at `request.form`
        file_to_convert : str
            The key for the file in the `files` dict to convert
        """
        self.files = files
        self.form = form
        self.file_to_convert = file_to_convert

        self.f = self.files[self.file_to_convert]
        self.filename_base = self.f.filename.split(".")[0]  # E.g. ethane.mol --> ethane

        self.in_filename = f"{UPLOAD_DIR}/{self.f.filename}"

        self.f.save(self.in_filename)

        # Retrieve 'from' and 'to' file formats
        self.from_format = self.form['from']
        self.to_format = self.form['to']

        self.converter = self.form['converter']

        self.out_filename = f"{DOWNLOAD_DIR}/{self.filename_base}.{self.to_format}"

        # Set up files to log to
        local_log_base = f"{DOWNLOAD_DIR}/{self.f.filename}-{self.filename_base}.{self.to_format}"
        local_log = f"{local_log_base}.log"
        local_error = f"{local_log_base}.err"
        output_log = f"{DOWNLOAD_DIR}/{self.filename_base}.log.txt"

        # If any previous local logs exist, delete them
        if os.path.exists(local_log):
            os.remove(local_log)
        if os.path.exists(local_error):
            os.remove(local_error)
        if os.path.exists(output_log):
            os.remove(output_log)

        # Set up loggers - one for general-purpose log_utility, and one just for what we want to output to the user
        self.logger = log_utility.setUpDataConversionLogger(local_log_file=local_log,
                                                            local_logger_level=log_utility.DEFAULT_LOCAL_LOGGER_LEVEL,
                                                            stdout_output_level=logging.INFO)
        self.output_logger = log_utility.setUpDataConversionLogger(name="output",
                                                                   local_log_file=output_log,
                                                                   local_logger_raw_output=True,
                                                                   extra_loggers=[(
                                                                       local_log,
                                                                       log_utility.DEFAULT_LOCAL_LOGGER_LEVEL,
                                                                       False)])

    def __call__(self):
        """Run the file conversion
        """

        if self.converter == 'Open Babel':
            self._convert_ob()
        elif self.converter == 'Atomsk':
            self._convert_ato()
        else:
            self.output_logger.error(f"ERROR: Unknown logger '{self.converter}' requested")
            abort(405)

        self.in_size, self.out_size = check_file_size(self.in_filename, self.out_filename)

        self._append_to_log_file("conversions")

        return ('\nConverting from ' + self.filename_base + '.' + self.from_format + ' to ' + self.filename_base +
                '.' + self.to_format + '\n')

    def _convert_ob(self):
        stdouterr_ob = py.io.StdCaptureFD(in_=False)

        ob_conversion = openbabel.OBConversion()
        ob_conversion.SetInAndOutFormats(self.from_format, self.to_format)

        # Retrieve 'from' and 'to' option flags and arguments
        from_flags = self.form['from_flags']
        to_flags = self.form['to_flags']
        from_arg_flags = self.form['from_arg_flags']
        to_arg_flags = self.form['to_arg_flags']
        from_args = self.form['from_args']
        to_args = self.form['to_args']

        # Add option flags and arguments as appropriate
        for char in from_flags:
            ob_conversion.AddOption(char, ob_conversion.INOPTIONS)

        for char in to_flags:
            ob_conversion.AddOption(char, ob_conversion.OUTOPTIONS)

        read_flags_args = []
        write_flags_args = []

        for char in from_arg_flags:
            index = from_args.find('£')
            arg = from_args[0:index]
            from_args = from_args[index + 1:len(from_args)]
            read_flags_args.append(char + "  " + arg)
            ob_conversion.AddOption(char, ob_conversion.INOPTIONS, arg)

        for char in to_arg_flags:
            index = to_args.find('£')
            arg = to_args[0:index]
            to_args = to_args[index + 1:len(to_args)]
            write_flags_args.append(char + "  " + arg)
            ob_conversion.AddOption(char, ob_conversion.OUTOPTIONS, arg)

        # Read the file to be converted
        mol = openbabel.OBMol()
        ob_conversion.ReadFile(mol, self.in_filename)

        # Retrieve coordinate calculation type (Gen2D, Gen3D, neither)
        calc_type = self.form['coordinates']

        option = 'N/A'

        # Calculate atomic coordinates
        if calc_type != 'neither':
            # Retrieve coordinate calculation option (fastest, fast, medium, better, best)
            option = self.form['coordOption']

            gen = openbabel.OBOp.FindType(calc_type)
            gen.Do(mol, option)

        # Write the converted file
        ob_conversion.WriteFile(mol, self.out_filename)

        self.out, self.err = stdouterr_ob.reset()   # Grab stdout and stderr

        if self.file != 'file':  # Website only (i.e., not command line option)
            os.remove(self.in_filename)
            self.from_format = self.form['from_full']
            self.to_format = self.form['to_full']
            self.quality = self.form['success']
        else:
            self.quality = self.get_quality(self.from_format, self.to_format)

        if self.err.find('Error') > -1:
            self._log_error_ob()
            stdouterr_ob.done()
            abort(405)  # return http status code 405
        else:
            self._log_success_ob()
            stdouterr_ob.done()

    def _convert_ato(self):
        atomsk = subprocess.run(['sh', 'atomsk.sh', self.f.filename, self.filename_base,
                                self.to_format], capture_output=True, text=True)

        self.out = atomsk.stdout
        self.err = atomsk.stderr

        if self.file != 'file':   # Website only (i.e., not command line option)
            os.remove(self.in_filename)
            self.from_format = self.form['from_full']
            self.to_format = self.form['to_full']
            self.quality = self.form['success']
        else:
            self.quality = self.get_quality(self.from_format, self.to_format)

        if self.err.find('Error') > -1:
            self._log_error_ato()
            abort(405)   # return http status code 405
        else:
            self._log_success_ato()

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
