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
from flask import request, abort

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


def convert_file(file):
    """Convert the uploaded file to the required format, generating atomic coordinates if required
    """

    f = request.files[file]
    filename_base = f.filename.split(".")[0]  # E.g. ethane.mol --> ethane

    in_filename = 'static/uploads/' + f.filename

    f.save(in_filename)

    # Retrieve 'from' and 'to' file formats
    from_format = request.form['from']
    to_format = request.form['to']

    converter = request.form['converter']

    out_filename = f"{DOWNLOAD_DIR}/{filename_base}.{to_format}"

    # Set up files to log to
    local_log_base = f"{DOWNLOAD_DIR}/{f.filename}-{filename_base}.{to_format}"
    local_log = f"{local_log_base}.log"
    local_error = f"{local_log_base}.err"
    output_log = f"{DOWNLOAD_DIR}/{filename_base}.log.txt"

    # If any previous local logs exist, delete them
    if os.path.exists(local_log):
        os.remove(local_log)
    if os.path.exists(local_error):
        os.remove(local_error)
    if os.path.exists(output_log):
        os.remove(output_log)

    # Set up loggers - one for general-purpose log_utility, and one just for what we want to output to the user
    log_utility.setUpDataConversionLogger(local_log_file=local_log,
                                          local_logger_level=log_utility.DEFAULT_LOCAL_LOGGER_LEVEL,
                                          stdout_output_level=logging.INFO)
    log_utility.setUpDataConversionLogger(name="output",
                                          local_log_file=output_log,
                                          local_logger_raw_output=True,
                                          extra_loggers=[(local_log, log_utility.DEFAULT_LOCAL_LOGGER_LEVEL, False)])

    if converter == 'Open Babel':
        stdouterr_ob = py.io.StdCaptureFD(in_=False)

        ob_conversion = openbabel.OBConversion()
        ob_conversion.SetInAndOutFormats(from_format, to_format)

        # Retrieve 'from' and 'to' option flags and arguments
        from_flags = request.form['from_flags']
        to_flags = request.form['to_flags']
        from_arg_flags = request.form['from_arg_flags']
        to_arg_flags = request.form['to_arg_flags']
        from_args = request.form['from_args']
        to_args = request.form['to_args']

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
        ob_conversion.ReadFile(mol, in_filename)

        # Retrieve coordinate calculation type (Gen2D, Gen3D, neither)
        calc_type = request.form['coordinates']

        option = 'N/A'

        # Calculate atomic coordinates
        if calc_type != 'neither':
            # Retrieve coordinate calculation option (fastest, fast, medium, better, best)
            option = request.form['coordOption']

            gen = openbabel.OBOp.FindType(calc_type)
            gen.Do(mol, option)

        # Write the converted file
        ob_conversion.WriteFile(mol, out_filename)

        out, err = stdouterr_ob.reset()   # Grab stdout and stderr

        # Determine file sizes for logging purposes and check output isn't too large
        in_size, out_size = check_file_size(in_filename, out_filename)

        if file != 'file':  # Website only (i.e., not command line option)
            os.remove(in_filename)
            from_format = request.form['from_full']
            to_format = request.form['to_full']
            quality = request.form['success']
        else:
            quality = get_quality(from_format, to_format)

        if err.find('Error') > -1:
            log_utility.log_error(from_format, to_format, converter, filename_base, calc_type, option, from_flags,
                                  to_flags, read_flags_args, write_flags_args, err)
            stdouterr_ob.done()
            abort(405)  # return http status code 405
        else:
            log_utility.log(from_format, to_format, converter, filename_base, calc_type, option, from_flags,
                            to_flags, read_flags_args, write_flags_args, quality, out, err)

        stdouterr_ob.done()
    elif request.form['converter'] == 'Atomsk':
        atomsk = subprocess.run(['sh', 'atomsk.sh', f.filename, filename_base,
                                to_format], capture_output=True, text=True)

        out = atomsk.stdout
        err = atomsk.stderr

        # Determine file sizes for logging purposes and check output isn't too large
        in_size, out_size = check_file_size(in_filename, out_filename)

        if file != 'file':   # Website only (i.e., not command line option)
            os.remove(in_filename)
            from_format = request.form['from_full']
            to_format = request.form['to_full']
            quality = request.form['success']
        else:
            quality = get_quality(from_format, to_format)

        if err.find('Error') > -1:
            log_utility.log_error_ato(from_format, to_format, converter, filename_base, err)
            abort(405)   # return http status code 405
        else:
            log_utility.log_ato(from_format, to_format, converter, filename_base, quality, out, err)

    log_utility.append_to_log_file("conversions", {
        "datetime": log_utility.get_date_time(),
        "fromFormat": from_format,
        "toFormat": to_format,
        "converter": converter,
        "fname": filename_base,
        "inSize": in_size,
        "outSize": out_size,
        "quality": quality,
    })

    return '\nConverting from ' + filename_base + '.' + from_format + ' to ' + filename_base + '.' + to_format + '\n'


# Query the JSON file to obtain conversion quality
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
