"""app.py

Version 1.0, 8th November 2024

This script acts as a server for the PSDI Data Conversion Service website.
"""

import hashlib
import os
import py.io
import json
import subprocess
from datetime import datetime
from openbabel import openbabel
from flask import Flask, request, render_template, abort, Response

from psdi_data_conversion import logging

logger = logging.getLogger()

# Maximum output file size in bytes
MEGABYTE = 1024*1024
MAX_FILE_SIZE = 1*MEGABYTE

# Create a token by hashing the current date and time.
dt = str(datetime.now())
token = hashlib.md5(dt.encode('utf8')).hexdigest()

# Create directory 'uploads' if not extant.
UPLOAD_DIR = './static/uploads'
if not os.path.exists(UPLOAD_DIR):
    os.mkdir(UPLOAD_DIR)

# Create directory 'downloads' if not extant.
DOWNLOAD_DIR = './static/downloads'
if not os.path.exists(DOWNLOAD_DIR):
    os.mkdir(DOWNLOAD_DIR)

app = Flask(__name__)


@app.route('/')
def website():
    """Return the web page along with the token
    """
    data = [{'token': token}]
    return render_template("index.htm", data=data)


@app.route('/convert/', methods=['POST'])
def convert():
    """Convert file to a different format and save to folder 'downloads'. Delete original file. Note that downloading is
    achieved in format.js
    """
    if request.form['token'] == token and token != '':
        return convert_file('fileToUpload')
    else:
        # return http status code 405
        abort(405)


@app.route('/conv/', methods=['POST'])
def conv():
    """Convert file (cURL)
    """
    return convert_file('file')


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
        logger.error(f"ERROR converting {in_filename} to {out_filename}: Output file exceeds maximum size.\n"
                     f"Input file size is {in_size/MEGABYTE:.2f} MB; Output file size is {out_size/MEGABYTE:.2f} "
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

    out_filename = 'static/downloads/' + filename_base + '.' + to_format

    local_error_log = f"{DOWNLOAD_DIR}/{f.filename}-{filename_base}.{to_format}.err"

    # If any previous error log exists, delete it
    if os.path.exists(local_error_log):
        os.remove(local_error_log)

    global logger
    logger = logging.getLogger(__name__, local_error_log)

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
            logError(from_format, to_format, converter, filename_base, calc_type, option, from_flags,
                     to_flags, read_flags_args, write_flags_args, err)
            stdouterr_ob.done()
            abort(405)  # return http status code 405
        else:
            log(from_format, to_format, converter, filename_base, calc_type, option, from_flags,
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
            log_error_ato(from_format, to_format, converter, filename_base, err)
            abort(405)   # return http status code 405
        else:
            log_ato(from_format, to_format, converter, filename_base, quality, out, err)

    append_to_log_file("conversions", {
        "datetime": get_date_time(),
        "fromFormat": from_format,
        "toFormat": to_format,
        "converter": converter,
        "fname": filename_base,
        "inSize": in_size,
        "outSize": out_size,
        "quality": quality,
    })

    return '\nConverting from ' + filename_base + '.' + from_format + ' to ' + filename_base + '.' + to_format + '\n'


@app.route('/feedback/', methods=['POST'])
def feedback():
    """Take feedback data from the web app and log it
    """

    try:

        entry = {
            "datetime": get_date_time(),
        }

        report = json.loads(request.form['data'])

        for key in ["type", "missing", "reason", "from", "to"]:
            if key in report:
                entry[key] = str(report[key])

        append_to_log_file("feedback", entry)

        return Response(status=201)

    except Exception:

        return Response(status=400)


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


@app.route('/delete/', methods=['POST'])
def delete():
    """Delete files in folder 'downloads'
    """
    os.remove('static/downloads/' + request.form['filename'])
    os.remove('static/downloads/' + request.form['logname'])

    return 'okay'


@app.route('/del/', methods=['POST'])
def delete_file():
    """Delete file (cURL)
    """
    os.remove(request.form['filepath'])
    return 'Server-side file ' + request.form['filepath'] + ' deleted\n'


def get_date():
    """Retrieve current date as a string

    Returns
    -------
    str
        Current date in the format YYYY-MM-DD
    """
    today = datetime.today()
    return str(today.year) + '-' + format(today.month) + '-' + format(today.day)


def get_time():
    """Retrieve current time as a string

    Returns
    -------
    str
        Current time in the format HH:MM:SS
    """
    today = datetime.today()
    return format(today.hour) + ':' + format(today.minute) + ':' + format(today.second)


def get_date_time():
    """Retrieve current date and time as a string

    Returns
    -------
    str
        Current date and time in the format YYYY-MM-DD HH:MM:SS
    """
    return get_date() + ' ' + get_time()


def log(from_format, to_format, converter, fname, calc_type, option, from_flags, to_flags, read_flags_args,
        write_flags_args, quality, out, err):
    """Write Open Babel conversion information to server-side file, ready for downloading to user

    Parameters
    ----------
    from_format : _type_
        _description_
    to_format : _type_
        _description_
    converter : _type_
        _description_
    fname : _type_
        _description_
    calc_type : _type_
        _description_
    option : _type_
        _description_
    from_flags : _type_
        _description_
    to_flags : _type_
        _description_
    read_flags_args : _type_
        _description_
    write_flags_args : _type_
        _description_
    quality : _type_
        _description_
    out : _type_
        _description_
    err : _type_
        _description_
    """

    message = (create_message(fname, from_format, to_format, converter, calc_type, option, from_flags, to_flags,
                              read_flags_args, write_flags_args) +
               'Quality:           ' + quality + '\n'
               'Success:           Assuming that the data provided was of the correct format, the conversion\n'
               '                   was successful (to the best of our knowledge) subject to any warnings below.\n' +
               out + '\n' + err + '\n')

    f = open('static/downloads/' + fname + '.log.txt', 'w')
    f.write(message)
    f.close()


def log_ato(from_format, to_format, converter, fname, quality, out, err):
    """Write Atomsk conversion information to server-side file, ready for downloading to user

    Parameters
    ----------
    from_format : _type_
        _description_
    to_format : _type_
        _description_
    converter : _type_
        _description_
    fname : _type_
        _description_
    quality : _type_
        _description_
    out : _type_
        _description_
    err : _type_
        _description_
    """

    message = (create_message_start(fname, from_format, to_format, converter) +
               'Quality:           ' + quality + '\n'
               'Success:           Assuming that the data provided was of the correct format, the conversion\n'
               '                   was successful (to the best of our knowledge) subject to any warnings below.\n' +
               out + '\n' + err + '\n')

    f = open('static/downloads/' + fname + '.log.txt', 'w')
    f.write(message)
    f.close()


def logError(from_format, to_format, converter, fname, calc_type, option, from_flags, to_flags, read_flags_args,
             write_flags_args, err):
    """Write Open Babel conversion error information to server-side log file

    Parameters
    ----------
    from_format : _type_
        _description_
    to_format : _type_
        _description_
    converter : _type_
        _description_
    fname : _type_
        _description_
    calc_type : _type_
        _description_
    option : _type_
        _description_
    from_flags : _type_
        _description_
    to_flags : _type_
        _description_
    read_flags_args : _type_
        _description_
    write_flags_args : _type_
        _description_
    err : _type_
        _description_
    """
    message = create_message(fname, from_format, to_format, converter, calc_type, option,
                             from_flags, to_flags, read_flags_args, write_flags_args) + err + '\n'
    logger.error(message)


def log_error_ato(from_format, to_format, converter, fname, err):
    """Write Atomsk conversion error information to server-side log file

    Parameters
    ----------
    from_format : _type_
        _description_
    to_format : _type_
        _description_
    converter : _type_
        _description_
    fname : _type_
        _description_
    err : _type_
        _description_
    """
    message = create_message(fname, from_format, to_format, converter) + err + '\n'
    logger.error(message)


def create_message(fname, from_format, to_format, converter, calc_type, option, from_flags, to_flags, read_flags_args,
                   write_flags_args):
    """Create message for log files

    Parameters
    ----------
    fname : _type_
        _description_
    from_format : _type_
        _description_
    to_format : _type_
        _description_
    converter : _type_
        _description_
    calc_type : _type_
        _description_
    option : _type_
        _description_
    from_flags : _type_
        _description_
    to_flags : _type_
        _description_
    read_flags_args : _type_
        _description_
    write_flags_args : _type_
        _description_

    Returns
    -------
    message : str
        The message for log files
    """
    message = ''

    if calc_type == 'neither':
        message = 'Coord. gen.:       none\n'
    else:
        message += 'Coord. gen.:       ' + calc_type + '\n'

    message += 'Coord. option:     ' + option + '\n'

    if from_flags == '':
        message += 'Read options:      none\n'
    else:
        message += 'Read options:      ' + from_flags + '\n'

    if to_flags == '':
        message += 'Write options:     none\n'
    else:
        message += 'Write options:     ' + to_flags + '\n'

    if len(read_flags_args) == 0:
        message += 'Read opts + args:  none\n'
    else:
        heading_added = False

        for pair in read_flags_args:
            if not heading_added:
                message += 'Read opts + args:  ' + pair + '\n'
                heading_added = True
            else:
                message += '                   ' + pair + '\n'

    if len(write_flags_args) == 0:
        message += 'Write opts + args: none\n'
    else:
        heading_added = False

        for pair in write_flags_args:
            if not heading_added:
                message += 'Write opts + args: ' + pair + '\n'
                heading_added = True
            else:
                message += '                   ' + pair + '\n'

    return create_message_start(fname, from_format, to_format, converter) + message


def create_message_start(fname, from_format, to_format, converter):
    """Create beginning of message for log files

    Parameters
    ----------
    fname : _type_
        _description_
    from_format : _type_
        _description_
    to_format : _type_
        _description_
    converter : _type_
        _description_

    Returns
    -------
    str
        The beginning of a message for log files, containing generic information about what was trying to be done
    """
    return ('Date:              ' + get_date() + '\n'
            'Time:              ' + get_time() + '\n'
            'File name:         ' + fname + '\n'
            'From:              ' + from_format + '\n'
            'To:                ' + to_format + '\n'
            'Converter:         ' + converter + '\n')


def append_to_log_file(log_name, data):
    """Append data to a log file

    Parameters
    ----------
    log_name : _type_
        _description_
    data : _type_
        _description_
    """

    if (os.environ.get('ENABLE_DCS_LOG') is not None):
        print(json.dumps(data))


@app.route('/data/', methods=['GET'])
def data():
    """Check that the incoming token matches the one sent to the user (should mostly prevent spambots). Write date- and
    time-stamped user input to server-side file 'user_responses'.

    $$$$$$$$$$ Retained in case direct logging is required in the future. $$$$$$$$$$

    Returns
    -------
    str
        Output status - 'okay' if exited successfuly
    """
    if request.args['token'] == token and token != '':
        message = '[' + get_date_time() + '] ' + request.args['data'] + '\n'

        f = open("user_responses", "a")
        f.write(message)
        f.close()

        return 'okay'
    else:
        # return http status code 405
        abort(405)


def format(time):
    """Ensure that an element of date or time (month, day, hours, minutes or seconds) always has two digits.

    Parameters
    ----------
    time : str or int
        Digit(s) indicating date or month

    Returns
    -------
    str
        2-digit value indicating date or month
    """
    num = str(time)

    if len(num) == 1:
        return '0' + num
    else:
        return num
