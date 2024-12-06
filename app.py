"""app.py

Version 1.0, 8th November 2024

This script acts as a server for the PSDI Data Conversion Service website.
"""

import hashlib
import os
import py.io
import json
import subprocess
from multiprocessing import Lock
from datetime import datetime
from openbabel import openbabel
from flask import Flask, request, render_template, abort, Response

# Maximum output file size in bytes
MEGABYTE = 1024*1024
MAX_FILE_SIZE = 1*MEGABYTE

# A lock to prevent multiple threads logging at the same time.
logLock = Lock()

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

# File to log any errors that occur
ERROR_LOG_FILENAME = "error_log.txt"
GLOBAL_ERROR_LOG = f"./{ERROR_LOG_FILENAME}"

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
        return convertFile('fileToUpload')
    else:
        # return http status code 405
        abort(405)


@app.route('/conv/', methods=['POST'])
def conv():
    """Convert file (cURL)
    """
    return convertFile('file')


def logErrorMessage(message, localErrorLog):
    """Report an error message in both the global and local error logs

    Parameters
    ----------
    message : str
        The error message
    localErrorLog : str
        Fully-qualified name of error log local to this process
    """
    for error_log in (GLOBAL_ERROR_LOG, localErrorLog):
        with open(error_log, 'a') as f:
            f.write(message)


def checkFileSize(inFilename, outFilename, localErrorLog):
    """Get file sizes, checking that output file isn't too large

    Parameters
    ----------
    inFilename : str
        Fully-qualified name of input file
    outFilename : str
        Fully-qualified name of output file
    localErrorLog : str
        Fully-qualified name of error log local to this process

    Returns
    -------
    inSize : int
        Size of input file in bytes
    outSize : int
        Size of output file in bytes
    """
    inSize = os.path.getsize(inFilename)
    outSize = os.path.getsize(outFilename)

    # Check that the output file doesn't exceed the maximum allowed size
    if outSize > MAX_FILE_SIZE:
        logErrorMessage(f"ERROR converting {inFilename} to {outFilename}: Output file exceeds maximum size.\n" +
                        f"Input file size is {inSize/MEGABYTE:.2f} MB; Output file size is {outSize/MEGABYTE:.2f} " +
                        f"MB; maximum output file size is {MAX_FILE_SIZE/MEGABYTE:.2f} MB.\n",
                        localErrorLog)

        # Delete output and input files
        os.remove(inFilename)
        os.remove(outFilename)

        abort(405)   # return http status code 405

    return inSize, outSize


def convertFile(file):
    """Convert the uploaded file to the required format, generating atomic coordinates if required
    """

    f = request.files[file]
    fname = f.filename.split(".")[0]  # E.g. ethane.mol --> ethane

    inFilename = 'static/uploads/' + f.filename

    f.save(inFilename)

    # Retrieve 'from' and 'to' file formats
    fromFormat = request.form['from']
    toFormat = request.form['to']

    converter = request.form['converter']

    outFilename = 'static/downloads/' + fname + '.' + toFormat

    localErrorLog = f"{DOWNLOAD_DIR}/{f.filename}-{fname}.{toFormat}.err"

    # If any previous error log exists, delete it
    if os.path.exists(localErrorLog):
        os.remove(localErrorLog)

    if converter == 'Open Babel':
        stdouterrOB = py.io.StdCaptureFD(in_=False)

        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats(fromFormat, toFormat)

        # Retrieve 'from' and 'to' option flags and arguments
        fromFlags = request.form['from_flags']
        toFlags = request.form['to_flags']
        fromArgFlags = request.form['from_arg_flags']
        toArgFlags = request.form['to_arg_flags']
        fromArgs = request.form['from_args']
        toArgs = request.form['to_args']

        # Add option flags and arguments as appropriate
        for char in fromFlags:
            obConversion.AddOption(char, obConversion.INOPTIONS)

        for char in toFlags:
            obConversion.AddOption(char, obConversion.OUTOPTIONS)

        readFlagsArgs = []
        writeFlagsArgs = []

        for char in fromArgFlags:
            index = fromArgs.find('£')
            arg = fromArgs[0:index]
            fromArgs = fromArgs[index + 1:len(fromArgs)]
            readFlagsArgs.append(char + "  " + arg)
            obConversion.AddOption(char, obConversion.INOPTIONS, arg)

        for char in toArgFlags:
            index = toArgs.find('£')
            arg = toArgs[0:index]
            toArgs = toArgs[index + 1:len(toArgs)]
            writeFlagsArgs.append(char + "  " + arg)
            obConversion.AddOption(char, obConversion.OUTOPTIONS, arg)

        # Read the file to be converted
        mol = openbabel.OBMol()
        obConversion.ReadFile(mol, inFilename)

        # Retrieve coordinate calculation type (Gen2D, Gen3D, neither)
        calcType = request.form['coordinates']

        option = 'N/A'

        # Calculate atomic coordinates
        if calcType != 'neither':
            # Retrieve coordinate calculation option (fastest, fast, medium, better, best)
            option = request.form['coordOption']

            gen = openbabel.OBOp.FindType(calcType)
            gen.Do(mol, option)

        # Write the converted file
        obConversion.WriteFile(mol, outFilename)

        out, err = stdouterrOB.reset()   # Grab stdout and stderr

        # Determine file sizes for logging purposes and check output isn't too large
        inSize, outSize = checkFileSize(inFilename, outFilename, localErrorLog)

        if file != 'file':  # Website only (i.e., not command line option)
            os.remove(inFilename)
            fromFormat = request.form['from_full']
            toFormat = request.form['to_full']
            quality = request.form['success']
        else:
            quality = getQuality(fromFormat, toFormat)

        if err.find('Error') > -1:
            logError(fromFormat, toFormat, converter, fname, calcType, option, fromFlags,
                     toFlags, readFlagsArgs, writeFlagsArgs, err, localErrorLog)
            stdouterrOB.done()
            abort(405)  # return http status code 405
        else:
            log(fromFormat, toFormat, converter, fname, calcType, option, fromFlags,
                toFlags, readFlagsArgs, writeFlagsArgs, quality, out, err)

        stdouterrOB.done()
    elif request.form['converter'] == 'Atomsk':
        atomsk = subprocess.run(['sh', 'atomsk.sh', f.filename, fname, toFormat], capture_output=True, text=True)

        out = atomsk.stdout
        err = atomsk.stderr

        # Determine file sizes for logging purposes and check output isn't too large
        inSize, outSize = checkFileSize(inFilename, outFilename, localErrorLog)

        if file != 'file':   # Website only (i.e., not command line option)
            os.remove(inFilename)
            fromFormat = request.form['from_full']
            toFormat = request.form['to_full']
            quality = request.form['success']
        else:
            quality = getQuality(fromFormat, toFormat)

        if err.find('Error') > -1:
            logErrorAto(fromFormat, toFormat, converter, fname, err, localErrorLog)
            abort(405)   # return http status code 405
        else:
            logAto(fromFormat, toFormat, converter, fname, quality, out, err)

    appendToLogFile("conversions", {
        "datetime": getDateTime(),
        "fromFormat": fromFormat,
        "toFormat": toFormat,
        "converter": converter,
        "fname": fname,
        "inSize": inSize,
        "outSize": outSize,
        "quality": quality,
    })

    return '\nConverting from ' + fname + '.' + fromFormat + ' to ' + fname + '.' + toFormat + '\n'


@app.route('/feedback/', methods=['POST'])
def feedback():
    """Take feedback data from the web app and log it
    """

    try:

        entry = {
            "datetime": getDateTime(),
        }

        report = json.loads(request.form['data'])

        for key in ["type", "missing", "reason", "from", "to"]:
            if key in report:
                entry[key] = str(report[key])

        appendToLogFile("feedback", entry)

        return Response(status=201)

    except Exception:

        return Response(status=400)


# Query the JSON file to obtain conversion quality
def getQuality(fromExt, toExt):
    """Query the JSON file to obtain conversion quality
    """

    try:

        # Load JSON file.
        with open("static/data/data.json") as datafile:
            data = json.load(datafile)

        fromFormat = [d for d in data["formats"] if d["extension"] == fromExt]
        toFormat = [d for d in data["formats"] if d["extension"] == toExt]
        openBabel = [d for d in data["converters"] if d["name"] == "Open Babel"]

        openBabelId = openBabel[0]["id"]
        fromId = fromFormat[0]["id"]
        toId = toFormat[0]["id"]

        convertsTo = [d for d in data["converts_to"] if
                      d["converters_id"] == openBabelId and d["in_id"] == fromId and d["out_id"] == toId]

        return convertsTo[0]["degree_of_success"]

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
def deleteFile():
    """Delete file (cURL)
    """
    os.remove(request.form['filepath'])
    return 'Server-side file ' + request.form['filepath'] + ' deleted\n'


def getDate():
    """Retrieve current date as a string

    Returns
    -------
    str
        Current date in the format YYYY-MM-DD
    """
    today = datetime.today()
    return str(today.year) + '-' + format(today.month) + '-' + format(today.day)


def getTime():
    """Retrieve current time as a string

    Returns
    -------
    str
        Current time in the format HH:MM:SS
    """
    today = datetime.today()
    return format(today.hour) + ':' + format(today.minute) + ':' + format(today.second)


def getDateTime():
    """Retrieve current date and time as a string

    Returns
    -------
    str
        Current date and time in the format YYYY-MM-DD HH:MM:SS
    """
    return getDate() + ' ' + getTime()


def log(fromFormat, toFormat, converter, fname, calcType, option, fromFlags, toFlags, readFlagsArgs, writeFlagsArgs,
        quality, out, err):
    """Write Open Babel conversion information to server-side file, ready for downloading to user

    Parameters
    ----------
    fromFormat : _type_
        _description_
    toFormat : _type_
        _description_
    converter : _type_
        _description_
    fname : _type_
        _description_
    calcType : _type_
        _description_
    option : _type_
        _description_
    fromFlags : _type_
        _description_
    toFlags : _type_
        _description_
    readFlagsArgs : _type_
        _description_
    writeFlagsArgs : _type_
        _description_
    quality : _type_
        _description_
    out : _type_
        _description_
    err : _type_
        _description_
    """

    message = (createMessage(fname, fromFormat, toFormat, converter, calcType, option, fromFlags, toFlags,
                             readFlagsArgs, writeFlagsArgs) +
               'Quality:           ' + quality + '\n' +
               'Success:           Assuming that the data provided was of the correct format, the conversion\n' +
               '                   was successful (to the best of our knowledge) subject to any warnings below.\n' +
               out + '\n' + err + '\n')

    f = open('static/downloads/' + fname + '.log.txt', 'w')
    f.write(message)
    f.close()


def logAto(fromFormat, toFormat, converter, fname, quality, out, err):
    """Write Atomsk conversion information to server-side file, ready for downloading to user

    Parameters
    ----------
    fromFormat : _type_
        _description_
    toFormat : _type_
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

    message = createMessageStart(fname, fromFormat, toFormat, converter) + \
        'Quality:           ' + quality + '\n' \
        'Success:           Assuming that the data provided was of the correct format, the conversion\n' \
        '                   was successful (to the best of our knowledge) subject to any warnings below.\n' + \
        out + '\n' + err + '\n'

    f = open('static/downloads/' + fname + '.log.txt', 'w')
    f.write(message)
    f.close()


def logError(fromFormat, toFormat, converter, fname, calcType, option, fromFlags, toFlags, readFlagsArgs,
             writeFlagsArgs, err, localErrorLog):
    """Write Open Babel conversion error information to server-side log file

    Parameters
    ----------
    fromFormat : _type_
        _description_
    toFormat : _type_
        _description_
    converter : _type_
        _description_
    fname : _type_
        _description_
    calcType : _type_
        _description_
    option : _type_
        _description_
    fromFlags : _type_
        _description_
    toFlags : _type_
        _description_
    readFlagsArgs : _type_
        _description_
    writeFlagsArgs : _type_
        _description_
    err : _type_
        _description_
    localErrorLog : _type_
        _description_
    """
    message = createMessage(fname, fromFormat, toFormat, converter, calcType, option,
                            fromFlags, toFlags, readFlagsArgs, writeFlagsArgs) + err + '\n'
    logErrorMessage(message, localErrorLog)


def logErrorAto(fromFormat, toFormat, converter, fname, err, localErrorLog):
    """Write Atomsk conversion error information to server-side log file

    Parameters
    ----------
    fromFormat : _type_
        _description_
    toFormat : _type_
        _description_
    converter : _type_
        _description_
    fname : _type_
        _description_
    err : _type_
        _description_
    localErrorLog : _type_
        _description_
    """
    message = createMessage(fname, fromFormat, toFormat, converter) + err + '\n'
    logErrorMessage(message, localErrorLog)


def createMessage(fname, fromFormat, toFormat, converter, calcType, option, fromFlags, toFlags, readFlagsArgs,
                  writeFlagsArgs):
    """Create message for log files

    Parameters
    ----------
    fname : _type_
        _description_
    fromFormat : _type_
        _description_
    toFormat : _type_
        _description_
    converter : _type_
        _description_
    calcType : _type_
        _description_
    option : _type_
        _description_
    fromFlags : _type_
        _description_
    toFlags : _type_
        _description_
    readFlagsArgs : _type_
        _description_
    writeFlagsArgs : _type_
        _description_

    Returns
    -------
    str
        The message for log files
    """
    str = ''

    if calcType == 'neither':
        str = 'Coord. gen.:       none\n'
    else:
        str += 'Coord. gen.:       ' + calcType + '\n'

    str += 'Coord. option:     ' + option + '\n'

    if fromFlags == '':
        str += 'Read options:      none\n'
    else:
        str += 'Read options:      ' + fromFlags + '\n'

    if toFlags == '':
        str += 'Write options:     none\n'
    else:
        str += 'Write options:     ' + toFlags + '\n'

    if len(readFlagsArgs) == 0:
        str += 'Read opts + args:  none\n'
    else:
        headingAdded = False

        for pair in readFlagsArgs:
            if not headingAdded:
                str += 'Read opts + args:  ' + pair + '\n'
                headingAdded = True
            else:
                str += '                   ' + pair + '\n'

    if len(writeFlagsArgs) == 0:
        str += 'Write opts + args: none\n'
    else:
        headingAdded = False

        for pair in writeFlagsArgs:
            if not headingAdded:
                str += 'Write opts + args: ' + pair + '\n'
                headingAdded = True
            else:
                str += '                   ' + pair + '\n'

    return createMessageStart(fname, fromFormat, toFormat, converter) + str


def createMessageStart(fname, fromFormat, toFormat, converter):
    """Create beginning of message for log files

    Parameters
    ----------
    fname : _type_
        _description_
    fromFormat : _type_
        _description_
    toFormat : _type_
        _description_
    converter : _type_
        _description_

    Returns
    -------
    str
        The beginning of a message for log files, containing generic information about what was trying to be done
    """
    return 'Date:              ' + getDate() + '\n' \
           'Time:              ' + getTime() + '\n' \
           'File name:         ' + fname + '\n' \
           'From:              ' + fromFormat + '\n' \
           'To:                ' + toFormat + '\n' \
           'Converter:         ' + converter + '\n'


def appendToLogFile(logName, data):
    """Append data to a log file

    Parameters
    ----------
    logName : _type_
        _description_
    data : _type_
        _description_
    """

    return

    # logLock.acquire()

    # try:
    #     if re.match(r"^[a-z]+$", logName):
    #         with open(f"var/{logName}.log", "a") as logFile:
    #             logFile.write(f"{json.dumps(data)}\n")

    # finally:
    #     logLock.release()


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
        message = '[' + getDateTime() + '] ' + request.args['data'] + '\n'

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
