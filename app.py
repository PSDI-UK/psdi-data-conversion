
#   app.py
#   Version 1.0, 8th November 2024

#   This script acts as a server for the PSDI Data Conversion Service website.

import hashlib, os, py.io, json, re, subprocess, time
#import hashlib, os, glob, psycopg2, py.io, json, re $$$$$$$$$$$$$$$$$$$$$$$$$$$ DELETE $$$$$$$$$$$$$$$$$$$$$$$$$$$$
from multiprocessing import Lock
#from psycopg2 import sql $$$$$$$$$$$$$$$$$$$$$$$$$$$$ DELETE $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
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
upDir = './static/uploads'
if not os.path.exists(upDir) :
    os.mkdir(upDir)

# Create directory 'downloads' if not extant.
downDir = './static/downloads'
if not os.path.exists(downDir) :
    os.mkdir(downDir)

app = Flask(__name__)

# Return the web page along with the token.
@app.route('/')
def website() :
    data = [{'token': token}]
    return render_template("index.htm", data=data)

# Convert file to a different format and save to folder 'downloads'. Delete original file.
# Note that downloading is achieved in format.js
@app.route('/convert/', methods=['POST'])
def convert() :
    if request.form['token'] == token and token != '' :
        return convertFile('fileToUpload')
    else :
        # return http status code 405
        abort(405)

# Convert file (cURL)
@app.route('/conv/', methods=['POST'])
def conv() :
    return convertFile('file')

# Get file sizes, checking that output file isn't too large
def checkFileSize(inFilename, outFilename):
    inSize = os.path.getsize(inFilename)
    outSize = os.path.getsize(outFilename)

    # Check that the output file doesn't exceed the maximum allowed size
    if outSize > MAX_FILE_SIZE:
        with open('error_log.txt', 'a') as f:
            f.write("ERROR: Output file exceeds maximum size.\n" +
                f"File size is {outSize/MEGABYTE:.2f} MB; maximum size is {MAX_FILE_SIZE/MEGABYTE:.2f} MB.")

        # Delete output and input files
        os.remove(inFilename)
        os.remove(outFilename)

        abort(405)   # return http status code 405

    return inSize, outSize

# Convert the uploaded file to the required format, generating atomic coordinates if required
def convertFile(file) :

    f = request.files[file]
        
    inFilename = 'static/uploads/' + f.filename
    outFilename = 'static/downloads/' + fname + '.' + toFormat
    
    f.save(inFilename)
    fname = f.filename.split(".")[0]  # E.g. ethane.mol --> ethane

    # Retrieve 'from' and 'to' file formats
    fromFormat = request.form['from']
    toFormat = request.form['to']

    converter = request.form['converter']

    if converter == 'Open Babel' :
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
        for char in fromFlags :
            obConversion.AddOption(char, obConversion.INOPTIONS)

        for char in toFlags :
            obConversion.AddOption(char, obConversion.OUTOPTIONS)

        readFlagsArgs = []
        writeFlagsArgs = []

        for char in fromArgFlags :
            index = fromArgs.find('£')
            arg = fromArgs[0:index]
            fromArgs = fromArgs[index + 1:len(fromArgs)]
            readFlagsArgs.append(char + "  " + arg)
            obConversion.AddOption(char, obConversion.INOPTIONS, arg)

        for char in toArgFlags :
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
        if calcType != 'neither' :
            # Retrieve coordinate calculation option (fastest, fast, medium, better, best)
            option = request.form['coordOption']

            gen = openbabel.OBOp.FindType(calcType)
            gen.Do(mol, option)

        # Write the converted file
        obConversion.WriteFile(mol, outFilename)

        out,err = stdouterrOB.reset()   # Grab stdout and stderr

        # Determine file sizes for logging purposes and check output isn't too large
        inSize, outSize = checkFileSize(inFilename, outFilename)

        if file != 'file' : # Website only (i.e., not command line option)
            os.remove(inFilename)
            fromFormat = request.form['from_full']
            toFormat = request.form['to_full']
            quality = request.form['success']
        else :
            quality = getQuality(fromFormat, toFormat)

        if err.find('Error') > -1 :
            error_log(fromFormat, toFormat, converter, fname, calcType, option, fromFlags, toFlags, readFlagsArgs, writeFlagsArgs, err)
            stdouterrOB.done()
            abort(405) # return http status code 405
        else :
            log(fromFormat, toFormat, converter, fname, calcType, option, fromFlags, toFlags, readFlagsArgs, writeFlagsArgs, quality, out, err)

        stdouterrOB.done()
    elif request.form['converter'] == 'Atomsk' :
        atomsk = subprocess.run(['sh', 'atomsk.sh', f.filename, fname, toFormat], capture_output = True, text = True)

        out = atomsk.stdout
        err = atomsk.stderr

        # Determine file sizes for logging purposes and check output isn't too large
        inSize, outSize = checkFileSize(inFilename, outFilename)

        if file != 'file' :   # Website only (i.e., not command line option)
            os.remove(inFilename)
            fromFormat = request.form['from_full']
            toFormat = request.form['to_full']
            quality = request.form['success']
        else :
            quality = getQuality(fromFormat, toFormat)

        if err.find('Error') > -1 :
            error_log_ato(fromFormat, toFormat, converter, fname, err)
            abort(405)   # return http status code 405
        else :
            log_ato(fromFormat, toFormat, converter, fname, quality, out, err)

    appendToLogFile("conversions", {
        "datetime": get_date_time(),
        "fromFormat": fromFormat,
        "toFormat": toFormat,
        "converter": converter,
        "fname": fname,
        "inSize": inSize,
        "outSize": outSize,
        "quality": quality,
    })

    return '\nConverting from ' + fname + '.' + fromFormat + ' to ' + fname + '.' + toFormat +'\n'

# Take feedback data from the web app and log it.
@app.route('/feedback/', methods=['POST'])
def feedback() :

    try:

        entry = {
            "datetime": get_date_time(),
        }

        report = json.loads(request.form['data'])

        for key in ["type", "missing", "reason", "from", "to"]:
            if key in report:
                entry[key] = str(report[key])

        appendToLogFile("feedback", entry)

        return Response(status=201)

    except:

        return Response(status=400)


# Query the JSON file to obtain conversion quality
def getQuality(fromExt, toExt):

    try:

        # Load JSON file.
        with open("static/data/data.json") as datafile:
            data = json.load(datafile)

        from_format = [d for d in data["formats"] if d["extension"] == fromExt]
        to_format = [d for d in data["formats"] if d["extension"] == toExt]
        open_babel = [d for d in data["converters"] if d["name"] == "Open Babel"]

        open_babel_id = open_babel[0]["id"]
        from_id = from_format[0]["id"]
        to_id = to_format[0]["id"]

        converts_to = [d for d in data["converts_to"] if
            d["converters_id"] == open_babel_id and d["in_id"] == from_id and d["out_id"] == to_id]

        return converts_to[0]["degree_of_success"]

    except:

        return "unknown"

# Delete files in folder 'downloads'
@app.route('/delete/', methods=['POST'])
def delete() :
    os.remove('static/downloads/' + request.form['filename'])
    os.remove('static/downloads/' + request.form['logname'])

    return 'okay'

# Delete file (cURL)
@app.route('/del/', methods=['POST'])
def deleteFile() :
    os.remove(request.form['filepath'])
    return 'Server-side file ' + request.form['filepath'] + ' deleted\n'

# Retrieve current date as a string.
def get_date() :
    today = datetime.today()    
    return str(today.year) + '-' + format(today.month) + '-' + format(today.day)

# Retrieve current time as a string.
def get_time() :
    today = datetime.today()    
    return format(today.hour) + ':' + format(today.minute) + ':' + format(today.second)

# Retrieve current date and time as a string.
def get_date_time() :
    return get_date() + ' ' + get_time()

# Write Open Babel conversion information to server-side file, ready for downloading to user.
def log(fromFormat, toFormat, converter, fname, calcType, option, fromFlags, toFlags, readFlagsArgs, writeFlagsArgs, quality, out, err) :
    message = create_message(fname, fromFormat, toFormat, converter, calcType, option, fromFlags, toFlags, readFlagsArgs, writeFlagsArgs) + \
              'Quality:           ' + quality + '\n' \
              'Success:           Assuming that the data provided was of the correct format, the conversion\n' \
              '                   was successful (to the best of our knowledge) subject to any warnings below.\n' + \
              out + '\n' + err + '\n'

    f = open('static/downloads/' + fname + '.log.txt', 'w')
    f.write(message)
    f.close()

# Write Atomsk conversion information to server-side file, ready for downloading to user.
def log_ato(fromFormat, toFormat, converter, fname, quality, out, err) :
    message = create_message_start(fname, fromFormat, toFormat, converter) + \
              'Quality:           ' + quality + '\n' \
              'Success:           Assuming that the data provided was of the correct format, the conversion\n' \
              '                   was successful (to the best of our knowledge) subject to any warnings below.\n' + \
              out + '\n' + err + '\n'

    f = open('static/downloads/' + fname + '.log.txt', 'w')
    f.write(message)
    f.close()

# Write Open Babel conversion error information to server-side log file.
def error_log(fromFormat, toFormat, converter, fname, calcType, option, fromFlags, toFlags, readFlagsArgs, writeFlagsArgs, err) :
    message = create_message(fname, fromFormat, toFormat, converter, calcType, option, fromFlags, toFlags, readFlagsArgs, writeFlagsArgs) + err + '\n'

    f = open('error_log.txt', 'a')
    f.write(message)
    f.close()

# Write Atomsk conversion error information to server-side log file.
def error_log_ato(fromFormat, toFormat, converter, fname, err) :
    message = create_message(fname, fromFormat, toFormat, converter) + err + '\n'

    f = open('error_log.txt', 'a')
    f.write(message)
    f.close()

# Create message for log files.
def create_message(fname, fromFormat, toFormat, converter, calcType, option, fromFlags, toFlags, readFlagsArgs, writeFlagsArgs) :
    str = ''

    if calcType == 'neither' :
        str = 'Coord. gen.:       none\n'
    else :
        str += 'Coord. gen.:       ' + calcType + '\n'

    str += 'Coord. option:     ' + option + '\n'

    if fromFlags == '' :
        str += 'Read options:      none\n'
    else :
        str += 'Read options:      ' + fromFlags  + '\n'

    if toFlags == '' :
        str += 'Write options:     none\n'
    else :
        str += 'Write options:     ' + toFlags  + '\n'

    if len(readFlagsArgs) == 0 :
        str += 'Read opts + args:  none\n'
    else :
        headingAdded = False

        for pair in readFlagsArgs :
            if not headingAdded :
                str += 'Read opts + args:  ' + pair + '\n'
                headingAdded = True
            else :
                str += '                   ' + pair + '\n'

    if len(writeFlagsArgs) == 0 :
        str += 'Write opts + args: none\n'
    else :
        headingAdded = False

        for pair in writeFlagsArgs :
            if not headingAdded :
                str += 'Write opts + args: ' + pair + '\n'
                headingAdded = True
            else :
                str += '                   ' + pair + '\n'

    return create_message_start(fname, fromFormat, toFormat, converter) + str

# Create beginning of message for log files
def create_message_start(fname, fromFormat, toFormat, converter) :
    return 'Date:              ' + get_date() + '\n' \
           'Time:              ' + get_time() + '\n' \
           'File name:         ' + fname + '\n' \
           'From:              ' + fromFormat + '\n' \
           'To:                ' + toFormat + '\n' \
           'Converter:         ' + converter + '\n'

# Append data to a log file.
def appendToLogFile(log_name, data):

    return

    # logLock.acquire()

    # try:
    #     if re.match(r"^[a-z]+$", log_name):
    #         with open(f"var/{log_name}.log", "a") as log_file:
    #             log_file.write(f"{json.dumps(data)}\n")

    # finally:
    #     logLock.release()

# Check that the incoming token matches the one sent to the user (should mostly prevent spambots).
# Write date- and time-stamped user input to server-side file 'user_responses'.
# $$$$$$$$$$ Retained in case direct logging is required in the future. $$$$$$$$$$
@app.route('/data/', methods=['GET'])
def data() :
    if request.args['token'] == token and token != '' :
        message = '[' + get_date_time() + '] ' + request.args['data'] + '\n'

        f = open("user_responses", "a")
        f.write(message)
        f.close()

        return 'okay'
    else :
        # return http status code 405
        abort(405)

# Ensure that an element of date or time (month, day, hours, minutes or seconds) always has two digits.
def format(time) :
    num = str(time)

    if len(num) == 1 :
        return '0' + num
    else :
        return num
