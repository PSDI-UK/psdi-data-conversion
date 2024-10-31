
#   app.py
#   Version 1.0, 14th October 2024

#   This script acts as a server for the PSDI Data Conversion Service website.

import hashlib, os, glob, psycopg2, py.io, json
from psycopg2 import sql
from datetime import datetime
from openbabel import openbabel
from flask import Flask, request, render_template, abort

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

# Convert the uploaded file to the required format, generating atomic coordinates if required
def convertFile(file) :
    f = request.files[file]
    f.save('static/uploads/' + f.filename)
    fname = f.filename.split(".")[0]  # E.g. ethane.mol --> ethane

    # Retrieve 'from' and 'to' file formats
    fromFormat = request.form['from']
    toFormat = request.form['to']

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
    obConversion.ReadFile(mol, 'static/uploads/' + f.filename)

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
    obConversion.WriteFile(mol, 'static/downloads/' + fname + '.' + toFormat)

    # Determine values for logging purposes from here
    inSize = os.path.getsize('static/uploads/' + f.filename)
    outSize = os.path.getsize('static/downloads/' + fname + '.' + toFormat)

    # Website only (i.e., not command line option)
    if file != 'file' :
        os.remove('static/uploads/' + f.filename)
        fromFormat = request.form['from_full']
        toFormat = request.form['to_full']
        quality = request.form['success']
    else :
        quality = getQuality(fromFormat, toFormat)

    converter = 'Open Babel'           # $$$$$$$$$$ TODO: Replace hard coding when more than one converter option. $$$$$$$$$$
    out,err = stdouterrOB.reset()

    if err.find('Error') > -1 :
        error_log(fromFormat, toFormat, converter, fname, calcType, option, fromFlags, toFlags, readFlagsArgs, writeFlagsArgs, err)
        stdouterrOB.done()

        # return http status code 405
        abort(405)
    else :
        log(fromFormat, toFormat, converter, fname, calcType, option, fromFlags, toFlags, readFlagsArgs, writeFlagsArgs, quality, out, err)

    stdouterrOB.done()

    query = sql.SQL("INSERT INTO {table} ({ident}, {datetime}, {convFrom}, {convTo}, {conv}, {filename}, {fromSize}, {toSize}, {success}) " \
                    "VALUES ((SELECT COALESCE(MAX({ident}), 0) FROM {table}) + 1, %s, %s, %s, %s, %s, %s, %s, %s)").format(
                table=sql.Identifier('conversion_log'),
                ident=sql.Identifier('id'),
                datetime=sql.Identifier('date'),
                convFrom=sql.Identifier('convert_from'),
                convTo=sql.Identifier('convert_to'),
                conv=sql.Identifier('converter'),
                filename=sql.Identifier('file_name'),
                fromSize=sql.Identifier('from_file_size'),
                toSize=sql.Identifier('to_file_size'),
                success=sql.Identifier('success_status'))
    values = [get_date_time(), fromFormat, toFormat, converter, fname, inSize, outSize, quality]
    logConversion(query, values)

    return '\nConverting from ' + fname + '.' + fromFormat + ' to ' + fname + '.' + toFormat +'\n'

# Query the PostgreSQL database
@app.route('/query/', methods=['POST'])
def query() :
    if request.form['token'] == token and token != '' :
        # Establish a connection with the PostgreSQL database
        try:
            #db_conn = psycopg2.connect(database="format")
            db_conn = psycopg2.connect(dbname="psdi", user="psdi", password="SharkCat1", \
                                       host="psdi.postgres.database.azure.com", port=5432)
        except psycopg2.DatabaseError as Error:
            print(f"Connection to database failed. {Error}")

        # Query database
        with db_conn.cursor() as cursor:
            cursor.execute(request.form['data'])

            if request.form['data'].startswith('SELECT') :
                results = cursor.fetchall()
            else :
                db_conn.commit()

        # Close connection to database
        if db_conn:
            db_conn.close()
 
        if request.form['data'].startswith('SELECT') :
            # Construct a string from the array and return it
            ansArray = []

            for row in results :
                line = ''

                for el in row :
                    line += "£" + str(el)

                ansArray.append(line)

            return '$'.join(ansArray)
        else :
            return 'true'
    else :
        # return http status code 405
        abort(405)

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

# Write conversion information to server-side file, ready for downloading to user.
def log(fromFormat, toFormat, converter, fname, calcType, option, fromFlags, toFlags, readFlagsArgs, writeFlagsArgs, quality, out, err) :
    message = create_message(fname, fromFormat, toFormat, converter, calcType, option, fromFlags, toFlags, readFlagsArgs, writeFlagsArgs) + \
              'Quality:           ' + quality + '\n' \
              'Success:           Assuming that the data provided was of the correct format, the conversion\n' \
              '                   was successful (to the best of our knowledge) subject to any warnings below.\n' + \
              out + '\n' + err + '\n'

    f = open('static/downloads/' + fname + '.log.txt', 'w')
    f.write(message)
    f.close()

# Write conversion error information to server-side log file.
def error_log(fromFormat, toFormat, converter, fname, calcType, option, fromFlags, toFlags, readFlagsArgs, writeFlagsArgs, err) :
    message = create_message(fname, fromFormat, toFormat, converter, calcType, option, fromFlags, toFlags, readFlagsArgs, writeFlagsArgs) + err + '\n'

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

    return 'Date:              ' + get_date() + '\n' \
           'Time:              ' + get_time() + '\n' \
           'File name:         ' + fname + '\n' \
           'From:              ' + fromFormat + '\n' \
           'To:                ' + toFormat + '\n' \
           'Converter:         ' + converter + '\n' + str

# Write to database table Conversion_Log.
def logConversion(query, values) :
    # Establish a connection with the PostgreSQL database
    try:
        #db_conn = psycopg2.connect(database="format")
        db_conn = psycopg2.connect(dbname="psdi", user="psdi", password="SharkCat1", \
                                   host="psdi.postgres.database.azure.com", port=5432)
    except psycopg2.DatabaseError as Error:
        print(f"Connection to database failed. {Error}")

    # Query database
    with db_conn.cursor() as cursor:
        cursor.execute(query, values)
        db_conn.commit()

    # Close connection to database
    if db_conn:
        db_conn.close()
 
    return 'true'

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
