
#   app.py
#   Version 1.0, 21st May 2024

#   This script acts as a server for the PSDI Data Conversion Tool website.

import hashlib, os, glob, psycopg2, py.io
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

# Delete any files remaining in folder 'downloads'
@app.route('/delete/', methods=['POST'])
def delete() :
    os.remove('static/downloads/' + request.form['filename'])
    os.remove('static/downloads/' + request.form['logname'])

    return 'okay'

# Convert file to a different format and save to folder 'downloads'. Delete original file.
# Note that downloading is achieved in format.js
@app.route('/convert/', methods=['POST'])
def convert() :
    if request.form['token'] == token and token != '' :
        f = request.files['fileToUpload']
        f.save('static/uploads/' + f.filename)
        fname = f.filename.split(".")[0]  # E.g. ethane.mol --> ethane

        # Retrieve 'from' and 'to' file formats.
        fromFormat = request.form['from']
        toFormat = request.form['to']

        stdouterrOB = py.io.StdCaptureFD(in_=False)

        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats(fromFormat, toFormat)

        # Retrieve 'from' and 'to' option flags.
        fromFlags = request.form['from_flags']
        toFlags = request.form['to_flags']

        for char in fromFlags :
            obConversion.AddOption(char, obConversion.INOPTIONS)

        for char in toFlags :
            obConversion.AddOption(char, obConversion.OUTOPTIONS)

        mol = openbabel.OBMol()
        obConversion.ReadFile(mol, 'static/uploads/' + f.filename)
        obConversion.WriteFile(mol, 'static/downloads/' + fname + '.' + toFormat)

        inSize = os.path.getsize('static/uploads/' + f.filename)
        outSize = os.path.getsize('static/downloads/' + fname + '.' + toFormat)
        os.remove('static/uploads/' + f.filename)

        fromFull = request.form['from_full']
        toFull = request.form['to_full']
        quality = request.form['success']
        converter = 'Open Babel'             # $$$$$$$$$$ TODO: Replace hard coding when more than one converter option. $$$$$$$$$$
        out,err = stdouterrOB.reset()

        if err.find('Error') > -1 :
            error_log(fromFull, toFull, converter, fname, err)
            stdouterrOB.done()

            # return http status code 405
            abort(405)
        else :
            log(fromFull, toFull, converter, fname, quality, out, err)

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
        values = [get_date_time(), fromFull, toFull, converter, fname, inSize, outSize, quality]

        logConversion(query, values)
        return 'okay'
    else :
        # return http status code 405
        abort(405)

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
def log(fromFormat, toFormat, converter, fname, quality, out, err) :
    message = create_message(fname, fromFormat, toFormat, converter) + \
              'Quality:   ' + quality + '\n' \
              'Success:   Assuming that the data provided was of the correct format, the conversion\n' \
              '           was successful (to the best of our knowledge) subject to any warnings below.\n' + \
              out + '\n' + err + '\n'

    f = open('static/downloads/' + fname + '.log.txt', 'w')
    f.write(message)
    f.close()

# Write conversion error information to server-side log file.
def error_log(fromFormat, toFormat, converter, fname, err) :
    message = create_message(fname, fromFormat, toFormat, converter) + err + '\n'

    f = open('error_log.txt', 'a')
    f.write(message)
    f.close()

# Create message for log files.
def create_message(fname, fromFormat, toFormat, converter) :
    return 'Date:      ' + get_date() + '\n' \
           'Time:      ' + get_time() + '\n' \
           'File name: ' + fname + '\n' \
           'From:      ' + fromFormat + '\n' \
           'To:        ' + toFormat + '\n' \
           'Converter: ' + converter + '\n'

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
