
#   app.py
#   Version 1.0, 17th April 2024

#   This script acts as a server for the Chemistry File Format Conversion Database website.

import hashlib, os, glob, psycopg2
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
            #cursor.execute("SELECT * FROM Converters")
            cursor.execute(request.form['data'])
            results = cursor.fetchall()

        # Close connection to database
        if db_conn:
            db_conn.close()
 
        # Construct a string from the array and return it
        ansArray = []

        for row in results :
            line = ''

            for el in row :
                line += "£" + str(el)

            ansArray.append(line)

        return '$'.join(ansArray)
    else :
        # return http status code 405
        abort(405)

# Convert file to a different format and save to folder 'downloads'. Delete original file.
# Note that downloading is achieved in format.js
@app.route('/convert/', methods=['POST'])
def convert() :
    if request.form['token'] == token and token != '' :
        # Delete any files remaining in folder 'downloads'.
        fileList = glob.glob('static/downloads/*.*')
        for file in fileList :
            os.remove(file)

        f = request.files['fileToUpload']
        f.save('static/uploads/' + f.filename)
        fname = f.filename.split(".")[0]  # E.g. ethane.mol --> ethane

        # Retrieve 'from' and 'to' file formats.
        fromFormat = request.form['from']
        toFormat = request.form['to']

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

        os.remove('static/uploads/' + f.filename)

        return 'okay'
    else :
        # return http status code 405
        abort(405)

# Check that the incoming token matches the one sent to the user (should mostly prevent spambots).
# Write date- and time-stamped user input to server-side file 'user_responses'.
@app.route('/data/', methods=['GET'])
def data() :
    if request.args['token'] == token and token != '' :
        today = datetime.today()
        message = '[' + str(today.year) + '-' + format(today.month) + '-' + format(today.day) + ' ' + format(today.hour) + \
                  ':' + format(today.minute) + ':' + format(today.second) + '] ' + request.args['data'] + '\n'

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
