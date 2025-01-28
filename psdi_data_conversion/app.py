"""app.py

Version 1.0, 8th November 2024

This script acts as a server for the PSDI Data Conversion Service website.
"""

import hashlib
import os
import json
from datetime import datetime
import sys
from flask import Flask, request, render_template, abort, Response

from psdi_data_conversion import log_utility
from psdi_data_conversion import constants as const
from psdi_data_conversion.converter import run_converter

# Create a token by hashing the current date and time.
dt = str(datetime.now())
token = hashlib.md5(dt.encode('utf8')).hexdigest()

# Check the authorisation envvar to see if we're checking auth
auth_ev = os.environ.get(const.AUTH_ENVVAR)
check_auth = (auth_ev is not None) and (auth_ev.lower() == "true")

# Get the logging style from the envvar for it
ev_logging = os.environ.get(const.LOG_MODE_ENVVAR)
if ev_logging is None:
    log_mode = const.LOG_DEFAULT
else:
    ev_logging = ev_logging.lower()
    if ev_logging not in const.L_ALLOWED_LOG_MODES:
        print(f"ERROR: Unrecognised logging option: {ev_logging}. Allowed options are: {const.L_ALLOWED_LOG_MODES}",
              file=sys.stderr)
        exit(1)
    log_mode = ev_logging

app = Flask(__name__)


@app.route('/')
def website():
    """Return the web page along with the token
    """
    # Get the maximum allowed size from the envvar for it
    ev_max_file_size = os.environ.get(const.MAX_FILESIZE_ENVVAR)
    if ev_max_file_size is not None:
        max_file_size = float(ev_max_file_size)*const.MEGABYTE
    else:
        max_file_size = const.DEFAULT_MAX_FILE_SIZE

    data = [{'token': token,
             'max_file_size': max_file_size}]
    return render_template("index.htm", data=data)


@app.route('/convert/', methods=['POST'])
def convert():
    """Convert file to a different format and save to folder 'downloads'. Delete original file. Note that downloading is
    achieved in format.js
    """
    if (not check_auth) or (request.form['token'] == token and token != ''):
        return run_converter(name=request.form['converter'],
                             filename=request.files[const.FILE_TO_UPLOAD_KEY].filename,
                             form=request.form,
                             log_mode=log_mode,
                             abort_callback=abort)
    else:
        # return http status code 405
        abort(405)


@app.route('/conv/', methods=['POST'])
def conv():
    """Convert file (cURL)
    """
    return run_converter(name=request.form['converter'],
                         files=request.files[const.FILE_KEY].filename,
                         form=request.form,
                         log_mode=log_mode,
                         abort_callback=abort)


@app.route('/feedback/', methods=['POST'])
def feedback():
    """Take feedback data from the web app and log it
    """

    try:

        entry = {
            "datetime": log_utility.get_date_time(),
        }

        report = json.loads(request.form['data'])

        for key in ["type", "missing", "reason", "from", "to"]:
            if key in report:
                entry[key] = str(report[key])

        log_utility.append_to_log_file("feedback", entry)

        return Response(status=201)

    except Exception:

        return Response(status=400)


@app.route('/delete/', methods=['POST'])
def delete():
    """Delete files in folder 'downloads'
    """
    os.remove(f"{const.DEFAULT_DOWNLOAD_DIR}/{request.form['filename']}")
    os.remove(f"{const.DEFAULT_DOWNLOAD_DIR}/{request.form['logname']}")

    return 'okay'


@app.route('/del/', methods=['POST'])
def delete_file():
    """Delete file (cURL)
    """
    os.remove(request.form['filepath'])
    return 'Server-side file ' + request.form['filepath'] + ' deleted\n'


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
    if check_auth and request.args['token'] == token and token != '':
        message = '[' + log_utility.get_date_time() + '] ' + request.args['data'] + '\n'

        with open("user_responses", "a") as f:
            f.write(message)

        return 'okay'
    else:
        # return http status code 405
        abort(405)
