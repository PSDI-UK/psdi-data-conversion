"""app.py

Version 1.0, 8th November 2024

This script acts as a server for the PSDI Data Conversion Service website.
"""

import hashlib
import os
import json
from datetime import datetime
from subprocess import run
import sys
import traceback
from flask import Flask, request, render_template, abort, Response

from psdi_data_conversion import log_utility
from psdi_data_conversion import constants as const
from psdi_data_conversion.converter import run_converter
from psdi_data_conversion.file_io import split_archive_ext

# Env var for the SHA of the latest commit
SHA_EV = "SHA"

# Env var for whether this is running in service mode or locally
SERVICE_MODE_EV = "SERVICE_MODE"

# Env var for whether this is a production release or development
PRODUCTION_EV = "PRODUCTION_MODE"

# Key for the label given to the file uploaded in the web interface
FILE_TO_UPLOAD_KEY = 'fileToUpload'

# Create a token by hashing the current date and time.
dt = str(datetime.now())
token = hashlib.md5(dt.encode('utf8')).hexdigest()

# Get the service and production modes from their envvars
service_mode_ev = os.environ.get(SERVICE_MODE_EV)
service_mode = (service_mode_ev is not None) and (service_mode_ev.lower() == "true")
production_mode_ev = os.environ.get(PRODUCTION_EV)
production_mode = (production_mode_ev is not None) and (production_mode_ev.lower() == "true")

# Get the logging mode and level from their envvars
ev_log_mode = os.environ.get(const.LOG_MODE_EV)
if ev_log_mode is None:
    log_mode = const.LOG_MODE_DEFAULT
else:
    ev_log_mode = ev_log_mode.lower()
    if ev_log_mode not in const.L_ALLOWED_LOG_MODES:
        print(f"ERROR: Unrecognised logging option: {ev_log_mode}. Allowed options are: {const.L_ALLOWED_LOG_MODES}",
              file=sys.stderr)
        exit(1)
    log_mode = ev_log_mode

ev_log_level = os.environ.get(const.LOG_LEVEL_EV)
if ev_log_level is None:
    log_level = None
else:
    try:
        log_level = log_utility.get_log_level_from_str(ev_log_level)
    except ValueError as e:
        print(f"ERROR: {str(e)}")
        exit(1)

app = Flask(__name__)


def get_last_sha() -> str:
    """Get the SHA of the last commit
    """

    # First check if the SHA is provided through an environmental variable
    ev_sha = os.environ.get(SHA_EV)
    if ev_sha:
        return ev_sha

    try:
        # This bash command calls `git log` to get info on the last commit, uses `head` to trim it to one line, then
        # uses `gawk` to get just the second word of this line, which is the SHA of this commit
        cmd = "git log -n 1 | head -n 1 | gawk '{print($2)}'"

        out_bytes = run(cmd, shell=True, capture_output=True).stdout
        out_str = str(out_bytes.decode()).strip()

    except Exception:
        print("ERROR: Could not determine SHA of most recent commit. Error was:\n" + traceback.format_exc(),
              file=sys.stderr)
        out_str = "N/A"

    return out_str


@app.route('/')
def website():
    """Return the web page along with the token
    """
    # Get the maximum allowed size from the envvar for it
    ev_max_file_size = os.environ.get(const.MAX_FILESIZE_EV)
    if ev_max_file_size is not None:
        max_file_size = float(ev_max_file_size)*const.MEGABYTE
    else:
        max_file_size = const.DEFAULT_MAX_FILE_SIZE

    data = [{'token': token,
             'max_file_size': max_file_size,
             'service_mode': service_mode,
             'production_mode': production_mode,
             'sha': get_last_sha()}]
    return render_template("index.htm", data=data)


@app.route('/convert/', methods=['POST'])
def convert():
    """Convert file to a different format and save to folder 'downloads'. Delete original file. Note that downloading is
    achieved in format.js
    """

    # Make sure the upload directory exists
    os.makedirs(const.DEFAULT_UPLOAD_DIR, exist_ok=True)

    # Save the file in the upload directory
    file = request.files[FILE_TO_UPLOAD_KEY]
    filename = filename = file.filename

    qualified_filename = os.path.join(const.DEFAULT_UPLOAD_DIR, filename)
    file.save(qualified_filename)
    qualified_output_log = os.path.join(const.DEFAULT_DOWNLOAD_DIR,
                                        split_archive_ext(filename)[0] + const.OUTPUT_LOG_EXT)

    if (not service_mode) or (request.form['token'] == token and token != ''):
        try:
            conversion_output = run_converter(name=request.form['converter'],
                                              filename=qualified_filename,
                                              data=request.form,
                                              to_format=request.form['to'],
                                              from_format=request.form['from'],
                                              strict=request.form['check_ext'],
                                              log_mode=log_mode,
                                              log_level=log_level,
                                              delete_input=True,
                                              abort_callback=abort)
        except Exception as e:

            # Check for anticipated exceptions, and write a simpler message for them
            for err_message in (const.ERR_CONVERSION_FAILED, const.ERR_CONVERTER_NOT_RECOGNISED,
                                const.ERR_EMPTY_ARCHIVE, const.ERR_WRONG_EXTENSIONS):
                if log_utility.string_with_placeholders_matches(err_message, str(e)):
                    with open(qualified_output_log, "w") as fo:
                        fo.write(str(e))
                    abort(const.STATUS_CODE_GENERAL)

            # Failsafe exception message
            msg = "The following unexpected exception was raised by the converter:\n" + traceback.format_exc()+"\n"
            with open(qualified_output_log, "w") as fo:
                fo.write(msg)
            abort(const.STATUS_CODE_GENERAL)

        return repr(conversion_output)
    else:
        # return http status code 405
        abort(405)


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

    realbase = os.path.realpath(const.DEFAULT_DOWNLOAD_DIR)

    realfilename = os.path.realpath(os.path.join(const.DEFAULT_DOWNLOAD_DIR, request.form['filename']))
    reallogname = os.path.realpath(os.path.join(const.DEFAULT_DOWNLOAD_DIR, request.form['logname']))

    if realfilename.startswith(realbase + os.sep) and reallogname.startswith(realbase + os.sep):

        os.remove(realfilename)
        os.remove(reallogname)

        return 'okay'

    else:

        return Response(status=400)


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
    if service_mode and request.args['token'] == token and token != '':
        message = '[' + log_utility.get_date_time() + '] ' + request.args['data'] + '\n'

        with open("user_responses", "a") as f:
            f.write(message)

        return 'okay'
    else:
        # return http status code 405
        abort(405)
