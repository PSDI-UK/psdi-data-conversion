"""app.py

Version 1.0, 8th November 2024

This script acts as a server for the PSDI Data Conversion Service website.
"""

from argparse import ArgumentParser
import hashlib
import os
import json
from datetime import datetime
from subprocess import run
import sys
import traceback
from flask import Flask, request, render_template, abort, Response

import psdi_data_conversion
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

# Get the maximum allowed size from the envvar for it
ev_max_file_size = os.environ.get(const.MAX_FILESIZE_EV)
if ev_max_file_size is not None:
    max_file_size = float(ev_max_file_size)*const.MEGABYTE
else:
    max_file_size = const.DEFAULT_MAX_FILE_SIZE

# And same for the Open Babel maximum file size
ev_max_file_size_ob = os.environ.get(const.MAX_FILESIZE_OB_EV)
if ev_max_file_size_ob is not None:
    max_file_size_ob = float(ev_max_file_size_ob)*const.MEGABYTE
else:
    max_file_size_ob = const.DEFAULT_MAX_FILE_SIZE_OB

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

    data = [{'token': token,
             'max_file_size': max_file_size,
             'max_file_size_ob': max_file_size_ob,
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
                                              strict=(request.form['check_ext'] != "false"),
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

            # If the exception provides a status code, get it
            status_code: int
            if hasattr(e, "status_code"):
                status_code = e.status_code
            else:
                status_code = const.STATUS_CODE_GENERAL

            # If the exception provides a message, report it
            if hasattr(e, "message"):
                msg = f"An unexpected exception was raised by the converter, with error message:\n{e.message}\n"
            else:
                # Failsafe exception message
                msg = ("The following unexpected exception was raised by the converter:\n" +
                       traceback.format_exc()+"\n")
            with open(qualified_output_log, "w") as fo:
                fo.write(msg)
            abort(status_code)

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


def start_app():
    """Start the Flask app - this requires being run from the base directory of the project, so this changes the
    current directory to there. Anything else which changes it while the app is running may interfere with its proper
    execution.
    """

    os.chdir(os.path.join(psdi_data_conversion.__path__[0], ".."))
    app.run()


def main():
    """Standard entry-point function for this script.
    """

    parser = ArgumentParser()

    parser.add_argument("--use-env-vars", action="store_true",
                        help="If set, all other arguments and defaults for this script are ignored, and environmental "
                        "variables and their defaults will instead control execution. These defaults will result in "
                        "the app running in production server mode.")

    parser.add_argument("--max-file-size", type=float, default=const.DEFAULT_MAX_FILE_SIZE,
                        help="The maximum allowed filesize in MB - 0 (default) indicates no maximum")

    parser.add_argument("--max-file-size-ob", type=float, default=const.DEFAULT_MAX_FILE_SIZE_OB,
                        help="The maximum allowed filesize in MB for the Open Babel converter, taking precendence over "
                        "the general maximum file size when Open Babel is used - 0 indicates no maximum. Default 1 MB.")

    parser.add_argument("--service-mode", action="store_true",
                        help="If set, will run as if deploying a service rather than the local GUI")

    parser.add_argument("--dev-mode", action="store_true",
                        help="If set, will expose development elements")

    parser.add_argument("--log-mode", type=str, default=const.LOG_FULL,
                        help="How logs should be stored. Allowed values are: \n"
                        "- 'full' - Multi-file logging, not recommended for the CLI, but allowed for a compatible "
                        "interface with the public web app"
                        "- 'simple' - Logs saved to one file"
                        "- 'stdout' - Output logs and errors only to stdout"
                        "- 'none' - Output only errors to stdout")

    parser.add_argument("--log-level", type=str, default=None,
                        help="The desired level to log at. Allowed values are: 'DEBUG', 'INFO', 'WARNING', 'ERROR, "
                             "'CRITICAL'. Default: 'INFO' for logging to file, 'WARNING' for logging to stdout")

    # Set global variables for settings based on parsed arguments, unless it's set to use env vars
    args = parser.parse_args()

    if not args.use_env_vars:

        global max_file_size
        max_file_size = args.max_file_size*const.MEGABYTE

        global max_file_size_ob
        max_file_size_ob = args.max_file_size_ob*const.MEGABYTE

        global service_mode
        service_mode = args.service_mode

        global production_mode
        production_mode = not args.dev_mode

        global log_mode
        log_mode = args.log_mode

        global log_level
        log_level = args.log_level

    start_app()


if __name__ == "__main__":

    main()
