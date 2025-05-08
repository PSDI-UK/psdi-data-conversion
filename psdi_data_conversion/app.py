"""app.py

Version 1.0, 8th November 2024

This script acts as a server for the PSDI Data Conversion Service website.
"""

import json
import os
import sys
from argparse import ArgumentParser
from collections.abc import Callable
from datetime import datetime
from functools import wraps
from hashlib import md5
from multiprocessing import Lock
from traceback import format_exc
from typing import Any

import werkzeug.serving
from flask import Flask, Response, abort, cli, render_template, request
from werkzeug.utils import secure_filename

import psdi_data_conversion
from psdi_data_conversion import constants as const
from psdi_data_conversion import log_utility
from psdi_data_conversion.converter import run_converter
from psdi_data_conversion.database import get_format_info
from psdi_data_conversion.file_io import split_archive_ext
from psdi_data_conversion.main import print_wrap
from psdi_data_conversion.website.env import get_env

# Key for the label given to the file uploaded in the web interface
FILE_TO_UPLOAD_KEY = 'fileToUpload'

# A lock to prevent multiple threads from logging at the same time
logLock = Lock()

# Create a token by hashing the current date and time.
dt = str(datetime.now())
token = md5(dt.encode('utf8')).hexdigest()


def suppress_warning(func: Callable[..., Any]) -> Callable[..., Any]:
    """Since we're using the development server as the user GUI, we monkey-patch Flask to disable the warnings that
    would otherwise appear for this so they don't confuse the user
    """
    @wraps(func)
    def wrapper(*args, **kwargs) -> Any:
        if args and isinstance(args[0], str) and args[0].startswith('WARNING: This is a development server.'):
            return ''
        return func(*args, **kwargs)
    return wrapper


werkzeug.serving._ansi_style = suppress_warning(werkzeug.serving._ansi_style)
cli.show_server_banner = lambda *_: None

app = Flask(__name__)


def limit_upload_size():
    """Impose a limit on the maximum file that can be uploaded before Flask will raise an error"""

    env = get_env()

    # Determine the largest possible file size that can be uploaded, keeping in mind that 0 indicates unlimited
    larger_max_file_size = env.max_file_size
    if (env.max_file_size > 0) and (env.max_file_size_ob > env.max_file_size):
        larger_max_file_size = env.max_file_size_ob

    if larger_max_file_size > 0:
        app.config['MAX_CONTENT_LENGTH'] = larger_max_file_size
    else:
        app.config['MAX_CONTENT_LENGTH'] = None


# Set the upload limit based on env vars to start with
limit_upload_size()


@app.route('/')
@app.route('/index.htm')
def website():
    """Return the web page along with relevant data
    """
    return render_template("index.htm",
                           token=token,
                           **get_env().d_kwargs)


@app.route('/accessibility.htm')
def accessibility():
    """Return the accessibility page
    """
    return render_template("accessibility.htm",
                           **get_env().d_kwargs)


@app.route('/documentation.htm')
def documentation():
    """Return the documentation page
    """
    return render_template("documentation.htm",
                           **get_env().d_kwargs)


@app.route('/convert/', methods=['POST'])
def convert():
    """Convert file to a different format and save to folder 'downloads'. Delete original file. Note that downloading is
    achieved in format.js
    """

    env = get_env()

    # Make sure the upload directory exists
    os.makedirs(const.DEFAULT_UPLOAD_DIR, exist_ok=True)

    file = request.files[FILE_TO_UPLOAD_KEY]
    filename = secure_filename(file.filename)

    qualified_filename = os.path.join(const.DEFAULT_UPLOAD_DIR, filename)
    file.save(qualified_filename)
    qualified_output_log = os.path.join(const.DEFAULT_DOWNLOAD_DIR,
                                        split_archive_ext(filename)[0] + const.OUTPUT_LOG_EXT)

    # Determine the input and output formats
    d_formats = {}
    for format_label in "to", "from":
        name = request.form[format_label]
        full_note = request.form[format_label+"_full"]

        l_possible_formats = get_format_info(name, which="all")

        # If there's only one possible format, use that
        if len(l_possible_formats) == 1:
            d_formats[format_label] = l_possible_formats[0]
            continue

        # Otherwise, find the format with the matching note
        for possible_format in l_possible_formats:
            if possible_format.note in full_note:
                d_formats[format_label] = possible_format
                break
        else:
            print(f"Format '{name}' with full description '{full_note}' could not be found in database.",
                  file=sys.stderr)
            abort(const.STATUS_CODE_GENERAL)

    if (not env.service_mode) or (request.form['token'] == token and token != ''):
        try:
            conversion_output = run_converter(name=request.form['converter'],
                                              filename=qualified_filename,
                                              data=request.form,
                                              to_format=d_formats["to"],
                                              from_format=d_formats["from"],
                                              strict=(request.form['check_ext'] != "false"),
                                              log_mode=env.log_mode,
                                              log_level=env.log_level,
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
                       format_exc()+"\n")
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

        # Write data in JSON format and send to stdout
        logLock.acquire()
        sys.stdout.write(f"{json.dumps(entry) + '\n'}")
        logLock.release()

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
    if get_env().service_mode and request.args['token'] == token and token != '':
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
    app.run(debug=get_env().debug_mode)


def main():
    """Standard entry-point function for this script.
    """

    parser = ArgumentParser()

    parser.add_argument("--use-env-vars", action="store_true",
                        help="If set, all other arguments and defaults for this script are ignored, and environmental "
                        "variables and their defaults will instead control execution. These defaults will result in "
                        "the app running in production server mode.")

    parser.add_argument("--max-file-size", type=float, default=const.DEFAULT_MAX_FILE_SIZE/const.MEGABYTE,
                        help="The maximum allowed filesize in MB - 0 (default) indicates no maximum")

    parser.add_argument("--max-file-size-ob", type=float, default=const.DEFAULT_MAX_FILE_SIZE_OB/const.MEGABYTE,
                        help="The maximum allowed filesize in MB for the Open Babel converter, taking precendence over "
                        "the general maximum file size when Open Babel is used - 0 indicates no maximum. Default 1 MB.")

    parser.add_argument("--service-mode", action="store_true",
                        help="If set, will run as if deploying a service rather than the local GUI")

    parser.add_argument("--dev-mode", action="store_true",
                        help="If set, will expose development elements, such as the SHA of the latest commit")

    parser.add_argument("--debug", action="store_true",
                        help="If set, will run the Flask server in debug mode, which will cause it to automatically "
                        "reload if code changes and show an interactive debugger in the case of errors")

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
        # Use `get_env` to overwrite the values from environmental variables with the values from the command-line
        # arguments
        get_env(args=args)

    # Set the upload limit based on provided arguments now
    limit_upload_size()

    print_wrap("Starting the PSDI Data Conversion GUI. This GUI is run as a webpage, which you can open by "
               "right-clicking the link below to open it in your default browser, or by copy-and-pasting it into your "
               "browser of choice.")

    start_app()


if __name__ == "__main__":

    main()
