"""
# post.py

This module defines the various web addresses which do something (the "POST" methods) provided by the website,
connecting them to relevant functions.
"""


import json
import os
import sys
from multiprocessing import Lock
from traceback import format_exc

from flask import Flask, Response, abort, request
from werkzeug.utils import secure_filename

from psdi_data_conversion import constants as const
from psdi_data_conversion import log_utility
from psdi_data_conversion.converter import run_converter
from psdi_data_conversion.database import get_format_info
from psdi_data_conversion.file_io import split_archive_ext
from psdi_data_conversion.gui.env import get_env

# Key for the label given to the file uploaded in the web interface
FILE_TO_UPLOAD_KEY = 'fileToUpload'

# A lock to prevent multiple threads from logging at the same time
logLock = Lock()


def convert():
    """Convert file to a different format and save to folder 'downloads'. Delete original file. Note that downloading is
    achieved in format.js
    """

    env = get_env()

    # Make sure the upload directory exists
    os.makedirs(const.DEFAULT_INPUT_DIR, exist_ok=True)

    file = request.files[FILE_TO_UPLOAD_KEY]
    filename = secure_filename(file.filename)

    qualified_filename = os.path.join(const.DEFAULT_INPUT_DIR, filename)
    file.save(qualified_filename)
    qualified_output_log = os.path.join(const.DEFAULT_OUTPUT_DIR,
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

    if (not env.service_mode) or (request.form['token'] == env.token and env.token != ''):
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


def delete():
    """Delete files in folder 'downloads'
    """

    realbase = os.path.realpath(const.DEFAULT_OUTPUT_DIR)

    realfilename = os.path.realpath(os.path.join(const.DEFAULT_OUTPUT_DIR, request.form['filename']))
    reallogname = os.path.realpath(os.path.join(const.DEFAULT_OUTPUT_DIR, request.form['logname']))

    if realfilename.startswith(realbase + os.sep) and reallogname.startswith(realbase + os.sep):

        os.remove(realfilename)
        os.remove(reallogname)

        return 'okay'

    else:

        return Response(status=400)


def init_post(app: Flask):
    """Connect the provided Flask app to each of the post methods
    """

    app.route('/convert/', methods=["POST"])(convert)

    app.route('/feedback/', methods=["POST"])(feedback)

    app.route('/delete/', methods=["POST"])(delete)
