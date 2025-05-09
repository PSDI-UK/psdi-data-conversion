"""
# get.py

This module defines the various webpages (the "GET" methods) provided by the website, connecting them to relevant
functions to return rendered templates.
"""


from flask import Flask, abort, render_template, request

from psdi_data_conversion import log_utility
from psdi_data_conversion.gui.env import get_env, get_env_kwargs


def index():
    """Return the web page along with relevant data
    """
    return render_template("index.htm",
                           **get_env_kwargs())


def accessibility():
    """Return the accessibility page
    """
    return render_template("accessibility.htm",
                           **get_env_kwargs())


def documentation():
    """Return the documentation page
    """
    return render_template("documentation.htm",
                           **get_env_kwargs())


def data():
    """Check that the incoming token matches the one sent to the user (should mostly prevent spambots). Write date- and
    time-stamped user input to server-side file 'user_responses'.

    $$$$$$$$$$ Retained in case direct logging is required in the future. $$$$$$$$$$

    Returns
    -------
    str
        Output status - 'okay' if exited successfuly
    """
    env = get_env()
    if env.service_mode and request.args['token'] == env.token and env.token != '':
        message = '[' + log_utility.get_date_time() + '] ' + request.args['data'] + '\n'

        with open("user_responses", "a") as f:
            f.write(message)

        return 'okay'
    else:
        # return http status code 405
        abort(405)


def init_get(app: Flask):
    """Connect the provided Flask app to each of the pages on the site
    """

    app.route('/')(index)
    app.route('/index.htm')(index)

    app.route('/accessibility.htm')(accessibility)

    app.route('/documentation.htm')(documentation)

    app.route('/data/')(data)
