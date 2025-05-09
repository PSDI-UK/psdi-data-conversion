"""
# post.py

This module defines the various web addresses which do something (the "POST" methods) provided by the website,
connecting them to relevant functions.
"""


import json

from flask import Flask, make_response, render_template, request

from psdi_data_conversion import constants as const
from psdi_data_conversion.gui.env import get_env_kwargs


def accessibility():
    """Return the accessibility page
    """
    return render_template("accessibility.htm",
                           **get_env_kwargs())


def save_accessibility():
    """Save the user's accessibility settings in a cookie
    """

    resp = make_response("Cookie saved successfully")

    d_settings: dict[str, str] = json.loads(request.form['data'])

    for key, val in d_settings.items():
        resp.set_cookie(key, val, max_age=const.YEAR)

    return resp


def load_accessibility():
    """Load the user's accessibility settings from the cookie
    """
    return json.dumps(request.cookies)


def init_accessibility(app: Flask):
    """Connect the provided Flask app to each of the post methods
    """

    app.route('/accessibility.htm')(accessibility)

    app.route('/save_accessibility/', methods=["POST"])(save_accessibility)
    app.route('/load_accessibility/', methods=["GET"])(load_accessibility)
