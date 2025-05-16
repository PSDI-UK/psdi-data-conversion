"""
# get.py

This module defines the various webpages (the "GET" methods) provided by the website, connecting them to relevant
functions to return rendered templates.
"""


from flask import Flask, render_template

from psdi_data_conversion.database import get_database_path
from psdi_data_conversion.gui.env import get_env_kwargs


def index():
    """Return the web page along with relevant data
    """
    return render_template("index.htm",
                           **get_env_kwargs())


def documentation():
    """Return the documentation page
    """
    return render_template("documentation.htm",
                           **get_env_kwargs())


def database():
    """Return the raw database JSON file
    """
    return open(get_database_path(), "r").read()


def report():
    """Return the report page
    """
    return render_template("report.htm",
                           **get_env_kwargs())


def init_get(app: Flask):
    """Connect the provided Flask app to each of the pages on the site
    """

    app.route('/')(index)
    app.route('/index.htm')(index)

    app.route('/documentation.htm')(documentation)

    app.route('/report.htm')(report)

    app.route('/database/')(database)
