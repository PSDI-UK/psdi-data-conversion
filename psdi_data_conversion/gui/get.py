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


def convert_ob():
    """Return the Open Babel convert page
    """
    return render_template("convert_ob.htm",
                           **get_env_kwargs())


def convert_ato():
    """Return the Atomsk convert page
    """
    return render_template("convert_ato.htm",
                           **get_env_kwargs())


def convert_c2x():
    """Return the c2x convert page
    """
    return render_template("convert_c2x.htm",
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


def download():
    """Return the download page
    """
    return render_template("download.htm",
                           **get_env_kwargs())


def feedback():
    """Return the feedback page
    """
    return render_template("feedback.htm",
                           **get_env_kwargs())


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

    app.route('/convert_ob.htm')(convert_ob)
    app.route('/convert_ato.htm')(convert_ato)
    app.route('/convert_c2x.htm')(convert_c2x)
    app.route('/database/')(database)
    app.route('/documentation.htm')(documentation)
    app.route('/download.htm')(download)
    app.route('/feedback.htm')(feedback)
    app.route('/report.htm')(report)
