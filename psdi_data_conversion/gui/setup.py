"""
# setup.py

This module handles setting up the Flask app
"""


import os
from collections.abc import Callable
from functools import wraps
from typing import Any

import werkzeug
from flask import Flask, cli

import psdi_data_conversion
from psdi_data_conversion import constants as const
from psdi_data_conversion.gui.accessibility import init_accessibility
from psdi_data_conversion.gui.env import get_env
from psdi_data_conversion.gui.get import init_get
from psdi_data_conversion.gui.post import init_post

_app: Flask | None = None


def _patch_flask_warning():
    """Monkey-patch Flask to disable the warnings that would otherwise appear for this so they don't confuse the user
    """

    def suppress_warning(func: Callable[..., Any]) -> Callable[..., Any]:
        @wraps(func)
        def wrapper(*args, **kwargs) -> Any:
            if args and isinstance(args[0], str) and args[0].startswith('WARNING: This is a development server.'):
                return ''
            return func(*args, **kwargs)
        return wrapper

    werkzeug.serving._ansi_style = suppress_warning(werkzeug.serving._ansi_style)
    cli.show_server_banner = lambda *_: None


def _init_app():
    """Create and return the Flask app with appropriate settings"""

    # Suppress Flask's warning, since we're using the dev server as a GUI
    _patch_flask_warning()

    app = Flask(const.APP_NAME)

    # Set the file upload limit based on env vars
    limit_upload_size(app)

    # Connect the app to the various pages and methods of the website
    init_get(app)
    init_post(app)
    init_accessibility(app)

    return app


def limit_upload_size(app: Flask | None = None):
    """Impose a limit on the maximum file that can be uploaded before Flask will raise an error"""

    if app is None:
        app = get_app()

    env = get_env()

    # Determine the largest possible file size that can be uploaded, keeping in mind that 0 indicates unlimited
    larger_max_file_size = env.max_file_size
    if (env.max_file_size > 0) and (env.max_file_size_ob > env.max_file_size):
        larger_max_file_size = env.max_file_size_ob

    if larger_max_file_size > 0:
        app.config['MAX_CONTENT_LENGTH'] = larger_max_file_size
    else:
        app.config['MAX_CONTENT_LENGTH'] = None


def get_app() -> Flask:
    """Get a reference to the global `Flask` app, creating it if necessary.
    """
    global _app
    if not _app:
        _app = _init_app()
    return _app


def start_app():
    """Start the Flask app - this requires being run from the base directory of the project, so this changes the
    current directory to there. Anything else which changes it while the app is running may interfere with its proper
    execution.
    """

    old_cwd = os.getcwd()

    try:
        os.chdir(os.path.join(psdi_data_conversion.__path__[0], ".."))
        get_app().run(debug=get_env().debug_mode)
    finally:
        # Return to the previous directory after running the app
        os.chdir(old_cwd)
