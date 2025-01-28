"""@file psdi_data_conversion/constants.py

Created 2025-01-23 by Bryan Gillis.

Miscellaneous constant values used within this project
"""

import logging

# Interface
# ---------

# Environmental variables
AUTH_ENVVAR = "AUTH"
LOG_MODE_ENVVAR = "LOGGING"
MAX_FILESIZE_ENVVAR = "MAX_FILESIZE"

# Allowed and default options for command-line arguments
L_ALLOWED_COORD_GENS = ["Gen2D", "Gen3D", "neither"]
DEFAULT_COORD_GEN = "neither"
L_ALLOWED_COORD_GEN_QUALS = ["fastest", "fast", "medium", "better", "best"]
DEFAULT_COORD_GEN_QUAL = "medium"

# Files and Folders
# -----------------

# Maximum output file size in bytes
MEGABYTE = 1024*1024
DEFAULT_MAX_FILE_SIZE = 1*MEGABYTE

DEFAULT_UPLOAD_DIR = './psdi_data_conversion/static/uploads'
DEFAULT_DOWNLOAD_DIR = './psdi_data_conversion/static/downloads'

# Logging
# -------

# Log formatting
LOG_FORMAT = r'[%(asctime)s] %(levelname)s [%(name)s.%(funcName)s:%(lineno)d] %(message)s'
TIMESTAMP_FORMAT = r"%Y-%m-%d %H:%M:%S"

DATE_RE_RAW = r"\d{4}-[0-1]\d-[0-3]\d"
TIME_RE_RAW = r"[0-2]\d:[0-5]\d:[0-5]\d"
DATETIME_RE_RAW = f"{DATE_RE_RAW} {TIME_RE_RAW}"

# Settings for global logger
GLOBAL_LOG_FILENAME = "./error_log.txt"
GLOBAL_LOGGER_LEVEL = logging.ERROR

# Log mode info and settings
LOG_FULL = "full"
LOG_SIMPLE = "simple"
LOG_STDOUT = "stdout"
LOG_NONE = "none"

LOG_DEFAULT = LOG_SIMPLE

L_ALLOWED_LOG_MODES = (LOG_FULL, LOG_SIMPLE, LOG_STDOUT, LOG_NONE)

LOG_EXT = ".log"
LOCAL_LOG_EXT = LOG_EXT
OUTPUT_LOG_EXT = f"{LOG_EXT}.txt"

# Settings for local logger
LOCAL_LOGGER_NAME = "data-conversion"
DEFAULT_LOCAL_LOGGER_LEVEL = logging.INFO
DEFAULT_LISTING_LOG_FILE = "data-convert-list" + LOG_EXT

# Converters
# ----------

# Converter names are determined based on the modules present in the 'converters' package by the 'converter' module
# This module contains constant dicts and lists of registered converters

# Default converter
CONVERTER_DEFAULT = 'Open Babel'

# Keys
# ----

# Keys used commonly by dicts
FILE_TO_UPLOAD_KEY = 'fileToUpload'

# Errors
# ------

# Status codes for various types of errors
STATUS_CODE_BAD_METHOD = 405
STATUS_CODE_SIZE = 413
STATUS_CODE_GENERAL = 422
