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
DEFAULT_MAX_FILE_SIZE = 0*MEGABYTE

DEFAULT_UPLOAD_DIR = './psdi_data_conversion/static/uploads'
DEFAULT_DOWNLOAD_DIR = './psdi_data_conversion/static/downloads'

# Filename of the database, relative to the base of the python package
DATABASE_FILENAME = "static/data/data.json"

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

# Key for the label given to the file uploaded in the web interface
FILE_TO_UPLOAD_KEY = 'fileToUpload'

# Keys for top-level and general items in the database
DB_FORMATS_KEY = "formats"
DB_CONVERTERS_KEY = "converters"
DB_CONVERTS_TO_KEY = "converts_to"
DB_ID_KEY = "id"
DB_NAME_KEY = "name"

# Keys for converter general info in the database
DB_DESC_KEY = "description"
DB_URL_KEY = "url"

# Keys for format general info in the database
DB_FORMAT_EXT_KEY = "extension"
DB_FORMAT_NOTE_KEY = "note"
DB_FORMAT_COMP_KEY = "composition"
DB_FORMAT_CONN_KEY = "connections"
DB_FORMAT_2D_KEY = "two_dim"
DB_FORMAT_3D_KEY = "three_dim"

# Keys for converts_to info in the database
DB_CONV_ID_KEY = "converters_id"
DB_IN_ID_KEY = "in_id"
DB_OUT_ID_KEY = "out_id"
DB_SUCCESS_KEY = "degree_of_success"

# Key bases for converter-specific items in the database
DB_IN_FLAGS_KEY_BASE = "flags_in"
DB_OUT_FLAGS_KEY_BASE = "flags_out"
DB_IN_OPTIONS_KEY_BASE = "argflags_in"
DB_OUT_OPTIONS_KEY_BASE = "argflags_out"
DB_IN_FLAGS_FORMATS_KEY_BASE = "format_to_flags_in"
DB_OUT_FLAGS_FORMATS_KEY_BASE = "format_to_flags_out"
DB_IN_OPTIONS_FORMATS_KEY_BASE = "format_to_argflags_in"
DB_OUT_OPTIONS_FORMATS_KEY_BASE = "format_to_argflags_out"

# Errors
# ------

# Status codes for various types of errors
STATUS_CODE_BAD_METHOD = 405
STATUS_CODE_SIZE = 413
STATUS_CODE_GENERAL = 422
