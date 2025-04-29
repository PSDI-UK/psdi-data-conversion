"""@file psdi_data_conversion/constants.py

Created 2025-01-23 by Bryan Gillis.

Miscellaneous constant values used within this project.

These values are stored as constants rather than hardcoded literals for various reasons, including:
- Better assurance of consistency that the same value is used every time
- If a value needs to be changed, this only needs to be done at one location
- Compatibility with IDE features - a constant can be checked for validity by an IDE, while e.g. a string key for a dict
  can't, allowing more errors to be caught and fixed by the IDE rather than at runtime
- The use of a constant may improve readability - e.g. `MEGABYTE = 1024*1024; max_file_size = 1*MEGABYTE` is more
  readable than `max_file_size = 1*1024*1024`, and so doesn't require a comment like the latter would

There are some known drawbacks to this approach which need to be considered though:
- Constants may obscure readability - it may be quite relevant to the reader exactly what a constant represents, which
  is obscured until they (at minimum) mouse over it
- More code is necessary to use a constant than a literal (at minimum it needs an extra line to define it, and when
  stored here, it also needs a line to import it or this module)

With these drawbacks in mind, we make the following recommendations for constant use:
- Messages for the user (print/logging messages, exceptions) should by default not be stored as constants. They should
  be made constants if it's necessary to reference their exact text elsewhere (either in the executable code or unit
  tests). In this case, the name of the constant should be descriptive, even if this means a rather long name
- If a value is only used in one file and only likely to ever be used in that file, it can be defined as a constant
  there (or if used only two or three times in quick succession, left as a literal)
- Of course, deviations from this should be made when necessary, such as to avoid circular imports
"""

import logging
import shutil

# Interface
# ---------

# The name of the command-line script
CL_SCRIPT_NAME = "psdi-data-convert"

# Environmental variables
LOG_MODE_EV = "LOG_MODE"
LOG_LEVEL_EV = "LOG_LEVEL"
MAX_FILESIZE_EV = "MAX_FILESIZE"
MAX_FILESIZE_OB_EV = "MAX_FILESIZE_OB"
SERVICE_MODE_EV = "SERVICE_MODE"

# Files and Folders
# -----------------

MEGABYTE = 1024*1024

# Maximum output file size in bytes
DEFAULT_MAX_FILE_SIZE = 0 * MEGABYTE
DEFAULT_MAX_FILE_SIZE_OB = 1 * MEGABYTE

DEFAULT_UPLOAD_DIR = './psdi_data_conversion/static/uploads'
DEFAULT_DOWNLOAD_DIR = './psdi_data_conversion/static/downloads'

# Filename of the database, relative to the base of the python package
DATABASE_FILENAME = "static/data/data.json"

# Archive extensions and formats ('format' here meaning the value expected by shutil's archive functions)

ZIP_EXTENSION = ".zip"
ZIP_FORMAT = "zip"

D_ZIP_FORMATS = {ZIP_EXTENSION: ZIP_FORMAT}

TAR_EXTENSION = ".tar"
TAR_FORMAT = "tar"
GZTAR_EXTENSION = ".tar.gz"
GZTAR_FORMAT = "gztar"
BZTAR_EXTENSION = ".tar.bz"
BZTAR_FORMAT = "bztar"
XZTAR_EXTENSION = ".tar.xz"
XZTAR_FORMAT = "xztar"

D_TAR_FORMATS = {TAR_EXTENSION: TAR_FORMAT,
                 GZTAR_EXTENSION: GZTAR_FORMAT,
                 BZTAR_EXTENSION: BZTAR_FORMAT,
                 XZTAR_EXTENSION: XZTAR_FORMAT}

# A list of specifically the extensions that are combinations of multiple different extensions
L_COMPOUND_EXTENSIONS = [GZTAR_EXTENSION, BZTAR_EXTENSION, XZTAR_EXTENSION]

# Formats which are supported by shutil's built-in archive utility
D_SUPPORTED_ARCHIVE_FORMATS = {**D_ZIP_FORMATS, **D_TAR_FORMATS}

L_UNSUPPORTED_ARCHIVE_EXTENSIONS = [".rar", ".7z"]

L_ALL_ARCHIVE_EXTENSIONS = [*D_SUPPORTED_ARCHIVE_FORMATS.keys(), *L_UNSUPPORTED_ARCHIVE_EXTENSIONS]


# Logging and Formatting
# ----------------------

# Number of character spaces allocated for flags/options

# Get the terminal width so we can prettily print help text - default to 80 chars by 20 lines
TERM_WIDTH, _ = shutil.get_terminal_size((80, 20))

# Log formatting
LOG_FORMAT = r'[%(asctime)s] %(levelname)s [%(name)s.%(funcName)s:%(lineno)d] %(message)s'
TIMESTAMP_FORMAT = r"%Y-%m-%d %H:%M:%S"

# Regex to match date/time format
DATE_RE_RAW = r"\d{4}-[0-1]\d-[0-3]\d"
TIME_RE_RAW = r"[0-2]\d:[0-5]\d:[0-5]\d"
DATETIME_RE_RAW = f"{DATE_RE_RAW} {TIME_RE_RAW}"

# Log mode info and settings
LOG_FULL = "full"
LOG_FULL_FORCE = "full-force"
LOG_SIMPLE = "simple"
LOG_STDOUT = "stdout"
LOG_NONE = "none"

LOG_MODE_DEFAULT = LOG_SIMPLE

L_ALLOWED_LOG_MODES = (LOG_FULL, LOG_FULL_FORCE,  LOG_SIMPLE, LOG_STDOUT, LOG_NONE)

LOG_EXT = ".log"
OUTPUT_LOG_EXT = f"{LOG_EXT}.txt"

# Settings for global logger
GLOBAL_LOG_FILENAME = "./error_log.txt"
GLOBAL_LOGGER_LEVEL = logging.ERROR

# Settings for local logger
LOCAL_LOGGER_NAME = "data-conversion"
DEFAULT_LOCAL_LOGGER_LEVEL = logging.INFO
DEFAULT_LISTING_LOG_FILE = "data-convert-list" + LOG_EXT

# Converters and Related
# ----------------------

# Converter names are determined based on the modules present in the 'converters' package by the 'converter' module
# This module contains constant dicts and lists of registered converters

# Default converter - this must match the name of one of the registered converters
CONVERTER_DEFAULT = 'Open Babel'

# File format properties which are used to judge conversion quality - KEY is the key for it in the database, and LABEL
# is how we want to print it out for the user
QUAL_COMP_KEY = "composition"
QUAL_COMP_LABEL = "Atomic composition is"
QUAL_CONN_KEY = "connections"
QUAL_CONN_LABEL = "Atomic connections are"
QUAL_2D_KEY = "two_dim"
QUAL_2D_LABEL = "2D atomic coordinates are"
QUAL_3D_KEY = "three_dim"
QUAL_3D_LABEL = "2D atomic coordinates are"

D_QUAL_LABELS = {QUAL_COMP_KEY: QUAL_COMP_LABEL,
                 QUAL_CONN_KEY: QUAL_CONN_LABEL,
                 QUAL_2D_KEY: QUAL_2D_LABEL,
                 QUAL_3D_KEY: QUAL_3D_LABEL}

# Notes for conversion quality
QUAL_NOTE_IN_UNKNOWN = ("Potential data extrapolation: {} represented in the output format but its representation in "
                        "the input format is unknown")
QUAL_NOTE_OUT_UNKNOWN = ("Potential data loss: {} represented in the input format, but its representation in the "
                         "output format is unknown")
QUAL_NOTE_BOTH_UNKNOWN = ("Potential data loss or extrapolation: {} potentially not supported in either or both of "
                          "the input and output formats")
QUAL_NOTE_IN_MISSING = "Potential data extrapolation: {} represented in the output format but not the input format"
QUAL_NOTE_OUT_MISSING = "Potential data loss: {} represented in the input format but not the output format"

# Conversion quality strings
QUAL_UNKNOWN = 'unknown'
QUAL_VERYGOOD = 'very good'
QUAL_GOOD = 'good'
QUAL_OKAY = 'okay'
QUAL_POOR = 'poor'
QUAL_VERYPOOR = 'very poor'

# Errors
# ------

# HTTP status codes for various types of errors
STATUS_CODE_BAD_METHOD = 405
STATUS_CODE_SIZE = 413
STATUS_CODE_GENERAL = 500

# Error messages
ERR_CONVERTER_NOT_RECOGNISED = "Converter {} not recognized. Allowed converters are: "
ERR_WRONG_EXTENSIONS = "Input file '{file}' does not have expected extension '{ext}'"
ERR_EMPTY_ARCHIVE = "No files to convert were contained in archive"
ERR_CONVERSION_FAILED = ("File conversion failed for one or more files. Lines from the output log "
                         "{} which indicate possible sources of error: ")
