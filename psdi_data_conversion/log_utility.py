"""@file psdi-data-conversion/psdi_data_conversion/logging.py

Created 2024-12-09 by Bryan Gillis.

Functions and classes related to logging and other messaging for the user
"""

import logging
import os
import re
import sys
from datetime import datetime

from psdi_data_conversion import constants as const

D_LOG_LEVELS = {"notset": logging.NOTSET,
                "debug": logging.DEBUG,
                "info": logging.INFO,
                "warn": logging.WARNING,
                "warning": logging.WARNING,
                "error": logging.ERROR,
                "critical": logging.CRITICAL,
                "fatal": logging.CRITICAL}


def get_log_level_from_str(log_level_str: str | None) -> int:
    """Gets a log level, as one of the literal ints defined in the `logging` module, from the string representation
    of it.
    """

    if not log_level_str:
        return logging.NOTSET
    try:
        return D_LOG_LEVELS[log_level_str.lower()]
    except KeyError:
        raise ValueError(f"Unrecognised logging level: '{log_level_str}'. Allowed levels are (case-insensitive): "
                         f"{list(D_LOG_LEVELS.keys())}")


def set_up_data_conversion_logger(name=const.LOCAL_LOGGER_NAME,
                                  local_log_file=None,
                                  local_logger_level=const.DEFAULT_LOCAL_LOGGER_LEVEL,
                                  local_logger_raw_output=False,
                                  extra_loggers=None,
                                  suppress_global_handler=False,
                                  stdout_output_level=None,
                                  mode="a"):
    """Registers a logger with the provided name and sets it up with the desired options

    Parameters
    ----------
    name : str | None
        The desired logging channel for this logger. Should be a period-separated string such as "input.files" etc.
        By default "data-conversion"
    local_log_file : str | None
        The file to log to for local logs. If None, will not set up local logging
    local_logger_level : int
        The logging level to set up for the local logger, using one of the levels defined in the base Python `logging`
        module, by default `logging.INFO`
    local_logger_raw_output : bool
        If set to True, output to the local logger will be logged with no formatting, exactly as input. Otherwise
        (default) it will include a timestamp and indicate the logging level
    extra_loggers : Iterable[Tuple[str, int, bool, str]]
        A list of one or more tuples of the format (`filename`, `level`, `raw_output`, `mode`) specifying these
        options (defined the same as for the local logger) for one or more additional logging channels.
    suppress_global_handler : bool
        If set to True, will not add the handler which sends all logs to the global log file, default False
    stdout_output_level : int | None
        The logging level (using one of the levels defined in the base Python `logging` module) at and above which to
        log output to stdout. If None (default), nothing will be sent to stdout
    mode : str
        Either "a" for append to existing log or "w" to overwrite existing log, default "a"

    Returns
    -------
    Logger
    """

    # Get a logger using the inherited method before setting up any file handling for it
    logger = logging.Logger(name)

    if extra_loggers is None:
        extra_loggers = []

    # Set up filehandlers for the global and local logging
    for (filename, level,
         raw_output, write_mode) in ((const.GLOBAL_LOG_FILENAME, const.GLOBAL_LOGGER_LEVEL, False, "a"),
                                     (local_log_file, local_logger_level, local_logger_raw_output, mode),
                                     *extra_loggers):
        if level is None or (suppress_global_handler and filename == const.GLOBAL_LOG_FILENAME):
            continue
        _add_filehandler_to_logger(logger, filename, level, raw_output, write_mode)

    # Set up stdout output if desired
    if stdout_output_level is not None:

        stream_handler = logging.StreamHandler(sys.stdout)

        # Check if stdout output is already handled, and update that handler if so
        handler_already_present = False
        for handler in logger.handlers:
            if (isinstance(handler, logging.StreamHandler) and handler.stream == sys.stdout):
                handler_already_present = True
                stream_handler = handler
                break

        if not handler_already_present:
            logger.addHandler(stream_handler)

        stream_handler.setLevel(stdout_output_level)
        if stdout_output_level < logger.level or logger.level == logging.NOTSET:
            logger.setLevel(stdout_output_level)

        stream_handler.setFormatter(logging.Formatter(const.LOG_FORMAT, datefmt=const.TIMESTAMP_FORMAT))

    return logger


def _add_filehandler_to_logger(logger, filename, level, raw_output, mode):
    """Private function to add a file handler to a logger only if the logger doesn't already have a handler for that
    file, and set the logging level for the handler
    """
    # Skip if filename is None
    if filename is None:
        return

    file_handler = logging.FileHandler(filename, mode)

    # Check if the file to log to is already in the logger's filehandlers
    handler_already_present = False
    for handler in logger.handlers:
        if (isinstance(handler, logging.FileHandler) and
                handler.baseFilename == os.path.abspath(filename)):
            handler_already_present = True
            file_handler = handler
            break

    # Add a FileHandler for the file if it's not already present, make sure the path to the log file exists,
    # and set the logging level
    if not handler_already_present:
        filename_loc = os.path.split(filename)[0]
        if filename_loc != "":
            os.makedirs(filename_loc, exist_ok=True)

        file_handler = logging.FileHandler(filename)

        logger.addHandler(file_handler)

    # Set or update the logging level and formatter for the handler
    if level is not None:
        file_handler.setLevel(level)
        if level < logger.level or logger.level == logging.NOTSET:
            logger.setLevel(level)
    if not raw_output:
        file_handler.setFormatter(logging.Formatter(const.LOG_FORMAT, datefmt=const.TIMESTAMP_FORMAT))

    return


def get_date():
    """Retrieve current date as a string

    Returns
    -------
    str
        Current date in the format YYYY-MM-DD
    """
    today = datetime.today()
    return str(today.year) + '-' + format(today.month) + '-' + format(today.day)


def get_time():
    """Retrieve current time as a string

    Returns
    -------
    str
        Current time in the format HH:MM:SS
    """
    today = datetime.today()
    return format(today.hour) + ':' + format(today.minute) + ':' + format(today.second)


def get_date_time():
    """Retrieve current date and time as a string

    Returns
    -------
    str
        Current date and time in the format YYYY-MM-DD HH:MM:SS
    """
    return get_date() + ' ' + get_time()


def format(time):
    """Ensure that an element of date or time (month, day, hours, minutes or seconds) always has two digits.

    Parameters
    ----------
    time : str or int
        Digit(s) indicating date or month

    Returns
    -------
    str
        2-digit value indicating date or month
    """
    num = str(time)

    if len(num) == 1:
        return '0' + num
    else:
        return num


def string_with_placeholders_matches(test_pattern: str, parent_str: str) -> bool:
    """An advanced version of "`test_pattern` in `parent_str`" for checking if a substring is in a containing string,
    which allows for `test_pattern` to contain placeholders (e.g. a string like "The file name is: {file}".).

    The test here splits `test_pattern` up into segments which exclude the placeholders, and checks that all segments
    are in `parent_string`. As such, it has the potential to return false positives in the segments are all in the
    parent string but split far apart or out of order, so this method should not be used if it is critical that false
    positives be avoided.

    Parameters
    ----------
    test_pattern : str
        The pattern to check if it's contained in `parent_str`
    parent_str : str
        The string to search in for `test_pattern`

    Returns
    -------
    bool
        True if `test_pattern` appears to be in `parent_str` with some placeholders filled in, False otherwise
    """

    l_test_segments = re.split(r"\{.*?\}", test_pattern)
    all_segments_in_parent = all([s in parent_str for s in l_test_segments])

    return all_segments_in_parent
