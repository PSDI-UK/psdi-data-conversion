"""@file psdi-data-conversion/psdi_data_conversion/logging.py

Created 2024-12-09 by Bryan Gillis.

Functions and classes related to logging
"""

from datetime import datetime
import json
import logging
import os
import sys

LOG_FORMAT = r'%(asctime)s - %(levelname)s - %(message)s'
TIMESTAMP_FORMAT = r"%Y-%m-%d %H:%M:%S"

# Settings for global logger
GLOBAL_LOG_FILENAME = "./error_log.txt"
GLOBAL_LOGGER_LEVEL = logging.ERROR

# Settings for local logger
NAME = "data-conversion"
DEFAULT_LOCAL_LOGGER_LEVEL = logging.INFO

# Set up the global logger when this module is first imported
global_handler = logging.FileHandler(GLOBAL_LOG_FILENAME)
global_handler.setLevel(GLOBAL_LOGGER_LEVEL)


def setUpDataConversionLogger(name=NAME,
                              local_log_file=None,
                              local_logger_level=DEFAULT_LOCAL_LOGGER_LEVEL,
                              local_logger_raw_output=False,
                              extra_loggers=None,
                              stdout_output_level=None):
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
    extra_loggers : Iterable[Tuple[str, int, bool]]
        A list of one or more tuples of the format (`filename`, `level`, `raw_output`) specifying these options
        (defined the same as for the local logger) for one or more additional logging channels.
    stdout_output_level : int | None
        The logging level (using one of the levels defined in the base Python `logging` module) at and above which to
        log output to stdout. If None (default), nothing will be sent to stdout

    Returns
    -------
    Logger
    """

    # Get a logger using the inherited method before setting up any file handling for it
    logger = logging.getLogger(name)

    if extra_loggers is None:
        extra_loggers = []

    # Set up filehandlers for the global and local logging
    for (filename, level, raw_output) in ((GLOBAL_LOG_FILENAME, GLOBAL_LOGGER_LEVEL, False),
                                          (local_log_file, local_logger_level, local_logger_raw_output),
                                          *extra_loggers):
        _add_filehandler_to_logger(logger, filename, level, raw_output)

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

        stream_handler.setFormatter(logging.Formatter(LOG_FORMAT, datefmt=TIMESTAMP_FORMAT))

    return logger


def _add_filehandler_to_logger(logger, filename, level, raw_output):
    """Private function to add a file handler to a logger only if the logger doesn't already have a handler for that
    file, and set the logging level for the handler
    """
    # Skip if filename is None
    if filename is None:
        return

    file_handler = logging.FileHandler(filename)

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
        os.makedirs(os.path.split(filename)[0], exist_ok=True)

        file_handler = logging.FileHandler(filename)

        logger.addHandler(file_handler)

    # Set or update the logging level and formatter for the handler
    if level is not None:
        file_handler.setLevel(level)
        if level < logger.level or logger.level == logging.NOTSET:
            logger.setLevel(level)
    if not raw_output:
        file_handler.setFormatter(logging.Formatter(LOG_FORMAT, datefmt=TIMESTAMP_FORMAT))

    return


def getDataConversionLogger(name=NAME):
    """A specialisation of getting a logger with `logging.getLogger` which uses a default name, to provide a bulwark
    against using the root logger and potentially breaking something with Flask.
    """

    # Get a logger using the inherited method before setting up any file handling for it
    return logging.getLogger(name)


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


def log(from_format, to_format, converter, fname, calc_type, option, from_flags, to_flags, read_flags_args,
        write_flags_args, quality, out, err):
    """Write Open Babel conversion information to server-side file, ready for downloading to user

    Parameters
    ----------
    from_format : _type_
        _description_
    to_format : _type_
        _description_
    converter : _type_
        _description_
    fname : _type_
        _description_
    calc_type : _type_
        _description_
    option : _type_
        _description_
    from_flags : _type_
        _description_
    to_flags : _type_
        _description_
    read_flags_args : _type_
        _description_
    write_flags_args : _type_
        _description_
    quality : _type_
        _description_
    out : _type_
        _description_
    err : _type_
        _description_
    """

    message = (create_message(fname, from_format, to_format, converter, calc_type, option, from_flags, to_flags,
                              read_flags_args, write_flags_args) +
               'Quality:           ' + quality + '\n'
               'Success:           Assuming that the data provided was of the correct format, the conversion\n'
               '                   was successful (to the best of our knowledge) subject to any warnings below.\n' +
               out + '\n' + err)

    getDataConversionLogger("output").info(message)


def log_ato(from_format, to_format, converter, fname, quality, out, err, logger):
    """Write Atomsk conversion information to server-side file, ready for downloading to user

    Parameters
    ----------
    from_format : _type_
        _description_
    to_format : _type_
        _description_
    converter : _type_
        _description_
    fname : _type_
        _description_
    quality : _type_
        _description_
    out : _type_
        _description_
    err : _type_
        _description_
    logger : Logger
        The logger to use to log the message
    """

    message = (create_message_start(fname, from_format, to_format, converter) +
               'Quality:           ' + quality + '\n'
               'Success:           Assuming that the data provided was of the correct format, the conversion\n'
               '                   was successful (to the best of our knowledge) subject to any warnings below.\n' +
               out + '\n' + err)

    getDataConversionLogger("output").info(message)


def log_error(from_format, to_format, converter, fname, calc_type, option, from_flags, to_flags, read_flags_args,
              write_flags_args, err):
    """Write Open Babel conversion error information to server-side log file

    Parameters
    ----------
    from_format : _type_
        _description_
    to_format : _type_
        _description_
    converter : _type_
        _description_
    fname : _type_
        _description_
    calc_type : _type_
        _description_
    option : _type_
        _description_
    from_flags : _type_
        _description_
    to_flags : _type_
        _description_
    read_flags_args : _type_
        _description_
    write_flags_args : _type_
        _description_
    err : _type_
        _description_
    """
    message = create_message(fname, from_format, to_format, converter, calc_type, option,
                             from_flags, to_flags, read_flags_args, write_flags_args) + err
    getDataConversionLogger().error(message)


def log_error_ato(from_format, to_format, converter, fname, err):
    """Write Atomsk conversion error information to server-side log file

    Parameters
    ----------
    from_format : _type_
        _description_
    to_format : _type_
        _description_
    converter : _type_
        _description_
    fname : _type_
        _description_
    err : _type_
        _description_
    """
    message = create_message(fname, from_format, to_format, converter) + err
    getDataConversionLogger().error(message)


def create_message(fname, from_format, to_format, converter, calc_type, option, from_flags, to_flags, read_flags_args,
                   write_flags_args):
    """Create message for log files

    Parameters
    ----------
    fname : _type_
        _description_
    from_format : _type_
        _description_
    to_format : _type_
        _description_
    converter : _type_
        _description_
    calc_type : _type_
        _description_
    option : _type_
        _description_
    from_flags : _type_
        _description_
    to_flags : _type_
        _description_
    read_flags_args : _type_
        _description_
    write_flags_args : _type_
        _description_

    Returns
    -------
    message : str
        The message for log files
    """
    message = ''

    if calc_type == 'neither':
        message = 'Coord. gen.:       none\n'
    else:
        message += 'Coord. gen.:       ' + calc_type + '\n'

    message += 'Coord. option:     ' + option + '\n'

    if from_flags == '':
        message += 'Read options:      none\n'
    else:
        message += 'Read options:      ' + from_flags + '\n'

    if to_flags == '':
        message += 'Write options:     none\n'
    else:
        message += 'Write options:     ' + to_flags + '\n'

    if len(read_flags_args) == 0:
        message += 'Read opts + args:  none\n'
    else:
        heading_added = False

        for pair in read_flags_args:
            if not heading_added:
                message += 'Read opts + args:  ' + pair + '\n'
                heading_added = True
            else:
                message += '                   ' + pair + '\n'

    if len(write_flags_args) == 0:
        message += 'Write opts + args: none\n'
    else:
        heading_added = False

        for pair in write_flags_args:
            if not heading_added:
                message += 'Write opts + args: ' + pair + '\n'
                heading_added = True
            else:
                message += '                   ' + pair + '\n'

    return create_message_start(fname, from_format, to_format, converter) + message


def create_message_start(fname, from_format, to_format, converter):
    """Create beginning of message for log files

    Parameters
    ----------
    fname : _type_
        _description_
    from_format : _type_
        _description_
    to_format : _type_
        _description_
    converter : _type_
        _description_

    Returns
    -------
    str
        The beginning of a message for log files, containing generic information about what was trying to be done
    """
    return ('Date:              ' + get_date() + '\n'
            'Time:              ' + get_time() + '\n'
            'File name:         ' + fname + '\n'
            'From:              ' + from_format + '\n'
            'To:                ' + to_format + '\n'
            'Converter:         ' + converter + '\n')


def append_to_log_file(log_name, data):
    """Append data to a log file

    Parameters
    ----------
    log_name : _type_
        _description_
    data : _type_
        _description_
    """

    if (os.environ.get('ENABLE_DCS_LOG') is not None):
        print(json.dumps(data))


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
