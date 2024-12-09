"""@file psdi-data-conversion/psdi_data_conversion/logging.py

Created 2024-12-09 by Bryan Gillis.

Functions and classes related to logging
"""

from datetime import datetime
import json
import logging
import os

# Global file to log any errors that occur
ERROR_LOG_FILENAME = "error_log.txt"
GLOBAL_ERROR_LOG = f"./{ERROR_LOG_FILENAME}"


class DataConversionLogger(logging.Logger):
    """Logger which logs to both a global file and a local file.

    Note: For compatibility with the parent class, this class uses `camelCase` for function names.
    """

    def __init__(self, name, level=logging.INFO):
        """Creates a logger as a child class of `logging.Logger`

        Parameters
        ----------
        name : str
            The desired logging channel for this logger. Should be a period-separated string such as "input.files" etc.
        level : int, optional
            The desired logging level for the local logger, using one of the constants defined in the `logging` module,
            by default logging.INFO
        """

        # Set up the global logger
        super().__init__(name, logging.NOTSET)
        self._global_handler = logging.FileHandler(ERROR_LOG_FILENAME)
        self.addHandler(self._global_handler)
        super().setLevel(logging.WARN)

        # Set up the local logger unless the local logging filename is `None`
        self._local_logger: logging.Logger | None = None
        self._local_handler: logging.FileHandler | None = None
        if _local_log_file is not None:
            old_logger_class = logging.getLoggerClass()
            try:
                logging.setLoggerClass(logging.Logger)
                self._local_logger = logging.getLogger(f"{name}-local")
            finally:
                logging.setLoggerClass(old_logger_class)
            self._local_handler = logging.FileHandler(_local_log_file)
            self._local_logger.addHandler(self._local_handler)
            self._local_logger.setLevel(level)

    def setLevel(self, level):
        """Override the `setLevel` method to set the same level with both the global and local loggers
        """
        self.setGlobalLevel(level)
        self.setLocalLevel(level)

    def setGlobalLevel(self, level):
        """Set the level for the global logger
        """
        super().setLevel(level)

    def setLocalLevel(self, level):
        """Set the level for the local logger
        """
        if self._local_logger is not None:
            self._local_logger.setLevel(level)

    def debug(self, *args, **kwargs):
        """Override the `debug` method to log with both the global and local loggers
        """
        super().debug(*args, **kwargs)
        if self._local_logger is not None:
            self._local_logger.debug(*args, **kwargs)

    def info(self, *args, **kwargs):
        """Override the `info` method to log with both the global and local loggers
        """
        super().info(*args, **kwargs)
        if self._local_logger is not None:
            self._local_logger.info(*args, **kwargs)

    def warning(self, *args, **kwargs):
        """Override the `warning` method to log with both the global and local loggers
        """
        super().warning(*args, **kwargs)
        if self._local_logger is not None:
            self._local_logger.warning(*args, **kwargs)

    def error(self, *args, **kwargs):
        """Override the `error` method to log with both the global and local loggers
        """
        super().error(*args, **kwargs)
        if self._local_logger is not None:
            self._local_logger.error(*args, **kwargs)

    def critical(self, *args, **kwargs):
        """Override the `critical` method to log with both the global and local loggers
        """
        super().critical(*args, **kwargs)
        if self._local_logger is not None:
            self._local_logger.critical(*args, **kwargs)

    def log(self, *args, **kwargs):
        """Override the `log` method to log with both the global and local loggers
        """
        super().log(*args, **kwargs)
        if self._local_logger is not None:
            self._local_logger.log(*args, **kwargs)

    def getGlobalFilename(self):
        """Get the filename which is used for the global logger

        Returns
        -------
        str
        """
        return self._global_handler.baseFilename

    def getLocalFilename(self):
        """Get the filename which is used for the local logger

        Returns
        -------
        str | None
        """
        if self._local_handler is None:
            return None
        return self._local_handler.baseFilename


def setLocalLoggerFilename(local_log_file):
    """Sets the filename to be used for the logs of the local logger. `None` can be used to indicate no local logging.

    Parameters
    ----------
    local_log_file : str | None
        The file to log to for local logs. `None` can be used to indicate no local logging.
    """
    global _local_log_file
    _local_log_file = local_log_file


_local_log_file: str | None = None


def getLocalLoggerFilename():
    """Get the filename which will be used for the local logs for any newly-instantiated `DataconversionLogger` objects.
    A value of `None` will indicate no local logging.

    Returns
    -------
    str | None
    """
    return _local_log_file


def getLogger(name=None, local_log_file=None) -> DataConversionLogger:
    """Gets the `DataConversionLogger` with the provided name if one exists, and creates one if not

    Parameters
    ----------
    name : str | None
        The desired logging channel for this logger. Should be a period-separated string such as "input.files" etc.
        By default None, indicating use the root logger
    local_log_file : _type_, optional
        The file to log to for local logs. If not provided, will use the filename last set with `setLocalLoggerFilename`

    Returns
    -------
    DataConversionLogger
    """

    old_logger_class = logging.getLoggerClass()
    old_local_filename = getLocalLoggerFilename()
    try:
        logging.setLoggerClass(DataConversionLogger)
        if local_log_file is not None:
            setLocalLoggerFilename(local_log_file)
        logger = logging.getLogger(name)
        return logger
    finally:
        logging.setLoggerClass(old_logger_class)
        setLocalLoggerFilename(old_local_filename)


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
               out + '\n' + err + '\n')

    getLogger("output").info(message)


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
               out + '\n' + err + '\n')

    getLogger("output").info(message)


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
                             from_flags, to_flags, read_flags_args, write_flags_args) + err + '\n'
    getLogger("error").error(message)


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
    message = create_message(fname, from_format, to_format, converter) + err + '\n'
    getLogger("error").error(message)


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
