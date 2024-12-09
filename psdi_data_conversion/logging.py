"""@file psdi-data-conversion/psdi_data_conversion/logging.py

Created 2024-12-09 by Bryan Gillis.

Functions and classes related to logging
"""

import logging

# Global file to log any errors that occur
ERROR_LOG_FILENAME = "error_log.txt"
GLOBAL_ERROR_LOG = f"./{ERROR_LOG_FILENAME}"

_local_log_file: str | None = None

logging.basicConfig(filename=GLOBAL_ERROR_LOG, level=logging.WARN)


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

        # Set up the local logger unless the local logging filename is `None`
        self._local_logger: logging.Logger | None = None
        self._local_handler: logging.FileHandler | None = None
        if _local_log_file is not None:
            old_logger_class = logging.getLoggerClass()
            try:
                logging.setLoggerClass(logging.Logger)
                self._local_logger = logging.getLogger(name)
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

    def log(self, level, msg, *args, **kwargs):
        """Override the `log` method to log with both the global and local loggers
        """
        super().log(level, msg, *args, **kwargs)
        if self._local_logger is not None:
            self._local_logger.log(level, msg, *args, **kwargs)

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


def getLocalLoggerFilename():
    """Get the filename which will be used for the local logs for any newly-instantiated `DataconversionLogger` objects.
    A value of `None` will indicate no local logging.

    Returns
    -------
    str | None
    """

    return _local_log_file


def getLogger(name, local_log_file=None) -> DataConversionLogger:
    """Gets the `DataConversionLogger` with the provided name if one exists, and creates one if not

    Parameters
    ----------
    name : str
        The desired logging channel for this logger. Should be a period-separated string such as "input.files" etc.
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
