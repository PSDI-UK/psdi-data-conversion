"""@file psdi-data-conversion/psdi_data_conversion/logging.py

Created 2024-12-09 by Bryan Gillis.

Functions and classes related to logging
"""

import logging

# Global file to log any errors that occur
ERROR_LOG_FILENAME = "error_log.txt"
GLOBAL_ERROR_LOG = f"./{ERROR_LOG_FILENAME}"
LOCAL_ERROR_LOG = "./static/downloads/local_error_log.txt"

logging.basicConfig(filename=ERROR_LOG_FILENAME)


class DataConversionLogger(logging.Logger):
    """Logger which logs to both a global file and a local file
    """

    def __init__(self, name, level=logging.NOTSET):
        """Creates a logger as a child class of `logging.Logger`

        Parameters
        ----------
        name : str
            The desired logging channel for this logger. Should be a period-separated string such as "input.files" etc.
        local_log_filename : str
            The filename of a file local to this process to use for logging.
        level : int, optional
            The desired logging level, using one of the constants defined in the `logging` module
        """

        # Set up the global logger
        super().__init__(name, level)
        global_handler = logging.FileHandler(ERROR_LOG_FILENAME)
        self.addHandler(global_handler)

        # Set up the local logger

        old_logger_class = logging.getLoggerClass()
        try:
            logging.setLoggerClass(logging.Logger)
            self.local_logger = logging.getLogger(name)
        finally:
            logging.setLoggerClass(old_logger_class)
        local_handler = logging.FileHandler(LOCAL_ERROR_LOG)
        self.local_logger.addHandler(local_handler)
        self.local_logger.setLevel(level)

    def setLevel(self, level):
        """Override the `setLevel` method to set the same level with both the global and local loggers
        """
        super().setLevel(level)
        self.local_logger.setLevel(level)

    def log(self, level, msg, *args, **kwargs):
        """Override the `log` method to log with both the global and local loggers
        """
        super().log(level, msg, *args, **kwargs)
        self.local_logger.log(level, msg, *args, **kwargs)


def getLogger(name):
    """Gets the `DataConversionLogger` with the provided name if one exists, and creates one if not
    """

    old_logger_class = logging.getLoggerClass()
    try:
        logging.setLoggerClass(DataConversionLogger)
        return logging.getLogger(name)
    finally:
        logging.setLoggerClass(old_logger_class)
