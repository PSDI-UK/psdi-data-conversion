"""@file psdi_data_conversion/database.py

Created 2025-02-03 by Bryan Gillis.

Python module provide utilities for accessing the converter database
"""

from __future__ import annotations

import json
import os
from typing import Any

from psdi_data_conversion import constants as const
from psdi_data_conversion.converter import D_REGISTERED_CONVERTERS
from psdi_data_conversion.converters.base import FileConverterException


class FileConverterDatabaseException(FileConverterException):
    """Class for any exceptions which arise from issues with the database classes and methods
    """
    pass


class ConverterInfo:
    """Class providing information on a converter stored in the PSDI Data Conversion database
    """

    def __init__(self,
                 name: str,
                 parent: DataConversionDatabase,
                 d_single_converter_info: dict[str, int | str],
                 d_data: dict[str, Any]):
        """Set up the class - this will be initialised within a `DataConversionDatabase`, which we set as the parent

        Parameters
        ----------
        name : str
            The name of the converter
        parent : DataConversionDatabase
            The database which this belongs to
        d_data : dict[str, Any]
            The loaded database dict
        """

        self.name = name
        self.parent = parent

        # Get info about the converter from the database
        self.id = d_single_converter_info.get(const.DB_ID_KEY, -1)
        self.description = d_single_converter_info.get(const.DB_DESC_KEY, "")
        self.url = d_single_converter_info.get(const.DB_URL_KEY, "")

        # Get necessary info about the converter from the class
        try:
            key_prefix = D_REGISTERED_CONVERTERS[name].database_key_prefix
        except KeyError:
            key_prefix = None

        self._arg_info = {}

        # If the converter class has no defined key prefix, don't add any extra info for it
        if key_prefix is None:
            return
        for key_base in (const.DB_FROM_FLAGS_IN_KEY_BASE,
                         const.DB_FROM_FLAGS_OUT_KEY_BASE,
                         const.DB_FROM_ARGFLAGS_IN_KEY_BASE,
                         const.DB_FROM_ARGFLAGS_OUT_KEY_BASE,
                         const.DB_TO_FLAGS_IN_KEY_BASE,
                         const.DB_TO_FLAGS_OUT_KEY_BASE,
                         const.DB_TO_ARGFLAGS_IN_KEY_BASE,
                         const.DB_TO_ARGFLAGS_OUT_KEY_BASE):
            self._arg_info[key_base] = d_data.get(key_prefix + key_base)


class FormatInfo:
    """Class providing information on a file format from the PSDI Data Conversion database
    """

    def __init__(self,
                 name: str,
                 parent: DataConversionDatabase,
                 d_single_format_info: dict[str, bool | int | str | None]):

        # Load attributes from input
        self.name = name
        self.parent = parent

        # Load attributes from the database
        self.id: int = d_single_format_info.get(const.DB_ID_KEY, -1)
        self.note: str = d_single_format_info.get(const.DB_FORMAT_NOTE_KEY, "")
        self.composition = bool(d_single_format_info.get(const.DB_FORMAT_COMP_KEY, False))
        self.connections = bool(d_single_format_info.get(const.DB_FORMAT_CONN_KEY, False))
        self.two_dim = bool(d_single_format_info.get(const.DB_FORMAT_2D_KEY, False))
        self.three_dim = bool(d_single_format_info.get(const.DB_FORMAT_3D_KEY, False))


class DataConversionDatabase:
    """Class providing interface for information contained in the PSDI Data Conversion database
    """

    def __init__(self, d_data: dict[str, Any]):
        """Initialise the DataConversionDatabase object

        Parameters
        ----------
        d_data : dict[str, Any]
            The dict of the database, as loaded in from the JSON file
        """

        # Store the database dict internally for debugging purposes
        self._d_data = d_data

        # Store top-level items not tied to a specific converter
        self.formats: list[dict[str, bool | int | str | None]] = d_data[const.DB_FORMATS_KEY]
        self.converters: list[dict[str, bool | int | str | None]] = d_data[const.DB_CONVERTERS_KEY]
        self.converts_to: list[dict[str, bool | int | str | None]] = d_data[const.DB_CONVERTS_TO_KEY]

        # Store information on each converter in a dict for it
        self.d_converter_info: dict[str, ConverterInfo] = {}

        for d_single_converter_info in self.converters:
            name: str = d_single_converter_info[const.DB_NAME_KEY]
            self.d_converter_info[name] = ConverterInfo(name, self, d_single_converter_info, d_data)

        # Store information on each format in a dict for it
        self.d_format_info: dict[str, ConverterInfo] = {}

        for d_single_format_info in self.formats:
            name: str = d_single_format_info[const.DB_FORMAT_EXT_KEY]
            self.d_format_info[name] = FormatInfo(name, self, d_single_format_info)


# The database will be loaded on demand when `get_database()` is called
_database: DataConversionDatabase | None = None


def load_database() -> DataConversionDatabase:
    """Load and return a new instance of the data conversion database from the JSON database file in this package. This
    function should not be called directly unless you specifically need a new instance of the database object and can't
    deepcopy the database returned by `get_database()`, as it's expensive to load it in.

    Returns
    -------
    DataConversionDatabase
    """

    # Find and load the database JSON file

    # For an interactive shell, __file__ won't be defined for this module, so use the constants module instead
    reference_file = os.path.realpath(const.__file__)

    qualified_database_filename = os.path.join(os.path.dirname(reference_file), const.DATABASE_FILENAME)
    d_data: dict = json.load(open(qualified_database_filename, "r"))

    return DataConversionDatabase(d_data)


def get_database() -> DataConversionDatabase:
    """Gets the global database object, loading it in first if necessary. Since it's computationally expensive to load
    the database, it's best treated as an immutable singleton.

    Returns
    -------
    DataConversionDatabase
        The global database object
    """
    global _database
    if _database is None:
        # Create the database object and store it globally
        _database = load_database()
    return _database


def get_converter_info(name: str) -> ConverterInfo:
    """Gets the information on a given converter stored in the database

    Parameters
    ----------
    name : str
        The name of the converter

    Returns
    -------
    ConverterInfo
        _description_
    """

    return get_database().d_converter_info[name]


data = get_database()
