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


class ConvertsToTable:
    """Class providing information on available file format conversions.

    Information on internal data handling of this class:

    The idea here is that we need to be able to get information on whether a converter can handle a conversion from one
    file format to another. This results in 3D data storage, with dimensions: Converter, Input Format, Output Format.
    The most important operations are (in roughly descending order of importance):

    - For a given Converter, Input Format, and Output Format, get whether or not the conversion is possible, and the
    degree of success if it is possible.
    - For a given Input Format and Output Format, list available Converters and their degrees of success
    - For a given Converter, list available Input Formats and Output Formats
    - For a given Input Format, list available Output Formats and Converters, and the degree of success of each

    At date of implementation, the data comprises 9 Converters and 280 Input/Output Formats, for 705,600 possibilities,
    increasing linearly with the number of converters and quadratically with the number of formats. (Self-to-self format
    conversions don't need to be stored, but this may not be a useful optimisation.)

    Conversion data is available for 23,013 Converter, Input, Output values, or ~3% of the total possible conversions.
    While this could currently work as a sparse array, it will likely be filled to become denser over time, so a dense
    representation makes the most sense.

    The present implementation uses a list-of-lists-of-lists approach, to avoid adding NumPy as a dependency
    until/unless efficiency concerns motivate it in the future.
    """

    def __init__(self,
                 l_converts_to: list[dict[str, bool | int | str | None]],
                 parent: DataConversionDatabase,
                 d_converter_info: dict[str, ConverterInfo],
                 d_format_info: dict[str, FormatInfo]):

        self.parent = parent

        # Store references to needed data
        self._l_converts_to = l_converts_to
        self._d_converter_info = d_converter_info
        self._d_format_info = d_format_info

        # To save data, we store the degree of success in the table as an integer indexed to possible strings
        # `self._l_dos` maps from integer index value to degree of success string, and `self._d_dos` maps from degree of
        # success string to index value
        self._l_dos: list[str | None] = [None] + list(set([x[const.DB_SUCCESS_KEY] for x in l_converts_to]))
        self._d_dos = {dos: index for index, dos in enumerate(self._l_dos)}

        # Build the conversion table, indexed Converter, Input Format, Output Format - note that each of these is
        # 1-indexed, so we add 1 to each of the lengths here
        num_converters = len(parent.converters)
        num_formats = len(parent.formats)

        self._table = [[[0 for k in range(num_formats+1)] for j in range(num_formats+1)]
                       for i in range(num_converters+1)]

        for possible_conversion in l_converts_to:

            conv_id: int = possible_conversion[const.DB_CONV_ID_KEY]
            in_id: int = possible_conversion[const.DB_IN_ID_KEY]
            out_id: int = possible_conversion[const.DB_OUT_ID_KEY]

            self._table[conv_id][in_id][out_id] = self._d_dos[possible_conversion[const.DB_SUCCESS_KEY]]

    def get_degree_of_success(self,
                              converter_name: str,
                              in_format: str,
                              out_format: str) -> str | None:
        """Get the degree of success for a desired conversion, represented as a string (or else None if not possible)

        Parameters
        ----------
        converter_name : str
            The name of the converter to use
        in_format : str
            The extension of the input file format
        out_format : str
            The extension of the output file format

        Returns
        -------
        str | None
            If the conversion is possible, returns a string describing the degree of success. If the conversion is not
            possible, returns None
        """

        conv_id: int = self._d_converter_info[converter_name].id
        in_id: int = self._d_format_info[in_format].id
        out_id: int = self._d_format_info[out_format].id

        return self._l_dos[self._table[conv_id][in_id][out_id]]


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

        # Placeholders for properties that are generated when needed
        self._d_converter_info: dict[str, ConverterInfo] | None = None
        self._d_format_info: dict[str, FormatInfo] | None = None

    @property
    def d_converter_info(self):
        """Generate the converter info dict when needed
        """
        if self._d_converter_info is None:
            self._d_converter_info: dict[str, ConverterInfo] = {}
            for d_single_converter_info in self.converters:
                name: str = d_single_converter_info[const.DB_NAME_KEY]
                self._d_converter_info[name] = ConverterInfo(name, self, d_single_converter_info, self._d_data)
        return self._d_converter_info

    @property
    def d_format_info(self):
        """Generate the format info dict when needed
        """
        if self._d_format_info is None:
            self._d_format_info: dict[str, FormatInfo] = {}

            for d_single_format_info in self.formats:
                name: str = d_single_format_info[const.DB_FORMAT_EXT_KEY]
                self._d_format_info[name] = FormatInfo(name, self, d_single_format_info)
        return self._d_format_info


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
    """

    return get_database().d_converter_info[name]


def get_format_info(name: str) -> FormatInfo:
    """Gets the information on a given file format stored in the database

    Parameters
    ----------
    name : str
        The name (extension) of the form

    Returns
    -------
    FormatInfo
    """

    return get_database().d_format_info[name]


data = get_database()
