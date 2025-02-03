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
from psdi_data_conversion.converters.base import FileConverter, FileConverterException


class FileConverterDatabaseException(FileConverterException):
    """Class for any exceptions which arise from issues with the database
    """
    pass


class ConverterInfo:
    """Class providing information on a converter stored in the database
    """

    def __init__(self,
                 converter_class: type[FileConverter],
                 parent: DataConversionDatabase,
                 d_data: dict[str, Any]):
        """Set up the class - this will be initialised within a `DataConversionDatabase`, which we set as the parent

        Parameters
        ----------
        converter_class : type[FileConverter]
            The class for running conversions with this converter
        parent : DataConversionDatabase
            The database which this belongs to
        d_data : dict[str, Any]
            The loaded database dict
        """

        self.parent = parent

        # Get info about the converter from the class
        self.name = converter_class.name

        # Placeholders for attributes which will be determined as needed
        self._id: int | None = None
        self._description: str | None = None
        self._url: str | None = None

        key_prefix = converter_class.database_key_prefix

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

    @property
    def id(self) -> int:
        """Get the converter's ID in the database, determining it first if needed.
        """
        if self._id is None:
            self._get_converter_general()
        return self._id

    def _get_converter_general(self):
        """Finds the converter's general info in the database and stores it in this class's member variables
        """
        # Search through the list of converters to find the one which has the name of this one
        l_matching_converters = [x for x in self.parent.converters if x['name'] == self.name]

        # Check we find exactly one
        if len(l_matching_converters) == 0:
            raise FileConverterDatabaseException(f"Converter {self.name} not found in database's list of converters")
        elif len(l_matching_converters) > 1:
            raise FileConverterDatabaseException(f"Converter {self.name} appears multiple times in database. Found "
                                                 "entries are: \n" + "\n".join(str(l_matching_converters)))

        d_converter_info = l_matching_converters[0]

        self._id = d_converter_info['id']
        self._description = d_converter_info['description']
        self._url = d_converter_info['url']


class DataConversionDatabase:
    """Class providing interface for information contained in the PSDI Data Conversion Database
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
        self.formats: list[dict[str, int | str | None]] = d_data[const.DB_FORMATS_KEY]
        self.converters: list[dict[str, int | str | None]] = d_data[const.DB_CONVERTERS_KEY]
        self.converts_to: list[dict[str, int | str | None]] = d_data[const.DB_CONVERTS_TO_KEY]

        # Store information on each converter in a dict for it
        self.converter_info = {}

        for converter_class in D_REGISTERED_CONVERTERS.values():
            self.converter_info[converter_class.name] = ConverterInfo(converter_class, self, d_data)


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


data = get_database()
