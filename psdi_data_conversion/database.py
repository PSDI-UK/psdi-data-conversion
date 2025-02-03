"""@file psdi_data_conversion/database.py

Created 2025-02-03 by Bryan Gillis.

Python module provide utilities for accessing the converter database
"""

import json
import os

from psdi_data_conversion import constants as const
from psdi_data_conversion.converter import D_REGISTERED_CONVERTERS


class DataConversionDatabase:
    """Class providing interface for information contained in the PSDI Data Conversion Database
    """

    def __init__(self, d_data: dict):
        """Initialise the DataConversionDatabase object

        Parameters
        ----------
        d_data : dict
            The dict of the database, as loaded in from the JSON file
        """

        # Store the database dict internally for debugging purposes
        self._d_data = d_data

        # Store top-level items not tied to a specific converter
        self.formats: dict = d_data[const.DB_FORMATS_KEY]
        self.converters: dict = d_data[const.DB_CONVERTERS_KEY]
        self.converts_to: dict = d_data[const.DB_CONVERTS_TO_KEY]

        # Dict to store info specific to converters
        self.converter_info = {}

        for converter_class in D_REGISTERED_CONVERTERS.values():

            key_prefix = converter_class.database_key_prefix

            # If the converter class has no defined key prefix, don't add any info for it
            if key_prefix is None:
                continue

            new_converter_info = {}

            for key_base in (const.DB_FROM_FLAGS_IN_KEY_BASE,
                             const.DB_FROM_FLAGS_OUT_KEY_BASE,
                             const.DB_FROM_ARGFLAGS_IN_KEY_BASE,
                             const.DB_FROM_ARGFLAGS_OUT_KEY_BASE,
                             const.DB_TO_FLAGS_IN_KEY_BASE,
                             const.DB_TO_FLAGS_OUT_KEY_BASE,
                             const.DB_TO_ARGFLAGS_IN_KEY_BASE,
                             const.DB_TO_ARGFLAGS_OUT_KEY_BASE):
                new_converter_info[key_base] = d_data.get(key_prefix + key_base)

            self.converter_info[converter_class.name] = new_converter_info


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
