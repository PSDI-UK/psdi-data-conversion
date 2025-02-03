"""@file psdi_data_conversion/database.py

Created 2025-02-03 by Bryan Gillis.

Python module provide utilities for accessing the converter database
"""

import json
import os

import psdi_data_conversion
from psdi_data_conversion import constants as const


class DataConversionDatabase(dict):
    pass


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
    qualified_data_base_filename = os.path.join(psdi_data_conversion.__path__, const.DATABASE_FILENAME)
    d_data = json.load(open(qualified_data_base_filename, "r"))

    return DataConversionDatabase(d_data)


def get_database() -> DataConversionDatabase:
    """Gets the global database object, loading it in first if necessary. Since it's computationally expensive to load
    the database, it's best treated as an immutable singleton.

    Returns
    -------
    DataConversionDatabase
        The global database object
    """
    if _database is None:
        # Create the database object and store it globally
        global _database
        _database = load_database()
    return _database
