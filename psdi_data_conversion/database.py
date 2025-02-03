"""@file psdi_data_conversion/database.py

Created 2025-02-03 by Bryan Gillis.

Python module provide utilities for accessing the converter database
"""

import json
import os

from psdi_data_conversion import constants as const


class DataConversionDatabase:
    def __init__(self, d_data: dict | None = None):
        if d_data is None:
            d_data = {}


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
