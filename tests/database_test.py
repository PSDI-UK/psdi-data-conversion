"""@file tests/database_test.py

Created 2025-02-03 by Bryan Gillis.

Unit tests relating to using the database
"""


from psdi_data_conversion.converter import L_REGISTERED_CONVERTERS
from psdi_data_conversion.database import get_database


def test_load():
    """Test that we can load and retrieve the database
    """

    db1 = get_database()
    db2 = get_database()

    # We should only get one database created, and any additional calls to `get_database()` should return the same
    assert db2 is db1


def test_open_babel_info():
    """Test that we can get the expected information on each converter
    """

    database = get_database()

    for name in L_REGISTERED_CONVERTERS:

        converter_info = database.d_converter_info[name]

        # Check database is properly set as parent
        assert converter_info.parent == database

        # Check name matches
        assert converter_info.name == name

        # Check ID is of proper type and an allowed value
        assert isinstance(converter_info.id, int)
        assert converter_info.id > 0

        # Check description has some text in it
        assert isinstance(converter_info.description, str)
        assert len(converter_info.description) > 0

        # Check URL appears reasonable
        assert isinstance(converter_info.url, str)
        assert "http" in converter_info.url
