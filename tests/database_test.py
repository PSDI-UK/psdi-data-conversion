"""@file tests/database_test.py

Created 2025-02-03 by Bryan Gillis.

Unit tests relating to using the database
"""


from psdi_data_conversion.database import get_database


def test_load():
    """Test that we can load and retrieve the database
    """

    db1 = get_database()
    db2 = get_database()

    # We should only get one database created, and any additional calls to `get_database()` should return the same
    assert db2 is db1
