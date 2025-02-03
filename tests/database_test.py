"""@file tests/database_test.py

Created 2025-02-03 by Bryan Gillis.

Unit tests relating to using the database
"""


from psdi_data_conversion.converter import L_REGISTERED_CONVERTERS
from psdi_data_conversion.database import get_converter_info, get_database, get_format_info


def test_load():
    """Test that we can load and retrieve the database
    """

    db1 = get_database()
    db2 = get_database()

    # We should only get one database created, and any additional calls to `get_database()` should return the same
    assert db2 is db1


def test_converter_info():
    """Test that we can get the expected information on each converter
    """

    database = get_database()

    for name in L_REGISTERED_CONVERTERS:

        converter_info = get_converter_info(name)

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


def test_format_info():
    """Test that we can get the expected information on a few test formats
    """

    database = get_database()

    for name in ("pdb", "cif", "mmcif", "inchi", "molreport"):

        format_info = get_format_info(name)

        # Check database is properly set as parent
        assert format_info.parent == database

        # Check name matches
        assert format_info.name == name

        # Check properties are as expected

        if name in ("mmcif", "inchi", "molreport"):
            assert format_info.composition, name
        else:
            assert not format_info.composition, name

        if name in ("inchi", "molreport"):
            assert format_info.connections, name
        else:
            assert not format_info.connections, name

        if name in ("mmcif", "molreport"):
            assert format_info.two_dim, name
        else:
            assert not format_info.two_dim, name

        if name == "mmcif":
            assert format_info.three_dim, name
        else:
            assert not format_info.three_dim, name
