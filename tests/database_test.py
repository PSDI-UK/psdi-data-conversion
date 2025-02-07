"""@file tests/database_test.py

Created 2025-02-03 by Bryan Gillis.

Unit tests relating to using the database
"""

from psdi_data_conversion import constants as const
from psdi_data_conversion.converter import L_REGISTERED_CONVERTERS
from psdi_data_conversion.converters.openbabel import CONVERTER_OB
from psdi_data_conversion.database import (get_conversion_quality, get_converter_info, get_database, get_format_info,
                                           get_in_format_args, get_out_format_args, get_possible_converters,
                                           get_possible_formats)


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


def test_format_args():
    """Test that we can get the flags and options allowed for specific formats for a given converter
    """
    converter_name = CONVERTER_OB
    in_format = "pdb"
    out_format = "cif"

    l_in_flags, _ = get_in_format_args(converter_name, in_format)
    l_out_flags, _ = get_out_format_args(converter_name, out_format)

    l_in_flag_names = [x.flag for x in l_in_flags]
    l_out_flag_names = [x.flag for x in l_out_flags]

    assert "b" in l_in_flag_names
    assert "c" in l_in_flag_names
    assert "s" in l_in_flag_names

    assert "g" in l_out_flag_names

    # Check that we can find a specific argument
    in_flag_info_0 = l_in_flags[0]
    assert get_in_format_args(converter_name, in_format, in_flag_info_0.flag) is in_flag_info_0
    out_flag_info_0 = l_out_flags[0]
    assert get_out_format_args(converter_name, out_format, out_flag_info_0.flag) is out_flag_info_0


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

        if name in ("pdb", "mmcif", "inchi", "molreport"):
            assert format_info.composition, name
        else:
            assert not format_info.composition, name

        if name in ("pdb", "inchi", "molreport"):
            assert format_info.connections, name
        else:
            assert not format_info.connections, name

        if name in ("mmcif", "molreport"):
            assert format_info.two_dim, name
        else:
            assert not format_info.two_dim, name

        if name in ("pdb", "mmcif"):
            assert format_info.three_dim, name
        else:
            assert not format_info.three_dim, name


def test_conversion_table():
    """Test that we can access data from the conversions table properly
    """

    database = get_database()

    # Check the conversions table parent is set properly
    conversions_table = database.conversions_table
    assert conversions_table.parent is database

    # Check we can get the correct degree of success
    assert get_conversion_quality(CONVERTER_OB, "pdb", "cif") == const.QUAL_UNKNOWN

    # Check we can get a list of possible converters for a given conversion
    l_possible_converters = get_possible_converters("pdb", "cif")
    assert CONVERTER_OB in [name for name in l_possible_converters]

    # Check that we can get a list of possible input/outpat formats for a given converter
    l_in_formats, l_out_formats = get_possible_formats(CONVERTER_OB)
    assert "pdb" in l_in_formats
    assert "cif" in l_out_formats
