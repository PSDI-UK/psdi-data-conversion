"""@file tests/database_test.py

Created 2025-02-03 by Bryan Gillis.

Unit tests relating to using the database
"""

import pytest

from psdi_data_conversion import constants as const
from psdi_data_conversion.converter import L_REGISTERED_CONVERTERS
from psdi_data_conversion.converters.atomsk import CONVERTER_ATO
from psdi_data_conversion.converters.c2x import CONVERTER_C2X
from psdi_data_conversion.converters.openbabel import CONVERTER_OB
from psdi_data_conversion.database import (FileConverterDatabaseException, disambiguate_formats,
                                           get_conversion_quality, get_converter_info, get_database, get_format_info,
                                           get_in_format_args, get_out_format_args, get_possible_conversions,
                                           get_possible_formats)
from psdi_data_conversion.utils import regularize_name


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

        format_info = get_format_info(name, which=0)

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


def test_format_info_options():
    """Test that we can get the expected information on a few test formats
    """

    # Check that we get an exception for an ambiguous format if we don't request which
    with pytest.raises(FileConverterDatabaseException):
        get_format_info("pdb")

    # Check that requesting all possibilities works as expected
    l_pdb_infos = get_format_info("pdb", which="all")
    assert l_pdb_infos[0] != l_pdb_infos[1]
    assert l_pdb_infos[0] == get_format_info("pdb", which=0)
    assert l_pdb_infos[1] == get_format_info("pdb", which=1)

    # Check that the shortcut for which format works
    assert get_format_info("pdb-0") == l_pdb_infos[0]
    assert get_format_info("pdb-1") == l_pdb_infos[1]

    # Check that the shortcut doesn't cause any problems even if the format is unambiguous
    assert get_format_info("cif-0") == get_format_info("cif")

    # Check that the format info provides the right disambiguated name
    assert get_format_info("cif").disambiguated_name == "cif"
    assert get_format_info("pdb-0").disambiguated_name == "pdb-0"
    assert get_format_info("pdb-1").disambiguated_name == "pdb-1"


def test_disambiguate_format():
    """Test that we can disambiguate formats when only one combination is possible for a conversion
    """

    # Test that we can disambiguate the right pdb format
    in_format, out_format = disambiguate_formats(CONVERTER_OB, "pdb", "cif")
    assert in_format == get_format_info("pdb-0")
    assert out_format == get_format_info("cif")

    # Test we get the expected exception if no conversion is possible
    with pytest.raises(FileConverterDatabaseException, match="is not supported"):
        disambiguate_formats(CONVERTER_C2X, "ins", "cml")

    # Test we get the expected exception if multiple conversions are possible
    with pytest.raises(FileConverterDatabaseException, match="is ambiguous"):
        disambiguate_formats(CONVERTER_C2X, "cif", "pdb")


def test_conversion_table():
    """Test that we can access data from the conversions table properly
    """

    database = get_database()

    # Check the conversions table parent is set properly
    conversions_table = database.conversions_table
    assert conversions_table.parent is database

    # Check we can get the correct conversion quality
    assert get_conversion_quality(CONVERTER_OB, "pdb", "cif").qual_str == const.QUAL_UNKNOWN
    assert get_conversion_quality(CONVERTER_ATO, "xyz", "inchi") is None

    # Do some detailed checks on one conversion
    qual = get_conversion_quality(CONVERTER_OB, "xyz", "inchi")

    assert qual.qual_str == const.QUAL_OKAY

    details = qual.details

    # Check the details are as expected
    assert const.QUAL_NOTE_OUT_MISSING.format(const.QUAL_2D_LABEL) in details
    assert const.QUAL_NOTE_OUT_MISSING.format(const.QUAL_3D_LABEL) in details
    assert const.QUAL_NOTE_IN_MISSING.format(const.QUAL_CONN_LABEL) in details
    assert const.QUAL_COMP_LABEL not in details

    # Check we don't have any extra lines in the details
    assert "\n\n" not in details
    assert not details.startswith("\n")
    assert not details.endswith("\n")

    # Check the property info dict is as expected (mostly covered by details check, so just a couple checks here)
    comp_prop_info = qual.d_prop_conversion_info[const.QUAL_COMP_KEY]
    assert comp_prop_info.input_supported is True
    assert comp_prop_info.output_supported is True
    assert comp_prop_info.label == const.QUAL_COMP_LABEL
    assert comp_prop_info.note == ""

    # Check we can get a list of possible converters for a given conversion
    l_possible_conversions = get_possible_conversions("pdb", "cif")
    assert (regularize_name(CONVERTER_OB), get_format_info("pdb", which=0),
            get_format_info("cif", which=0)) in l_possible_conversions

    # Check that we can get a list of possible input/outpat formats for a given converter
    l_in_formats, l_out_formats = get_possible_formats(CONVERTER_OB)
    assert get_format_info("pdb", which=0) in l_in_formats
    assert get_format_info("cif", which=0) in l_out_formats
