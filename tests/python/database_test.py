"""@file tests/database_test.py

Created 2025-02-03 by Bryan Gillis.

Unit tests relating to using the database
"""

from copy import deepcopy
from uuid import UUID

import pytest

from psdi_data_conversion import constants as const
from psdi_data_conversion import database as db
from psdi_data_conversion.converter import L_SUPPORTED_CONVERTERS
from psdi_data_conversion.testing import constants as tc
from psdi_data_conversion.utils import regularize_name


def test_load():
    """Test that we can load and retrieve the database
    """

    db1 = db.get_database()
    db2 = db.get_database()

    # We should only get one database created, and any additional calls to `get_database()` should return the same
    assert db2 is db1


@pytest.fixture(scope="module")
def database():
    return db.get_database()


def test_converter_info(database):
    """Test that we can get the expected information on each converter
    """

    l_converter_info = db.get_converter_info()
    l_converter_names = [x.name for x in l_converter_info]

    # Check that all supported converters are in this list
    for name in L_SUPPORTED_CONVERTERS:
        assert name in l_converter_names

    for converter_info in l_converter_info:

        if converter_info is None:
            continue

        # Check the name and pretty name are correct
        name = converter_info.name
        assert regularize_name(converter_info.pretty_name) == name

        # Check database is properly set as parent
        assert converter_info.parent == database

        # Check name matches
        assert converter_info.name == name

        # Check ID is of proper type and an allowed value
        assert isinstance(converter_info.id, int)
        assert converter_info.id > 0

        # Check that the UUID matches the ID
        assert converter_info.uuid == UUID(int=converter_info.id)

        # Check description has some text in it
        assert isinstance(converter_info.description, str)
        assert len(converter_info.description) > 0

        # Check URL appears reasonable
        assert isinstance(converter_info.url, str)
        assert "http" in converter_info.url

        # Check that this converter's info can be retrieved through all supported methods
        assert converter_info is db.get_converter_info(converter_info)
        assert converter_info is db.get_converter_info(converter_info.id)
        assert converter_info is db.get_converter_info(UUID(int=converter_info.id))
        assert converter_info is db.get_converter_info(str(UUID(int=converter_info.id)))
        assert converter_info is db.get_converter_info(UUID(int=converter_info.id).hex)
        assert converter_info is db.get_converter_info(converter_info.name)
        assert converter_info is db.get_converter_info(converter_info.pretty_name)


def test_format_args():
    """Test that we can get the flags and options allowed for specific formats for a given converter
    """
    converter_name = const.CONVERTER_OB
    in_format = tc.FORMAT_PDB_0
    out_format = tc.FORMAT_CIF

    l_in_flags, _ = db.get_in_format_args(converter_name, in_format)
    l_out_flags, _ = db.get_out_format_args(converter_name, out_format)

    l_in_flag_names = [x.flag for x in l_in_flags]
    l_out_flag_names = [x.flag for x in l_out_flags]

    assert "b" in l_in_flag_names
    assert "c" in l_in_flag_names
    assert "s" in l_in_flag_names

    assert "g" in l_out_flag_names

    # Check that we can find a specific argument
    in_flag_info_0 = l_in_flags[0]
    assert db.get_in_format_args(converter_name, in_format, in_flag_info_0.flag) is in_flag_info_0
    out_flag_info_0 = l_out_flags[0]
    assert db.get_out_format_args(converter_name, out_format, out_flag_info_0.flag) is out_flag_info_0

    # Check that the UUID is constructed appropriately for the info objects
    assert in_flag_info_0.uuid == UUID(int=in_flag_info_0.id)
    assert out_flag_info_0.uuid == UUID(int=out_flag_info_0.id)


def test_format_info(database):
    """Test that we can get the expected information on a few test formats
    """

    for name, id in (("pdb", tc.FORMAT_PDB_0), ("cif", tc.FORMAT_CIF), ("mmcif", tc.FORMAT_MMCIF),
                     ("inchi", tc.FORMAT_INCHI), ("molreport", tc.FORMAT_MOLREPORT)):

        format_info = db.get_format_info(id)

        # Check database is properly set as parent
        assert format_info.parent == database

        # Check name matches
        assert format_info.name == name

        # Check that the UUID is constructed appropriately
        assert format_info.uuid == UUID(int=format_info.id)

        # Check that this format's info can be retrieved through all supported methods
        assert format_info is db.get_format_info(format_info)
        assert format_info is db.get_format_info(format_info.id)
        assert format_info is db.get_format_info(UUID(int=format_info.id))
        assert format_info is db.get_format_info(str(UUID(int=format_info.id)))
        assert format_info is db.get_format_info(UUID(int=format_info.id).hex)
        assert format_info is db.get_format_info(format_info.name, which=0)
        assert format_info is db.get_format_info(format_info.disambiguated_name)

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
    with pytest.raises(db.FileConverterDatabaseException):
        db.get_format_info("pdb")

    # Check that requesting all possibilities works as expected
    l_pdb_infos = db.get_format_info("pdb", which="all")
    assert l_pdb_infos[0] != l_pdb_infos[1]
    assert l_pdb_infos[0] == db.get_format_info("pdb", which=0)
    assert l_pdb_infos[1] == db.get_format_info("pdb", which=1)

    # Check that the shortcut for which format works
    assert db.get_format_info("pdb-0") == l_pdb_infos[0]
    assert db.get_format_info("pdb-1") == l_pdb_infos[1]

    # Check that the shortcut doesn't cause any problems even if the format is unambiguous
    assert db.get_format_info("cif-0") == db.get_format_info("cif")

    # Check that the format info provides the right disambiguated name
    assert db.get_format_info("cif").disambiguated_name == "cif"
    assert db.get_format_info("pdb-0").disambiguated_name == "pdb-0"
    assert db.get_format_info("pdb-1").disambiguated_name == "pdb-1"

    # Check that formats are case-insensitive
    assert db.get_format_info("PDB-0") is db.get_format_info("pdb-0")


def test_disambiguate_format():
    """Test that we can disambiguate formats when only one combination is possible for a conversion
    """

    # Test that we can disambiguate the right pdb format
    in_format, out_format = db.disambiguate_formats(const.CONVERTER_OB, "pdb", "cif")
    assert in_format == db.get_format_info(tc.FORMAT_PDB_0)
    assert out_format == db.get_format_info(tc.FORMAT_CIF)

    # Test we get the expected exception if no conversion is possible
    with pytest.raises(db.FileConverterDatabaseException, match="is not supported"):
        db.disambiguate_formats(const.CONVERTER_C2X, "ins", "cml")

    # Test we get the expected exception if multiple conversions are possible
    with pytest.raises(db.FileConverterDatabaseException, match="is ambiguous"):
        db.disambiguate_formats(const.CONVERTER_C2X, "cif", "pdb")


def test_conversion_table(database):
    """Test that we can access data from the conversions table properly
    """

    # Check the conversions table parent is set properly
    conversions_table = database.conversions_table
    assert conversions_table.parent is database

    # Check we can get the correct conversion quality
    assert db.get_conversion_quality(const.CONVERTER_OB, tc.FORMAT_PDB_0, tc.FORMAT_CIF).qual_str == const.QUAL_UNKNOWN
    assert db.get_conversion_quality(const.CONVERTER_ATO, tc.FORMAT_XYZ_1, tc.FORMAT_INCHI) is None

    # Do some detailed checks on one conversion
    xyz_format_info = db.get_format_info(tc.FORMAT_XYZ_1)
    inchi_format_info = db.get_format_info(tc.FORMAT_INCHI)

    # "xyz" is ambiguous, but only one possibility has a valid conversion here, so check that we get that one
    qual = db.get_conversion_quality(const.CONVERTER_OB, "xyz", inchi_format_info)

    assert qual.qual_str == const.QUAL_OKAY
    assert qual.in_format is db.get_format_info(xyz_format_info)
    assert qual.out_format is db.get_format_info(inchi_format_info)

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
    l_possible_conversions = db.get_possible_conversions("pdb", "cif")
    assert (db.get_converter_info(const.CONVERTER_OB), db.get_format_info("pdb", which=0),
            db.get_format_info("cif", which=0)) in l_possible_conversions

    # Check that we can get a list of possible input/outpat formats for a given converter
    l_in_formats, l_out_formats = db.get_possible_formats(const.CONVERTER_OB)
    assert db.get_format_info(tc.FORMAT_PDB_0) in l_in_formats
    assert db.get_format_info(tc.FORMAT_CIF) in l_out_formats


def test_conversion_pathways():
    """Tests of determining conversion pathways between formats
    """

    # Check that we get `None` for converting from one format to itself
    assert db.get_conversion_pathway(tc.FORMAT_CIF, tc.FORMAT_CIF) is None

    # Check that we get `None` for an impossible conversion
    assert db.get_conversion_pathway(tc.FORMAT_CIF, tc.FORMAT_ABINIT) is None

    # Check that we get the expected single-step conversion for a known direct conversion
    cif_to_inchi_path = db.get_conversion_pathway(tc.FORMAT_CIF, tc.FORMAT_INCHI, only="registered")
    assert len(cif_to_inchi_path) == 1
    converter_info, in_format_info, out_format_info = cif_to_inchi_path[0]
    assert converter_info.name == regularize_name(const.CONVERTER_OB)
    assert in_format_info.id == tc.FORMAT_CIF
    assert out_format_info.id == tc.FORMAT_INCHI

    # Test getting a multi-step conversion - it's possible this will become direct in the future if a new converter is
    # added, so the test is a bit loose here
    inchi_to_moldy_path = db.get_conversion_pathway(tc.FORMAT_INCHI, tc.FORMAT_MOLDY)
    assert len(inchi_to_moldy_path) <= 2
    assert inchi_to_moldy_path[0][1].id == tc.FORMAT_INCHI
    assert inchi_to_moldy_path[-1][2].id == tc.FORMAT_MOLDY
    for i in range(len(inchi_to_moldy_path)-1):
        # Output format of each step should match input of next
        assert inchi_to_moldy_path[i][2] is inchi_to_moldy_path[i+1][1]
        # Each step should use a different converter
        assert inchi_to_moldy_path[i][0] != inchi_to_moldy_path[i+1][0]


@pytest.fixture()
def format_all(scope="module"):
    return db.FormatInfo("all", database, {key: True for key in db.D_PROP_BITS.keys()})


@pytest.fixture(scope="module")
def format_none(database):
    return db.FormatInfo("none", database, {key: False for key in db.D_PROP_BITS.keys()})


@pytest.fixture(scope="module")
def format_unknown(database):
    return db.FormatInfo("unknown", database, {key: None for key in db.D_PROP_BITS.keys()})


@pytest.fixture(scope="module")
def min_prop_weight():
    return db.STEP_WEIGHT


@pytest.fixture(scope="module")
def max_prop_weight():
    max_weight = db.STEP_WEIGHT
    for bit in db.D_PROP_BITS.values():
        max_weight |= 1 << bit
    return max_weight


@pytest.mark.parametrize("in_format, out_format, ex_weight", [("format_all", "format_all", "min_prop_weight"),
                                                              ("format_all", "format_none", "max_prop_weight"),
                                                              ("format_all", "format_unknown", "max_prop_weight"),
                                                              ("format_none", "format_all", "min_prop_weight"),
                                                              ("format_none", "format_none", "min_prop_weight"),
                                                              ("format_none", "format_unknown", "min_prop_weight"),
                                                              ("format_unknown", "format_all", "min_prop_weight"),
                                                              ("format_unknown", "format_none", "max_prop_weight"),
                                                              ("format_unknown", "format_unknown", "max_prop_weight"),])
def test_get_conversion_prop_weight(in_format, out_format, ex_weight, request):
    """Tests of get_conversion_prop_weight to ensure it calculates weight correctly for whether a format property is
    retained or not in a conversion, with all variations of all properties existing in input and output formats
    """
    assert db.get_conversion_prop_weight(request.getfixturevalue(in_format),
                                         request.getfixturevalue(out_format)) == request.getfixturevalue(ex_weight)


@pytest.mark.parametrize("prop", db.D_PROP_BITS.keys())
def test_get_conversion_prop_weight_prop_lost(format_none, prop):
    """Test each property individually when it's lost to ensure the right bit is set for each"""
    test_format = db.FormatInfo("test", database, {key: True if key == prop else False
                                                   for key in db.D_PROP_BITS.keys()})
    assert db.get_conversion_prop_weight(test_format, format_none) == db.STEP_WEIGHT | 1 << db.D_PROP_BITS[prop]


@pytest.mark.parametrize("in_prec, out_prec, ex_weight", [(None, None, 1 << db.PREC_MAX_DIGIT_LOSS*db.PREC_GAP_BITS),
                                                          (24, None, 1 << db.PREC_MAX_DIGIT_LOSS*db.PREC_GAP_BITS),
                                                          (None, 24, 1 << db.PREC_MAX_DIGIT_LOSS*db.PREC_GAP_BITS),
                                                          (24, 24, 1 << 0*db.PREC_GAP_BITS),
                                                          (24, 18, 1 << 6*db.PREC_GAP_BITS),
                                                          (24, 30, 1 << 0*db.PREC_GAP_BITS),
                                                          (24, 6, 1 << db.PREC_MAX_DIGIT_LOSS*db.PREC_GAP_BITS)])
def test_get_conversion_precision_weight(database, in_prec, out_prec, ex_weight):
    """Test that conversion precision weights are calculated correctly"""
    in_format = db.FormatInfo("in", database, {db.DB_FORMAT_PRECISION_KEY: in_prec})
    out_format = db.FormatInfo("in", database, {db.DB_FORMAT_PRECISION_KEY: out_prec})
    assert db.get_conversion_prec_weight(in_format, out_format) == ex_weight


def test_get_conversion_weight(format_all, format_none, max_prop_weight):
    """Test that getting the full conversion weight is calculated as expected"""
    in_format: db.FormatInfo = deepcopy(format_all)
    in_format.precision = 24
    out_format: db.FormatInfo = deepcopy(format_none)
    out_format.precision = 18

    assert db.get_conversion_weight(in_format, out_format) == ((max_prop_weight << db.PROP_WEIGHT_BIT_OFFSET) +
                                                               (1 << 6*db.PREC_GAP_BITS << db.PREC_WEIGHT_BIT_OFFSET))


@pytest.mark.parametrize("prop_weight, prec_weight, time_weight", [(0, 0, 0),
                                                                   (65535, 65535, 65535),
                                                                   (2788794, 1254542, 2148)])
def test_split_conversion_weight(prop_weight, prec_weight, time_weight):
    """Test that the function to split the conversion weight works as expected"""
    split_weight = db.split_conversion_weight(db.calc_conversion_weight(prop_weight, prec_weight, time_weight))
    assert split_weight.prop_weight == prop_weight
    assert split_weight.prec_weight == prec_weight
    assert split_weight.time_weight == time_weight
