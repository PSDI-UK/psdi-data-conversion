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


@pytest.fixture(scope="module", autouse=True)
def database():
    return db.get_database()


def test_load(database):
    """Test that we can load and retrieve the database
    """

    new_database = db.get_database()

    # We should only get one database created, and any additional calls to `get_database()` should return the same
    assert new_database is database


@pytest.fixture(scope="module")
def l_converter_info():
    return db.get_converter_info()


@pytest.mark.parametrize("converter", L_SUPPORTED_CONVERTERS)
def test_converter_info_supported(converter, l_converter_info: list[db.ConverterInfo]):
    """Test that info is available for all supported converters"""
    l_converter_names = [x.name for x in l_converter_info]
    assert converter in l_converter_names


@pytest.mark.parametrize("converter", L_SUPPORTED_CONVERTERS)
def test_converter_info_valid(converter,
                              database,
                              subtests):
    """Test that we can get the expected information on each converter
    """

    converter_info: db.ConverterInfo = db.get_converter_info(converter)

    with subtests.test("Test converter info is found"):
        assert converter_info is not None

    name = converter_info.name

    with subtests.test("Test pretty name matches name if regularized"):
        assert regularize_name(converter_info.pretty_name) == name

    with subtests.test("Test the name matches the name used in the input list"):
        assert converter == name

    with subtests.test("Test database is properly set as parent"):
        assert converter_info.parent == database

    with subtests.test("Test ID is of proper type and an allowed value"):
        assert isinstance(converter_info.id, int)
        assert converter_info.id > 0

    with subtests.test("Test that the UUID matches the ID"):
        assert converter_info.uuid == UUID(int=converter_info.id)

    with subtests.test("Test description has some text in it"):
        assert isinstance(converter_info.description, str)
        assert len(converter_info.description) > 0

    with subtests.test("Test URL appears reasonable"):
        assert isinstance(converter_info.url, str)
        assert "http" in converter_info.url

    with subtests.test("Test get info from itself"):
        assert converter_info is db.get_converter_info(converter_info)
    with subtests.test("Test get info from ID"):
        assert converter_info is db.get_converter_info(converter_info.id)
    with subtests.test("Test get info from UUID"):
        assert converter_info is db.get_converter_info(UUID(int=converter_info.id))
    with subtests.test("Test get info from UUID as string"):
        assert converter_info is db.get_converter_info(str(UUID(int=converter_info.id)))
    with subtests.test("Test get info from UUID as hex"):
        assert converter_info is db.get_converter_info(UUID(int=converter_info.id).hex)
    with subtests.test("Test get info from pretty name"):
        assert converter_info is db.get_converter_info(converter_info.pretty_name)


def test_format_args(subtests):
    """Test that we can get the flags and options allowed for specific formats for a given converter
    """
    converter_name = const.CONVERTER_OB
    in_format = tc.FORMAT_PDB_0
    out_format = tc.FORMAT_CIF

    with subtests.test("Test in flags are correct"):
        l_in_flags, _ = db.get_in_format_args(converter_name, in_format)
        l_in_flag_names = [x.name for x in l_in_flags]
        assert "b" in l_in_flag_names
        assert "c" in l_in_flag_names
        assert "s" in l_in_flag_names

    with subtests.test("Test out flags are correct"):
        l_out_flags, _ = db.get_out_format_args(converter_name, out_format)
        l_out_flag_names = [x.name for x in l_out_flags]
        assert "g" in l_out_flag_names

    with subtests.test("Test that we can find a specific in flag"):
        in_flag_info_0 = l_in_flags[0]
        assert db.get_in_format_args(converter_name, in_format, in_flag_info_0.name) is in_flag_info_0

    with subtests.test("Test that we can find a specific out flag"):
        out_flag_info_0 = l_out_flags[0]
        assert db.get_out_format_args(converter_name, out_format, out_flag_info_0.name) is out_flag_info_0

    with subtests.test("Test in flag UUID is correctly constructed"):
        assert in_flag_info_0.uuid == UUID(int=in_flag_info_0.id)
    with subtests.test("Test out flag UUID is correctly constructed"):
        assert out_flag_info_0.uuid == UUID(int=out_flag_info_0.id)


@pytest.mark.parametrize("name, id", (("pdb", tc.FORMAT_PDB_0), ("inchikey", tc.FORMAT_INCHIKEY),
                         ("mmcif", tc.FORMAT_MMCIF), ("inchi", tc.FORMAT_INCHI), ("molreport", tc.FORMAT_MOLREPORT)))
def test_format_info(name, id, database, subtests):
    """Test that we can get the expected information on a few test formats
    """

    format_info = db.get_format_info(id)

    with subtests.test("Test database is properly set as parent"):
        assert format_info.parent == database

    with subtests.test("Test name matches"):
        assert format_info.name == name

    with subtests.test("Test that the UUID is constructed appropriately"):
        assert format_info.uuid == UUID(int=format_info.id)

    # Check that this format's info can be retrieved through all supported methods
    with subtests.test("Test get format info from itself"):
        assert format_info is db.get_format_info(format_info)
    with subtests.test("Test get format info from ID"):
        assert format_info is db.get_format_info(format_info.id)
    with subtests.test("Test get format info from UUID"):
        assert format_info is db.get_format_info(UUID(int=format_info.id))
    with subtests.test("Test get format info from UUID as string"):
        assert format_info is db.get_format_info(str(UUID(int=format_info.id)))
    with subtests.test("Test get format info from UUID as hex"):
        assert format_info is db.get_format_info(UUID(int=format_info.id).hex)
    with subtests.test("Test get format info from name"):
        assert format_info is db.get_format_info(format_info.name, which=0)
    with subtests.test("Test get format info from disambiguated name"):
        assert format_info is db.get_format_info(format_info.disambiguated_name)

    # Check properties are as expected

    with subtests.test("Test composition property is correct"):
        if name in ("pdb", "mmcif", "inchi", "molreport"):
            assert format_info.composition, name
        else:
            assert not format_info.composition, name

    with subtests.test("Test connections property is correct"):
        if name in ("pdb", "inchi", "molreport"):
            assert format_info.connections, name
        else:
            assert not format_info.connections, name

    with subtests.test("Test 2D and 3D properties are correct"):
        if name in ("pdb", "mmcif"):
            assert format_info.two_dim, name
            assert format_info.three_dim, name
        else:
            assert not format_info.two_dim, name
            assert not format_info.three_dim, name

# "ent" is an alias of the PDB format info. Check various aspects of each to ensure they work correctly


@pytest.fixture(scope="module")
def pdb_format_info():
    return db.get_format_info(tc.FORMAT_PDB_0)


@pytest.fixture(scope="module")
def ent_format_info():
    return db.get_format_info(tc.FORMAT_ENT)


def test_format_alias_equality(pdb_format_info: db.FormatInfo, ent_format_info: db.FormatInfo, subtests):
    """Test that format aliases are separate but share common info"""

    with subtests.test("Test alias format infos aren't the same object"):
        assert ent_format_info != pdb_format_info

    with subtests.test("Test alias format infos share the same common format info"):
        assert ent_format_info.format_common_info is pdb_format_info.format_common_info


def test_format_alias_primary(pdb_format_info: db.FormatInfo, ent_format_info: db.FormatInfo, subtests):
    """Test that format aliases properly indicate which is the primary in all appropriate ways"""

    with subtests.test("Test the alias format info isn't labelled as primary"):
        assert not ent_format_info.is_primary
    with subtests.test("Test the primary format info is labelled as primary"):
        assert pdb_format_info.is_primary

    with subtests.test("Test the alias format references the primary as primary ID"):
        assert ent_format_info.primary_id != ent_format_info.id
        assert ent_format_info.primary_id == pdb_format_info.id

    with subtests.test("Test the alias format references the primary as primary name"):
        assert ent_format_info.primary_name != ent_format_info.name
        assert ent_format_info.primary_name == pdb_format_info.name

    with subtests.test("Test the primary format references itself as primary ID"):
        assert pdb_format_info.primary_id == pdb_format_info.id

    with subtests.test("Test the primary format references itself as primary name"):
        assert pdb_format_info.primary_name == pdb_format_info.name


def test_format_alias_dicts(pdb_format_info: db.FormatInfo, ent_format_info: db.FormatInfo, subtests):
    """Test that format aliases dicts of IDs and extensions behave appropriately"""

    with subtests.test("Test the alias dicts are the same"):
        assert ent_format_info.d_alias_exts is pdb_format_info.d_alias_exts

    with subtests.test("Test the alias dict includes alias format"):
        assert ent_format_info.d_alias_exts[ent_format_info.id] == ent_format_info.name
    with subtests.test("Test the alias dict includes primary format"):
        assert ent_format_info.d_alias_exts[pdb_format_info.id] == pdb_format_info.name


def test_format_alias_graph(database: db.DataConversionDatabase, subtests):
    """Test that format aliases are handled correctly in the graphs"""
    d_indices_from_uuids = database.conversions_table.d_indices_from_uuids
    d_uuids_from_indices = database.conversions_table.d_uuids_from_indices

    with subtests.test("Test both formats point to same vertex"):
        assert d_indices_from_uuids[tc.FORMAT_PDB_0] == d_indices_from_uuids[tc.FORMAT_ENT]

    with subtests.test("Test the vertex points back to the primary format"):
        assert d_uuids_from_indices[d_indices_from_uuids[tc.FORMAT_ENT]] == tc.FORMAT_PDB_0


def test_format_info_options(subtests):
    """Test that we can get the expected information on a few test formats
    """

    with subtests.test("Ambiguous format raises an error"), pytest.raises(db.FileConverterDatabaseException):
        db.get_format_info("pdb")

    # Check that requesting all possibilities works as expected
    l_pdb_infos = db.get_format_info("pdb", which="all")
    with subtests.test("Returned possibilities are different"):
        assert l_pdb_infos[0] != l_pdb_infos[1]
    for i in range(len(l_pdb_infos)):
        with subtests.test("Can get each returned possibility with `which`", i=i):
            assert l_pdb_infos[i] == db.get_format_info("pdb", which=i)
        with subtests.test("Disambiguated name for each format works"):
            assert db.get_format_info(f"pdb-{i}") == l_pdb_infos[i]

    with subtests.test("Disambiguated name doesn't cause any problems even if the format is unambiguous"):
        assert db.get_format_info("cif-0") == db.get_format_info("cif")

    with subtests.test("Right disambiguated name for unambiguous format"):
        assert db.get_format_info("cif").disambiguated_name == "cif"
    for i in range(len(l_pdb_infos)):
        with subtests.test("Right disambiguated name for ambiguous format", i=i):
            assert db.get_format_info(f"pdb-{i}").disambiguated_name == f"pdb-{i}"

    with subtests.test("Formats are case-insensitive"):
        assert db.get_format_info("PDB-0") is db.get_format_info("pdb-0")


def test_disambiguate_format(subtests):
    """Test that we can disambiguate formats when only one combination is possible for a conversion
    """

    with subtests.test("Can disambiguate only PDB format supported by Open Babel"):
        in_format, out_format = db.disambiguate_formats(const.CONVERTER_OB, "pdb", "cif")
        assert (in_format, out_format) == (db.get_format_info(tc.FORMAT_PDB_0), db.get_format_info(tc.FORMAT_CIF))

    with (subtests.test("Error if no conversion is possible"),
          pytest.raises(db.FileConverterDatabaseException, match="is not supported")):
        db.disambiguate_formats(const.CONVERTER_C2X, "ins", "cml")

    with (subtests.test("Error if multiple conversions are possible"),
          pytest.raises(db.FileConverterDatabaseException, match="is ambiguous")):
        db.disambiguate_formats(const.CONVERTER_C2X, "cif", "pdb")


def test_conversion_table(database: db.DataConversionDatabase, subtests):
    """Test that we can access data from the conversions table properly
    """

    conversions_table = database.conversions_table

    with subtests.test("Table parent is set up properly"):
        assert conversions_table.parent is database

    with subtests.test("Get quality for possible conversion"):
        assert db.get_conversion_quality(const.CONVERTER_OB, tc.FORMAT_PDB_0,
                                         tc.FORMAT_CIF).qual_str == const.QUAL_VERYGOOD
    with subtests.test("Get None quality for impossible conversion"):
        assert db.get_conversion_quality(const.CONVERTER_ATO, tc.FORMAT_XYZ_1, tc.FORMAT_INCHI) is None

    # Do some detailed checks on one conversion
    xyz_format_info = db.get_format_info(tc.FORMAT_XYZ_1)
    inchi_format_info = db.get_format_info(tc.FORMAT_INCHI)

    # "xyz" is ambiguous, but only one possibility has a valid conversion here, so check that we get that one
    qual = db.get_conversion_quality(const.CONVERTER_OB, "xyz", inchi_format_info)

    with subtests.test("Quality is okay as expected"):
        assert qual.qual_str == const.QUAL_OKAY
    with subtests.test("Quality references proper in format"):
        assert qual.in_format is db.get_format_info(xyz_format_info)
    with subtests.test("Quality references proper out format"):
        assert qual.out_format is db.get_format_info(inchi_format_info)

    details = qual.details

    # Check the details are as expected
    with subtests.test("Quality notes that 2D is missing in out format"):
        assert const.QUAL_NOTE_OUT_MISSING.format(const.QUAL_2D_LABEL) in details
    with subtests.test("Quality notes that 3D is missing in out format"):
        assert const.QUAL_NOTE_OUT_MISSING.format(const.QUAL_3D_LABEL) in details
    with subtests.test("Quality notes that connections is missing in in format"):
        assert const.QUAL_NOTE_IN_MISSING.format(const.QUAL_CONN_LABEL) in details
    with subtests.test("Quality notes don't mention composition"):
        assert const.QUAL_COMP_LABEL not in details

    with subtests.test("No double line breaks in details"):
        assert "\n\n" not in details
    with subtests.test("Details doesn't start with a line break"):
        assert not details.startswith("\n")
    with subtests.test("Details doesn't end with a line break"):
        assert not details.endswith("\n")

    # Check the property info dict is as expected (mostly covered by details check, so just a couple checks here)
    comp_prop_info = qual.d_prop_conversion_info[const.QUAL_COMP_KEY]
    with subtests.test("Quality info dict notes composition is supported in in format"):
        assert comp_prop_info.input_supported is True
    with subtests.test("Quality info dict notes composition is supported in out format"):
        assert comp_prop_info.output_supported is True
    with subtests.test("Quality info dict has correct label for composition"):
        assert comp_prop_info.label == const.QUAL_COMP_LABEL
    with subtests.test("Quality info dict has correct description for composition"):
        assert comp_prop_info.description == ""

    with subtests.test("Get possible conversion from database method for it"):
        l_possible_conversions = db.get_possible_conversions("pdb", "cif")
        assert (db.get_converter_info(const.CONVERTER_OB), db.get_format_info("pdb", which=0),
                db.get_format_info("cif", which=0)) in l_possible_conversions

    # Check that we can get a list of possible input/outpat formats for a given converter
    with subtests.test("Get expected input/output formats for converter"):
        l_in_formats, l_out_formats = db.get_possible_formats(const.CONVERTER_OB)
        with subtests.test("List contains an expected input format"):
            assert db.get_format_info(tc.FORMAT_PDB_0) in l_in_formats
        with subtests.test("List contains an expected output format"):
            assert db.get_format_info(tc.FORMAT_CIF) in l_out_formats


def test_conversion_pathway_to_self():
    """Test that we get `None` for converting from one format to itself"""
    assert db.get_conversion_pathway(tc.FORMAT_CIF, tc.FORMAT_CIF) is None


def test_conversion_pathway_impossible():
    """Test that we get `None` for an impossible conversion"""
    assert db.get_conversion_pathway(tc.FORMAT_CIF, tc.FORMAT_ABINIT) is None


def test_conversion_pathway_direct(subtests):
    """Test that we get the expected single-step conversion for a known direct conversion"""
    with subtests.test("Get a single-step path for a direct conversion"):
        cif_to_inchi_path = db.get_conversion_pathway(tc.FORMAT_CIF, tc.FORMAT_INCHI)
        assert len(cif_to_inchi_path) == 1

    # Check this step is a valid conversion
    with subtests.test("Step represents a valid conversion"):
        step = cif_to_inchi_path[0]
        assert step.is_valid()

    converter_info, in_format_info, out_format_info = step
    with subtests.test("Converter info in step is correct"):
        assert converter_info.name == regularize_name(const.CONVERTER_OB)
    with subtests.test("In format info in step is correct"):
        assert in_format_info.id == tc.FORMAT_CIF
    with subtests.test("Out format info in step is correct"):
        assert out_format_info.id == tc.FORMAT_INCHI


def _check_path_valid(path: db.ConversionPath, subtests):
    """Check that a path is valid and each step uses a different converter"""
    with subtests.test("Path is valid", name=path.get_name()):
        assert path.is_valid()
    with subtests.test("All converters in path unique", name=path.get_name()):
        s_converters = {step.converter for step in path}
        assert len(s_converters) == len(path)


@pytest.fixture(scope="module")
def inchi_to_moldy_path() -> db.ConversionPath:
    return db.get_conversion_pathway(tc.FORMAT_INCHI, tc.FORMAT_MOLDY)


def test_conversion_pathway_multistep(inchi_to_moldy_path: db.ConversionPath, subtests):
    """Test getting a multi-step conversion - it's possible this will become direct in the future if a new converter is
    added, so the test is a bit loose here"""

    with subtests.test("Path is as short as expected"):
        assert len(inchi_to_moldy_path) <= 2

    _check_path_valid(inchi_to_moldy_path, subtests)


def test_conversion_pathway_from_alias(subtests):
    """Test that if a conversion is requested from an alias, that alias is retained in the input path"""

    from_alias_path = db.get_conversion_pathway(tc.FORMAT_MOLD_ALIAS, tc.FORMAT_MOLDY)
    assert len(from_alias_path) > 1, "Test is only valid if path has at least 2 steps"

    with subtests.test("In format of path is the alias requested"):
        assert from_alias_path[0][1].id == tc.FORMAT_MOLD_ALIAS

    _check_path_valid(from_alias_path, subtests)


def test_conversion_pathway_to_alias(subtests):
    """Test that if a conversion is requested to an alias, that alias is retained in the output path"""

    to_alias_path = db.get_conversion_pathway(tc.FORMAT_MOLDY, tc.FORMAT_MOLD_ALIAS)
    assert len(to_alias_path) > 1, "Test is only valid if path has at least 2 steps"

    with subtests.test("Out format of path is the alias requested"):
        assert to_alias_path[-1][2].id == tc.FORMAT_MOLD_ALIAS

    _check_path_valid(to_alias_path, subtests)


@pytest.fixture(scope="module")
def l_best_inchi_to_moldy_paths():
    return db.get_possible_conversion_pathways(tc.FORMAT_INCHI, tc.FORMAT_MOLDY, include="best")


def test_conversion_pathways_best(l_best_inchi_to_moldy_paths: list[db.ConversionPath], subtests):
    """Test that we can successfully get a list of all equally-low-weight conversion pathways for a desired conversion
    """
    weight = None
    for path in l_best_inchi_to_moldy_paths:
        _check_path_valid(path, subtests)
        if weight is None:
            weight = path.get_weight()
        else:
            with subtests.test("All paths have the same weight", name=path.get_name()):
                assert path.get_weight() == weight


@pytest.fixture(scope="module")
def l_shortest_inchi_to_moldy_paths():
    return db.get_possible_conversion_pathways(tc.FORMAT_INCHI, tc.FORMAT_MOLDY, include="shortest")


def test_conversion_pathways_shortest(l_shortest_inchi_to_moldy_paths: list[db.ConversionPath],
                                      inchi_to_moldy_path: db.ConversionPath, subtests):
    """Test that we can successfully get a list of all equally-short conversion pathways for a desired conversion
    """
    path_length = None
    min_weight = inchi_to_moldy_path.get_weight()
    for path in l_shortest_inchi_to_moldy_paths:
        _check_path_valid(path, subtests)
        _check_path_valid(path, subtests)
        if path_length is None:
            path_length = len(path)
        else:
            with subtests.test("All paths have the same length", name=path.get_name()):
                assert path_length == len(path)
            with subtests.test("All paths have equal or more weight than best", name=path.get_name()):
                assert path.get_weight() >= min_weight


def test_conversion_pathways_different_amounts(inchi_to_moldy_path: db.ConversionPath,
                                               l_best_inchi_to_moldy_paths: list[db.ConversionPath],
                                               l_shortest_inchi_to_moldy_paths: list[db.ConversionPath],
                                               subtests):
    """Test that the different methods of getting paths give sane results - that the one path is one of the best paths,
    and the best paths are all included in the shortest paths"""

    with subtests.test("Single path in best paths"):
        assert inchi_to_moldy_path in l_best_inchi_to_moldy_paths

    lowest_weight = inchi_to_moldy_path.get_weight()
    shortest_len = len(inchi_to_moldy_path)

    for i, best_path in enumerate(l_best_inchi_to_moldy_paths):
        with subtests.test("Best paths have min weight", i=i):
            assert best_path.get_weight() == lowest_weight
        with subtests.test("Best paths in shortest paths", i=i):
            assert best_path in l_shortest_inchi_to_moldy_paths

    for i, shortest_path in enumerate(l_shortest_inchi_to_moldy_paths):
        with subtests.test("Shortest paths have shortest length", i=i):
            assert len(shortest_path) == shortest_len
        if shortest_path not in l_best_inchi_to_moldy_paths:
            with subtests.test("Shortest paths not in best have higher weight", i=i):
                assert shortest_path.get_weight() > lowest_weight


@pytest.fixture(scope="module")
def format_all(database):
    return db.FormatInfo.from_db(database, {db.DB_NAME_KEY: "all", **{key: True for key in db.D_PROP_BITS.keys()}})


@pytest.fixture(scope="module")
def format_none(database):
    return db.FormatInfo.from_db(database, {db.DB_NAME_KEY: "none", **{key: False for key in db.D_PROP_BITS.keys()}})


@pytest.fixture(scope="module")
def format_unknown(database):
    return db.FormatInfo.from_db(database, {db.DB_NAME_KEY: "unknown", **{key: None for key in db.D_PROP_BITS.keys()}})


@pytest.fixture(scope="module")
def min_prop_weight():
    return 0


@pytest.fixture(scope="module")
def max_prop_weight():
    max_weight = 0
    for bit in db.D_PROP_BITS.values():
        max_weight |= 1 << bit
    return max_weight


@pytest.fixture(scope="module")
def converter_ob():
    return db.get_converter_info("Open Babel")


@pytest.mark.parametrize("in_format, out_format, ex_weight", [("format_all", "format_all", "min_prop_weight"),
                                                              ("format_all", "format_none", "max_prop_weight"),
                                                              ("format_all", "format_unknown", "max_prop_weight"),
                                                              ("format_none", "format_all", "min_prop_weight"),
                                                              ("format_none", "format_none", "min_prop_weight"),
                                                              ("format_none", "format_unknown", "min_prop_weight"),
                                                              ("format_unknown", "format_all", "min_prop_weight"),
                                                              ("format_unknown", "format_none", "max_prop_weight"),
                                                              ("format_unknown", "format_unknown", "max_prop_weight"),])
def test_calc_conversion_prop_weight(in_format, out_format, ex_weight,
                                     converter_ob, request):
    """Tests of get_conversion_prop_weight to ensure it calculates weight correctly for whether a format property is
    retained or not in a conversion, with all variations of all properties existing in input and output formats
    """
    assert db.calc_conversion_prop_weight(converter_ob,
                                          request.getfixturevalue(in_format),
                                          request.getfixturevalue(out_format)) == request.getfixturevalue(ex_weight)


@pytest.mark.parametrize("prop", db.D_PROP_BITS.keys())
def test_calc_conversion_prop_weight_prop_lost(format_none, prop, converter_ob):
    """Test each property individually when it's lost to ensure the right bit is set for each"""
    test_format = db.FormatInfo.from_db(database, {db.DB_NAME_KEY: "unknown", **{key: True if key == prop else False
                                                                                 for key in db.D_PROP_BITS.keys()}})
    assert db.calc_conversion_prop_weight(converter_ob, test_format,
                                          format_none) == 1 << db.D_PROP_BITS[prop]


@pytest.mark.parametrize("in_prec, out_prec, ex_weight", [(None, None, 1 << db.PREC_MAX_DIGIT_LOSS*db.PREC_GAP_BITS),
                                                          (24, None, 1 << db.PREC_MAX_DIGIT_LOSS*db.PREC_GAP_BITS),
                                                          (None, 24, 1 << db.PREC_MAX_DIGIT_LOSS*db.PREC_GAP_BITS),
                                                          (24, 24, 1 << 0*db.PREC_GAP_BITS),
                                                          (24, 18, 1 << 6*db.PREC_GAP_BITS),
                                                          (24, 30, 1 << 0*db.PREC_GAP_BITS),
                                                          (24, 6, 1 << db.PREC_MAX_DIGIT_LOSS*db.PREC_GAP_BITS)])
def test_calc_conversion_precision_weight(database, in_prec, out_prec, ex_weight, converter_ob):
    """Test that conversion precision weights are calculated correctly"""
    in_format = db.FormatInfo.from_db(database, {db.DB_NAME_KEY: "in", db.DB_FORMAT_PRECISION_KEY: in_prec}, {})
    out_format = db.FormatInfo.from_db(database, {db.DB_NAME_KEY: "out", db.DB_FORMAT_PRECISION_KEY: out_prec}, {})
    assert db.calc_conversion_prec_weight(converter_ob, in_format, out_format) == ex_weight


def test_calc_conversion_weight(format_all, format_none, max_prop_weight, converter_ob):
    """Test that getting the full conversion weight is calculated as expected"""
    in_format: db.FormatInfo = deepcopy(format_all)
    in_format.format_common_info.precision = 24
    out_format: db.FormatInfo = deepcopy(format_none)
    out_format.format_common_info.precision = 18

    assert db.calc_conversion_weight(converter_ob, in_format,
                                     out_format,) == ((max_prop_weight << db.PROP_WEIGHT_BIT_OFFSET) +
                                                      (1 << 6*db.PREC_GAP_BITS << db.PREC_WEIGHT_BIT_OFFSET) +
                                                      (converter_ob.weight << db.CONV_WEIGHT_BIT_OFFSET))


@pytest.mark.parametrize("prop_weight, prec_weight, time_weight, conv_weight", [(0, 0, 0, 0),
                                                                                (65535, 65535, 255, 255),
                                                                                (2788794, 1254542, 122, 234)])
def test_split_conversion_weight(prop_weight, prec_weight, time_weight, conv_weight, subtests):
    """Test that the function to split the conversion weight works as expected"""
    split_weight = db.split_conversion_weight(db.combine_conversion_weight(
        prop_weight, prec_weight, time_weight, conv_weight))
    with subtests.test("Property weight is correctly split out from full weight"):
        assert split_weight.prop_weight == prop_weight
    with subtests.test("Precision weight is correctly split out from full weight"):
        assert split_weight.prec_weight == prec_weight
    with subtests.test("Time weight is correctly split out from full weight"):
        assert split_weight.time_weight == time_weight
    with subtests.test("Converter weight is correctly split out from full weight"):
        assert split_weight.conv_weight == conv_weight
