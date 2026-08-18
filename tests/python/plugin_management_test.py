"""@file tests/plugin_management_test.py

Created 2026-08-07 by Bryan Gillis.

Tests of the scripts to create and install converter plugins
"""

import json
import os
import shutil
import subprocess
from pathlib import Path

import pytest

from psdi_data_conversion import database as db
from psdi_data_conversion import utils
from psdi_data_conversion.main_create_plugin import TEST_DATA_KEY
from psdi_data_conversion.main_install_plugins import THRESHOLD_FORMAT_ID
from psdi_data_conversion.testing.constants import TEST_PATH_KEY

PROJECT_PATH = (Path(utils.__file__) / "../..").resolve()


if not utils.in_editable_mode():
    pytest.skip(reason="Plugin management scripts can only be tested when the project is installed in --editable mode",
                allow_module_level=True)


@pytest.fixture()
def mock_repo(tmp_path_factory) -> Path:
    """A fixture that provides a temporary copy of the project repo, ignoring various files and folders that are
    unneeded for this test
    """
    tmp_path: Path = tmp_path_factory.mktemp("mock-repo")
    mock_repo_path: Path = tmp_path / "psdi-data-conversion"
    shutil.copytree(PROJECT_PATH, mock_repo_path, symlinks=True,
                    ignore=shutil.ignore_patterns(".*", "bin", "gui", "templates", "scripts", "tests", "doc", "html",
                                                  "*.md", "*.html", "*.toml", "Dockerfile", "LICENSE", "__pycache__",
                                                  "testing", "output", "*.mmcif", "*.cif", "*.ins", "*.mol",
                                                  "*.molden", "*.tar", "*.tar.gz", "*.zip", "*.inchi", "*.cub", "*.xyz",
                                                  "*.pdb", "*.outmol", "*.cdxml", "*.cdjson", "*:Zone.Identifier"))
    return mock_repo_path


class PluginManagementBase:

    # Values set during execution
    mock_repo: str | None = None
    conv_path: Path | None = None
    _plugin_path: Path | None = None
    conv_module_path: Path | None = None
    data_path: Path | None = None
    _conv_module_text: str | None = None
    _conv_data: utils.JsonMainDict | None = None

    @pytest.fixture(autouse=True)
    def _setup(self, mock_repo):
        """Setup for each test"""
        self.mock_repo = mock_repo
        self.conv_path = self.mock_repo / "psdi_data_conversion/converters"

    def _reset_plugin_info(self):
        self._conv_module_text = None
        self._conv_data = None

    @property
    def plugin_path(self) -> Path:
        return self._plugin_path

    @plugin_path.setter
    def plugin_path(self, val: Path):
        self._plugin_path = val
        self.conv_module_path = self._plugin_path / "converter.py"
        self.data_path = self._plugin_path / "data.json"

        # Reset private variables that depend on this
        self._reset_plugin_info()

    @property
    def conv_module_text(self) -> str:
        if not self._conv_module_text:
            self._conv_module_text = open(self.conv_module_path).read()
        return self._conv_module_text

    @property
    def conv_data(self) -> utils.JsonMainDict:
        if not self._conv_data:
            self._conv_data = json.load(open(self.data_path))
        return self._conv_data


class TestCreatePlugin(PluginManagementBase):

    # Constants

    # Test names we'll use for the plugins
    PLUGIN_NAME = "Test Plugin Name"
    PLUGIN_NAME_2 = "Test Plugin Name 2"
    UNCONVERTABLE_NAME = "++"

    # Test label we'll use for the plugin, differing from that generated from the name
    TEST_PLUGIN_LABEL = "test_plugin_label"
    INVALID_LABEL_CAPS = "Label"
    INVALID_LABEL_CHAR = "label++"

    def _run_create_plugin(self, name: str, label: str | None = None, script=False, expect_fail=False, check=True):
        """Calls the script to create a plugin"""

        l_args = ["psdi-data-convert-create-plugin", name]
        if label:
            l_args += ["--label", label]
        if script:
            l_args.append("--script")
        env = {**os.environ, TEST_PATH_KEY: str(self.mock_repo)}

        process = subprocess.run(l_args, capture_output=True, text=True, env=env)

        if not expect_fail and process.returncode:
            pytest.fail(f"Plugin creation failed with return code {process.returncode} and stderr:\n{process.stderr}")
        elif expect_fail and not process.returncode:
            pytest.fail(f"Plugin creation succeeded when failure was expected with stdout:\n{process.stdout}")

        if check and not expect_fail:
            self._check_plugin(name, label=label, script=script)

        return process

    def _check_plugin(self, name: str, label: str | None = None, script=False):

        if label is None:
            label = "_".join(name.split()).lower()

        self.plugin_path: Path = self.conv_path / label
        assert self.plugin_path.is_dir()
        assert (self.plugin_path / "__init__.py").is_file()
        assert self.conv_module_path.is_file()
        assert self.data_path.is_file()

        # Check that the converter module was created as expected
        pascal_name = "".join([x.capitalize() for x in name.split(" ")]).replace(" ", "")
        assert f"{pascal_name} file converter" in self.conv_module_text
        if not script:
            assert f"class {pascal_name}FileConverter(FileConverter)" in self.conv_module_text
        else:
            assert f"class {pascal_name}FileConverter(ScriptFileConverter)" in self.conv_module_text
        assert f"File converter specialised to use {pascal_name} for conversions" in self.conv_module_text
        assert f"converter = {pascal_name}FileConverter" in self.conv_module_text

        # Check that the data was created as expected
        conv_meta: utils.JsonDict = self.conv_data[db.DB_CONVERTER_KEY]
        assert conv_meta[db.DB_ID_KEY] is None
        assert conv_meta[db.DB_NAME_KEY] == name
        assert conv_meta[db.DB_DESC_KEY] == f"{name} converter plugin"
        assert conv_meta[db.DB_INFO_KEY] == ""
        assert conv_meta[db.DB_URL_KEY] == ""
        assert conv_meta[db.DB_SUPPORT_AMBIG_EXT_KEY] is None
        assert conv_meta[db.DB_KEY_PREFIX_KEY] is None

        for key in (db.DB_EXTRA_FORMATS_KEY, db.DB_SUPPORTED_FORMATS_KEY, db.DB_IN_ONLY_FORMATS_KEY,
                    db.DB_OUT_ONLY_FORMATS_KEY, db.DB_SUPPORTED_CONVERSIONS_KEY, db.DB_UNSUPPORTED_CONVERSIONS_KEY,
                    db.DB_IN_FLAGS_KEY, db.DB_OUT_FLAGS_KEY, db.DB_IN_OPTIONS_KEY, db.DB_OUT_OPTIONS_KEY):
            val: list = self.conv_data[key]
            assert isinstance(self.conv_data[key], list)
            assert len(val) == 0

    def test_create_plugin_name_only(self):
        """Test creating a plugin in a simple case where just a name is supplied"""
        self._run_create_plugin(self.PLUGIN_NAME)

    def test_create_plugin_with_label(self):
        """Test a case where we provide a separate label that isn't generated from the name"""
        self._run_create_plugin(self.PLUGIN_NAME, label=self.TEST_PLUGIN_LABEL)

    def test_create_script_plugin_name_only(self):
        """Test creating a plugin in a simple case where just a name is supplied"""
        self._run_create_plugin(self.PLUGIN_NAME, script=True)

    def test_create_script_plugin_with_label(self):
        """Test a case where we provide a separate label that isn't generated from the name"""
        self._run_create_plugin(self.PLUGIN_NAME, label=self.TEST_PLUGIN_LABEL, script=True)

    def test_create_plugin_name_clash(self):
        """Test that we get a failure when there's a name clash"""
        self._run_create_plugin(self.PLUGIN_NAME, check=False)
        process = self._run_create_plugin(self.PLUGIN_NAME, label=self.TEST_PLUGIN_LABEL, expect_fail=True)
        assert f"Name '{self.PLUGIN_NAME}' clashes" in process.stderr

    def test_create_plugin_label_clash(self):
        """Test that we get a failure when there's a label clash"""
        self._run_create_plugin(self.PLUGIN_NAME, label=self.TEST_PLUGIN_LABEL, check=False)
        process = self._run_create_plugin(self.PLUGIN_NAME_2, label=self.TEST_PLUGIN_LABEL, expect_fail=True)
        assert f"Label '{self.TEST_PLUGIN_LABEL}' clashes" in process.stderr

    def test_create_plugin_unconvertable_name(self):
        """Test we get a failure if the name can't be converted to a valid label"""
        process = self._run_create_plugin(self.UNCONVERTABLE_NAME, expect_fail=True)
        assert f"A valid label could not be generated from converter name '{self.UNCONVERTABLE_NAME}'" in process.stderr

    def test_create_plugin_unconvertable_name_with_label(self):
        """Test that we can use an unconvertable name if we also provide a valid label"""
        self._run_create_plugin(self.UNCONVERTABLE_NAME, label=self.TEST_PLUGIN_LABEL)

    def test_create_plugin_invalid_label_caps(self):
        """Test we get a failure if the provided label is invalid due to including capital letters"""
        process = self._run_create_plugin(self.PLUGIN_NAME, label=self.INVALID_LABEL_CAPS, expect_fail=True)
        assert f"Label '{self.INVALID_LABEL_CAPS}' is invalid" in process.stderr

    def test_create_plugin_invalid_label_chars(self):
        """Test we get a failure if the provided label is invalid due to including invalid characters"""
        process = self._run_create_plugin(self.PLUGIN_NAME, label=self.INVALID_LABEL_CHAR, expect_fail=True)
        assert f"Label '{self.INVALID_LABEL_CHAR}' is invalid" in process.stderr


class TestInstallPlugins(PluginManagementBase):

    def _create_test_plugin(self, which: str):
        """Calls the plugin creation script in test mode """

        l_args: list[str] = ["psdi-data-convert-create-plugin", f"test_converter_{which}"]
        env = {**os.environ,
               TEST_PATH_KEY: str(self.mock_repo),
               TEST_DATA_KEY: "True"}

        process = subprocess.run(l_args, capture_output=True, text=True, env=env)
        if process.returncode:
            pytest.fail(
                f"Test plugin creation failed with return code {process.returncode} and stderr:\n{process.stderr}")

        return process

    def _run_install_plugins(self, force=False, expect_fail=False):
        """Calls the script to install plugins"""

        l_args: list[str] = ["psdi-data-convert-install-plugins"]

        if force:
            l_args.append("-f")

        env = {**os.environ, TEST_PATH_KEY: str(self.mock_repo)}

        process = subprocess.run(l_args, capture_output=True, text=True, env=env)

        if not expect_fail and process.returncode:
            pytest.fail(f"Plugin installation failed with return code {process.returncode}, stdout:\n{process.stdout}\n"
                        f"and stderr:\n{process.stderr}")
        elif expect_fail and not process.returncode:
            pytest.fail(f"Plugin installation succeeded when failure was expected with stdout:\n{process.stdout}")

        return process

    def _check_plugin_installed(self, which: str):
        """Checks that a test plugin has been successfully installed"""

        # Set the name of the plugin and get its JSON data
        self.plugin_path: Path = self.conv_path / f"test_converter_{which}"
        conv_data = self.conv_data

        # Look through the converter's data to check that it all appears to be installed
        conv_meta: utils.JsonDict = conv_data[db.DB_CONVERTER_KEY]
        assert conv_meta[db.DB_ID_KEY] > THRESHOLD_FORMAT_ID
        assert conv_meta[db.DB_SUPPORT_AMBIG_EXT_KEY] is not None

        # Check that no extra formats remain in the converter's data
        assert len(conv_data[db.DB_EXTRA_FORMATS_KEY]) == 0

        # Check that all referenced format IDs have been converted to UUIDs
        for key in (db.DB_SUPPORTED_FORMATS_KEY, db.DB_IN_ONLY_FORMATS_KEY, db.DB_OUT_ONLY_FORMATS_KEY):
            for format_id in conv_data[key]:
                assert format_id > THRESHOLD_FORMAT_ID

        for key in (db.DB_SUPPORTED_CONVERSIONS_KEY, db.DB_UNSUPPORTED_CONVERSIONS_KEY):
            for d_conversion in conv_data[key]:
                assert d_conversion[db.DB_IN_ID_KEY] > THRESHOLD_FORMAT_ID
                assert d_conversion[db.DB_OUT_ID_KEY] > THRESHOLD_FORMAT_ID

        for key in (db.DB_IN_FLAGS_KEY, db.DB_OUT_FLAGS_KEY, db.DB_IN_OPTIONS_KEY, db.DB_OUT_OPTIONS_KEY):
            for d_flag in conv_data[key]:
                for format_id in d_flag[db.DB_FORMAT_ID_LIST_KEY]:
                    assert format_id > THRESHOLD_FORMAT_ID

        # Check that all information about the converter is now in the proper database file
        os.environ[TEST_PATH_KEY] = str(self.mock_repo)
        database = db.load_database(prune=False)
        conv_info: db.ConverterInfo = database.get_converter_info(conv_meta[db.DB_ID_KEY])

        assert conv_info.id == conv_meta[db.DB_ID_KEY]
        assert conv_info.pretty_name == conv_meta[db.DB_NAME_KEY]
        assert conv_info.description == conv_meta[db.DB_DESC_KEY]
        assert conv_info.url == conv_meta[db.DB_URL_KEY]
        assert conv_info.supported is False
        assert conv_info.registered is False

    def test_simple_install(self):
        """Test installing a simple plugin"""
        self._create_test_plugin("simple")
        self._run_install_plugins()
        self._check_plugin_installed("simple")

    def test_complex_install(self):
        """Test installing a complex plugin"""
        self._create_test_plugin("complex")
        self._run_install_plugins()
        self._check_plugin_installed("complex")

    def test_questionable_install_fail(self):
        """Test that an installation with a questionable format fails if it isn't confirmed as new"""
        self._create_test_plugin("questionable")

        # Save a copy of the database from before installing to ensure it isn't changed
        self.plugin_path = self.conv_path / "test_converter_questionable"
        conv_data_init = self.conv_data

        process = self._run_install_plugins(expect_fail=True)
        assert process.returncode == 2

        # Check that the converter's data hasn't been changed

        # Set the plugin path again to force a refresh of cached data
        self.plugin_path = self.conv_path / "test_converter_questionable"
        conv_data_post_install = self.conv_data

        assert conv_data_init == conv_data_post_install

    def test_questionable_install_success(self):
        """Test that an installation with a questionable format succeeds if it is confirmed as new"""
        self._create_test_plugin("questionable")
        self.plugin_path = self.conv_path / "test_converter_questionable"

        # Set the questionable format as confirmed new, and check that we can install after doing so
        format_data: utils.JsonDict = self.conv_data[db.DB_EXTRA_FORMATS_KEY][4]

        # Check that we have the right format, in case the index of it changes, so we catch it explictly here
        assert format_data[db.DB_FORMAT_EXT_KEY] == "mol"

        # Set it as confirmed new and save it
        format_data[db.DB_FORMAT_CONFIRMED_NEW_KEY] = True
        json.dump(self.conv_data, open(self.data_path, "w"))

        # And check that it works if we run again
        self._run_install_plugins()
        self._check_plugin_installed("questionable")

    def test_questionable_install_force(self):
        """Test installing a plugin where it's questionable whether one of the formats in it might already be in the
        database"""
        self._create_test_plugin("questionable")
        self._run_install_plugins(force=True)
        self._check_plugin_installed("questionable")

    def test_conv_weight_calcs(self):
        """Test calculation of conversion weights"""
        self._create_test_plugin("complex")
        self._run_install_plugins()

        os.environ[TEST_PATH_KEY] = str(self.mock_repo)
        database = db.load_database(prune=False)
        l_d_converts_to: list[utils.JsonDict] = database.converts_to
        for d_converts_to in l_d_converts_to:
            in_format = database.get_format_info(d_converts_to[db.DB_IN_ID_KEY])
            out_format = database.get_format_info(d_converts_to[db.DB_OUT_ID_KEY])
            assert d_converts_to[db.DB_WEIGHT_KEY] == db.combine_conversion_weight(in_format, out_format)
