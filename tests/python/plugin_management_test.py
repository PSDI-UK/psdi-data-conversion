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

PROJECT_PATH = os.path.realpath(os.path.join(utils.__file__, "../.."))


@pytest.fixture(autouse=True, scope="module")
def check_editable_mode():
    """An autouse fixture for all tests in this module which causes them to fail if the package isn't installed in
    "editable" mode.
    """
    utils.confirm_editable_mode()


@pytest.fixture()
def mock_repo(tmp_path_factory):
    """A fixture that provides a temporary copy of the project repo, ignoring various files and folders that are
    unneeded for this test
    """
    tmp_path: Path = tmp_path_factory.mktemp("mock-repo")
    mock_repo_path = os.path.join(tmp_path, "psdi-data-conversion")
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
    conv_path: str | None = None
    _plugin_path: str | None = None
    conv_module_path: str | None = None
    data_path: str | None = None
    _conv_module_text: str | None = None
    _conv_data: utils.JsonMainDict | None = None

    @pytest.fixture(autouse=True)
    def _setup(self, mock_repo):
        """Setup for each test"""
        self.mock_repo = mock_repo
        self.conv_path = os.path.join(self.mock_repo, "psdi_data_conversion", "converters")

    def _reset_plugin_info(self):
        self._conv_module_text = None
        self._conv_data = None

    @property
    def plugin_path(self):
        return self._plugin_path

    @plugin_path.setter
    def plugin_path(self, val: str):
        self._plugin_path = val
        self.conv_module_path = os.path.join(self._plugin_path, "converter.py")
        self.data_path = os.path.join(self._plugin_path, "data.json")

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
    PLUGIN_NAME_3 = "Test Plugin Name 3"
    PLUGIN_NAME_4 = "Test Plugin Name 4"
    SCRIPT_PLUGIN_NAME = "Test Script Plugin Name"
    SCRIPT_PLUGIN_NAME_2 = "Test Script Plugin Name 2"
    UNCONVERTABLE_NAME = "++"

    # Test label we'll use for the plugin, differing from that generated from the name
    TEST_PLUGIN_LABEL = "test_plugin_label"
    TEST_PLUGIN_LABEL_2 = "test_plugin_label_2"
    TEST_SCRIPT_PLUGIN_LABEL = "test_script_plugin_label"
    INVALID_LABEL_CAPS = "Label"
    INVALID_LABEL_CHAR = "label++"

    def _run_create_plugin(self, name: str, label: str | None = None, script=False, expect_fail=False, check=True):
        """Calls the script to create a plugin"""

        l_args = ["psdi-data-convert-create-plugin", name]
        if label:
            l_args += ["--label", label]
        if script:
            l_args.append("--script")
        env = {**os.environ,
               TEST_PATH_KEY: self.mock_repo}

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

        self.plugin_path: str = os.path.join(self.conv_path, label)
        assert os.path.isdir(self.plugin_path)
        assert os.path.isfile(os.path.join(self.plugin_path, "__init__.py"))
        assert os.path.isfile(self.conv_module_path)
        assert os.path.isfile(self.data_path)

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

    def test_create_plugin_simple(self):
        """Test that the plugin creation script works as expected in simple cases"""

        # Test a simple case where just a name is supplied
        self._run_create_plugin(self.PLUGIN_NAME)

        # Test a case where we provide a separate label that isn't generated from the name
        self._run_create_plugin(self.PLUGIN_NAME_2, label=self.TEST_PLUGIN_LABEL)

        # Then repeat these with a script plugin
        self._run_create_plugin(self.SCRIPT_PLUGIN_NAME, script=True)
        self._run_create_plugin(self.SCRIPT_PLUGIN_NAME_2, label=self.TEST_SCRIPT_PLUGIN_LABEL, script=True)

    def test_create_plugin_fail_cases(self):
        """Test that the plugin creation script fails when expected to"""

        # Test that we get a failure when there's a name clash
        self._run_create_plugin(self.PLUGIN_NAME, check=False)
        process = self._run_create_plugin(self.PLUGIN_NAME, label=self.TEST_PLUGIN_LABEL, expect_fail=True)
        assert f"Name '{self.PLUGIN_NAME}' clashes" in process.stderr

        # Test that we get a failure when there's a label clash
        self._run_create_plugin(self.PLUGIN_NAME_2, label=self.TEST_PLUGIN_LABEL, check=False)
        process = self._run_create_plugin(self.PLUGIN_NAME_3, label=self.TEST_PLUGIN_LABEL, expect_fail=True)
        assert f"Label '{self.TEST_PLUGIN_LABEL}' clashes"

        # Test we get a failure if the name can't be converted to a valid label, but not if we provide a valid label
        process = self._run_create_plugin(self.UNCONVERTABLE_NAME, expect_fail=True)
        assert f"A valid label could not be generated from converter name '{self.UNCONVERTABLE_NAME}'" in process.stderr
        self._run_create_plugin(self.UNCONVERTABLE_NAME, label=self.TEST_PLUGIN_LABEL_2)

        # Test we get a failure if the provided label is invalid
        process = self._run_create_plugin(self.PLUGIN_NAME_4, label=self.INVALID_LABEL_CAPS, expect_fail=True)
        assert f"Label '{self.INVALID_LABEL_CAPS}' is invalid" in process.stderr
        process = self._run_create_plugin(self.PLUGIN_NAME_4, label=self.INVALID_LABEL_CHAR, expect_fail=True)
        assert f"Label '{self.INVALID_LABEL_CHAR}' is invalid" in process.stderr


class TestInstallPlugins(PluginManagementBase):

    def _create_test_plugin(self, which: str):
        """Calls the plugin creation script in test mode """

        l_args: list[str] = ["psdi-data-convert-create-plugin", f"test_converter_{which}"]
        env = {**os.environ,
               TEST_PATH_KEY: self.mock_repo,
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

        env = {**os.environ, TEST_PATH_KEY: self.mock_repo}

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
        self.plugin_path: str = os.path.join(self.conv_path, f"test_converter_{which}")
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

    def test_questionable_install(self):
        """Test installing a plugin where it's questionable whether one of the formats in it might already be in the
        database"""
        self._create_test_plugin("questionable")

        # Save a copy of the database from before installing to ensure it isn't changed
        self.plugin_path = os.path.join(self.conv_path, "test_converter_questionable")
        conv_data_init = self.conv_data

        process = self._run_install_plugins(expect_fail=True)
        assert process.returncode == 2

        # Check that the converter's data hasn't been changed
        self.plugin_path = os.path.join(self.conv_path, "test_converter_questionable")
        conv_data_post_install = self.conv_data

        assert conv_data_init == conv_data_post_install

        # Now try setting the questionable format as confirmed new, and check that we can install after doing so
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
