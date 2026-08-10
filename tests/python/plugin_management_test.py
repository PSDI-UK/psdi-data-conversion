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


class TestCreatePlugin:

    # Constants

    # Test name we'll use for the plugin
    PLUGIN_NAME = "Test Plugin Name"

    # Expected label generated for the plugin based on the above name
    EX_PLUGIN_LABEL = "test_plugin_name"

    # Test label we'll use for the plugin, differing from that generated from the name
    TEST_PLUGIN_LABEL = "test_plugin_label"

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

    @property
    def plugin_path(self):
        return self._plugin_path

    @plugin_path.setter
    def plugin_path(self, val: str):
        self._plugin_path = val
        self.conv_module_path = os.path.join(self._plugin_path, "converter.py")
        self.data_path = os.path.join(self._plugin_path, "data.json")

        # Reset private variables that depend on this
        self._conv_module_text = None
        self._conv_data = None

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

    def _run_create_plugin(self, name: str, label: str | None = None, script=False):
        """Calls the script to create a plugin"""

        l_args = ["psdi-data-convert-create-plugin", name, "--test-path", self.mock_repo]
        if label:
            l_args += ["--label", label]
        if script:
            l_args.append("--script")

        return subprocess.run(l_args, capture_output=True, text=True)

    def test_create_plugin(self):
        """Test that the plugin creation script works as expected"""
        process = self._run_create_plugin(self.PLUGIN_NAME)

        # Check there was no error in execution
        if process.returncode:
            pytest.fail(f"Plugin creation failed with return code {process.returncode} and stderr:\n{process.stderr}")

        # Check that the plugin was created as expected
        self.plugin_path: str = os.path.join(self.conv_path, self.EX_PLUGIN_LABEL)
        assert os.path.isdir(self.plugin_path)
        assert os.path.isfile(os.path.join(self.plugin_path, "__init__.py"))
        assert os.path.isfile(self.conv_module_path)
        assert os.path.isfile(self.data_path)

        # Check that the converter module was created as expected
        assert f"{self.PLUGIN_NAME.replace(" ", "")} file converter" in self.conv_module_text
        assert f"class {self.PLUGIN_NAME.replace(" ", "")}FileConverter(FileConverter)" in self.conv_module_text
        assert (f"File converter specialised to use {self.PLUGIN_NAME.replace(" ", "")} for conversions"
                in self.conv_module_text)
        assert f"converter = {self.PLUGIN_NAME.replace(" ", "")}FileConverter" in self.conv_module_text

        # Check that the data was created as expected
        conv_meta: utils.JsonDict = self.conv_data[db.DB_CONVERTER_KEY]
        assert conv_meta[db.DB_ID_KEY] is None
        assert conv_meta[db.DB_NAME_KEY] == self.PLUGIN_NAME
        assert conv_meta[db.DB_DESC_KEY] == f"{self.PLUGIN_NAME} converter plugin"
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
