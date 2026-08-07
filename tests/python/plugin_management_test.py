"""@file tests/plugin_management_test.py

Created 2026-08-07 by Bryan Gillis.

Tests of the scripts to create and install converter plugins
"""

import os
import shutil

import pytest

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
    """A fixture that provides a temporary copy of the project repo
    """
    with tmp_path_factory.mktemp("mock-repo") as tmp_path:
        mock_repo_path = os.path.join(tmp_path, "psdi-data-conversion")
        shutil.copytree(PROJECT_PATH, mock_repo_path)
        yield mock_repo_path
