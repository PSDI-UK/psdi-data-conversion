"""
# converter_test.py

Unit tests of the converter class. This module uses the common test specifications defined in
psdi_data_conversion/testing/conversion_test_specs.py so that a common set of conversion tests is performed through
the Python library (this module), the command-line application, and the GUI.
"""

import logging
import math
import os

import pytest

from psdi_data_conversion import constants as const
from psdi_data_conversion.converter import L_REGISTERED_CONVERTERS
from psdi_data_conversion.converters.c2x import C2xFileConverter
from psdi_data_conversion.converters.openbabel import OBFileConverter
from psdi_data_conversion.database import get_format_info
from psdi_data_conversion.testing.conversion_test_specs import l_library_test_specs
from psdi_data_conversion.testing.utils import run_test_conversion_with_library
from psdi_data_conversion.utils import regularize_name


@pytest.fixture(autouse=True)
def setup_test() -> None:
    """Reset global aspects before a test, so that different tests won't interfere with each other"""

    # Remove the global log file if one exists
    try:
        os.remove(const.GLOBAL_LOG_FILENAME)
    except FileNotFoundError:
        pass

    # Clear any existing loggers so new ones will be created fresh
    logging.Logger.manager.loggerDict.clear()


def test_default():
    """Test that the default converter is registered.
    """
    assert regularize_name(const.CONVERTER_DEFAULT) in L_REGISTERED_CONVERTERS


@pytest.mark.parametrize("test_spec", l_library_test_specs,
                         ids=lambda x: x.name)
def test_conversions(test_spec):
    """Run all conversion tests in the defined list of test specifications
    """
    run_test_conversion_with_library(test_spec)


def test_envvars():
    """Test that setting appropriate envvars will set them for a file converter
    """

    test_file_size = 1234
    os.environ[const.MAX_FILESIZE_EV] = str(test_file_size)

    pdb_format_id = get_format_info("pdb", which=0).id

    converter = C2xFileConverter(filename="1NE6.mmcif",
                                 to_format=pdb_format_id,
                                 use_envvars=True,)
    assert math.isclose(converter.max_file_size, test_file_size*const.MEGABYTE)

    # And also check it isn't applied if we don't ask it to use envvars
    converter_no_ev = C2xFileConverter(filename="1NE6.mmcif",
                                       to_format=pdb_format_id,
                                       use_envvars=False,)
    assert not math.isclose(converter_no_ev.max_file_size, test_file_size*const.MEGABYTE)

    # And check that OB uses its own EV
    converter = OBFileConverter(filename="1NE6.mmcif",
                                to_format=pdb_format_id,
                                use_envvars=True,)
    assert not math.isclose(converter.max_file_size, test_file_size*const.MEGABYTE)

    os.environ[const.MAX_FILESIZE_OB_EV] = str(test_file_size)
    converter = OBFileConverter(filename="1NE6.mmcif",
                                to_format=pdb_format_id,
                                use_envvars=True,)
    assert math.isclose(converter.max_file_size, test_file_size*const.MEGABYTE)
