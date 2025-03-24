"""@file psdi-data-conversion/tests/converter_test.py

Created 2024-12-15 by Bryan Gillis.

Unit tests of the converter class
"""

from copy import deepcopy
import logging
import math
import os
import shutil
from typing import Any
import pytest

from psdi_data_conversion import constants as const
from psdi_data_conversion.converter import L_REGISTERED_CONVERTERS
from psdi_data_conversion.converters.openbabel import CONVERTER_OB, OBFileConverter
from psdi_data_conversion.main import FileConverterInputException
from psdi_data_conversion.testing.utils import run_test_conversion_with_library
from psdi_data_conversion.testing import conversion_test_specs as specs


@pytest.fixture()
def base_data():
    """A fixture providing a default `data` object which can be used to instantiate a converter
    """
    return {'token': '1041c0a661d118d5f28e7c6830375dd0',
            'from_full': 'mmcif: Macromolecular Crystallographic Info',
            'to_full': 'pdb: Protein Data Bank',
            'success': 'good',
            'from_flags': '',
            'to_flags': '',
            'from_arg_flags': '',
            'from_args': '',
            'to_arg_flags': '',
            'to_args': '',
            'coordinates': 'neither',
            'coordOption': 'medium',
            'upload_file': 'true'}


@pytest.fixture()
def tmp_upload_path(tmp_path):
    """Fixture providing a temporary directory for upload files
    """
    upload_path = os.path.join(tmp_path, "uploads")
    os.makedirs(upload_path, exist_ok=True)
    return upload_path


@pytest.fixture()
def tmp_download_path(tmp_path):
    """Fixture providing a temporary directory for download files
    """
    download_path = os.path.join(tmp_path, "downloads")
    os.makedirs(download_path, exist_ok=True)
    return download_path


def test_default():
    """Test that the default converter is registered.
    """
    assert const.CONVERTER_DEFAULT in L_REGISTERED_CONVERTERS


class TestConverter:

    @pytest.fixture(autouse=True)
    def setup_test(self, base_data: dict[str, Any]) -> None:
        """Reset global aspects before a test, so that different tests won't interfere with each other,
        and save references to fixtures.
        """

        # Remove the global log file if one exists
        try:
            os.remove(const.GLOBAL_LOG_FILENAME)
        except FileNotFoundError:
            pass

        # Clear any existing loggers so new ones will be created fresh
        logging.Logger.manager.loggerDict.clear()

        # Save fixtures
        self.base_data = base_data

    def get_converter_kwargs(self,
                             data: dict | None = None,
                             **kwargs):
        """Get the keyword arguments to be passed to a FileConverter for testing
        """
        combined_data = deepcopy(self.base_data)
        if data:
            combined_data.update(data)
        standard_kwargs: dict[str, Any] = {"data": combined_data}
        standard_kwargs.update(kwargs)
        return standard_kwargs

    def check_file_status(self, input_exist=None, output_exist=None):
        """Common check for unit tests on whether the input/output files from a conversion exist or not

        Parameters
        ----------
        input_exist : bool | None
            If True, will check that the input file does exist, if False, will check that it doesn't, if None, won't
            check either way.
        output_exist : bool | None
            If True, will check that the output file does exist, if False, will check that it doesn't, if None, won't
            check either way.
        """

        ex_output_filename_base = os.path.splitext(self.filename_base)[0]
        ex_output_filename = os.path.join(self.tmp_download_path, f"{ex_output_filename_base}.{self.to_format}")

        for check_condition, filename in ((input_exist, self.source_filename),
                                          (output_exist, ex_output_filename)):
            if check_condition is None:
                continue
            if check_condition:
                assert os.path.isfile(filename), f"Expected file {filename} does not exist"
            else:
                assert not os.path.exists(filename), f"File {filename} exists, but should have been deleted"

    def get_logs(self):
        """Get the log filenames and text content after the converter has run, for each of the three log types
        """

        self.global_log_filename = const.GLOBAL_LOG_FILENAME
        self.output_log_filename = os.path.join(self.tmp_download_path,
                                                f"{self.filename_base}{const.OUTPUT_LOG_EXT}")

        self.global_log_text: str | None = None
        self.output_log_text: str | None = None

        for log_type in ("global", "output"):
            try:
                log_file = getattr(self, f"{log_type}_log_filename")
                log_text = open(log_file, "r").read()
            except FileNotFoundError:
                log_text = ""
            setattr(self, f"{log_type}_log_text", log_text)

    def test_basic_conversions(self):
        """Run a basic set of conversions with various converters and file formats which we expect to succeed without
        issue.
        """
        run_test_conversion_with_library(specs.basic_tests)

    def test_open_babel_warning(self):
        """Run a test that expected warnings from Open Babel are captured in the log
        """
        run_test_conversion_with_library(specs.open_babel_warning_test)

    def test_exceed_output_file_size(self):
        """Run a test of the converter to ensure it reports an error properly if the output file size is too large
        """
        run_test_conversion_with_library(specs.max_size_test)

    def test_invalid_converter(self):
        """Run a test of the converter to ensure it reports an error properly if an invalid converter is requested
        """

        self.get_input_info(filename="1NE6.mmcif")

        self.run_converter(expect_exception=FileConverterInputException,
                           name="INVALID")

        # Check that the input file remains but no output file has been created
        self.check_file_status(input_exist=True, output_exist=False)

    def test_quality_note(self):
        """Run a test of the converter on an `.xyz` to `.inchi` conversion which we expect to have warnings about data
        loss and extrapolation
        """
        run_test_conversion_with_library(specs.quality_note_test)

    def test_cleanup(self):
        """Test that input files are deleted if requested
        """

        self.get_input_info(filename="hemoglobin.pdb",
                            to_format="cif")

        # Make a copy of the source file in the uploads directory, and point the self.source_filename variable to it
        # so it'll be properly checked to be deleted later
        upload_filename = os.path.join(self.tmp_upload_path, self.local_filename)
        shutil.copyfile(self.source_filename, upload_filename)
        self.source_filename = upload_filename

        self.run_converter(name=CONVERTER_OB,
                           filename=upload_filename,
                           delete_input=True)

        # Check that the input file has been deleted and the output file exists where we expect it to
        self.check_file_status(input_exist=False, output_exist=True)

    def test_xyz_to_inchi_err(self):
        """Run a test of the converter on an `.xyz` to `.inchi` conversion we expect to fail
        """

        self.get_input_info(filename="quartz_err.xyz",
                            to_format="inchi")

        self.run_converter(expect_code=const.STATUS_CODE_GENERAL)

        # Check that the input and output files have properly been deleted
        self.check_file_status(input_exist=True, output_exist=False)

    def test_envvars(self):
        """Test that setting appropriate envvars will set them for a file converter
        """

        test_file_size = 1234
        os.environ[const.MAX_FILESIZE_EV] = str(test_file_size)

        self.get_input_info(filename="1NE6.mmcif")
        converter = OBFileConverter(use_envvars=True,
                                    **self.get_converter_kwargs())
        assert math.isclose(converter.max_file_size, test_file_size*const.MEGABYTE)

        # And also check it isn't applied if we don't ask it to use envvars
        converter_no_ev = OBFileConverter(use_envvars=False,
                                          **self.get_converter_kwargs())
        assert not math.isclose(converter_no_ev.max_file_size, test_file_size*const.MEGABYTE)
