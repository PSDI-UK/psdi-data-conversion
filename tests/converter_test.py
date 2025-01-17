"""@file psdi-data-conversion/tests/converter_test.py

Created 2024-12-15 by Bryan Gillis.

Unit tests of the converter class
"""

from copy import deepcopy
import logging
import os
import re
import pytest

from psdi_data_conversion.log_utility import DATETIME_RE_RAW, GLOBAL_LOG_FILENAME
from psdi_data_conversion.converter import (CONVERTER_ATO, CONVERTER_OB, get_file_storage, LOCAL_LOG_EXT,
                                            OUTPUT_LOG_EXT, FILE_TO_UPLOAD_KEY, STATUS_CODE_BAD_METHOD,
                                            STATUS_CODE_GENERAL, STATUS_CODE_SIZE, FileConverter,
                                            FileConverterAbortException)


@pytest.fixture()
def base_mock_form():
    """A fixture providing a default `form` object which can be used to instantiate a converter
    """
    return {'token': '1041c0a661d118d5f28e7c6830375dd0',
            'converter': CONVERTER_OB,
            'from': 'mmcif',
            'to': 'pdb',
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


class TestConverter:

    @pytest.fixture(autouse=True)
    def setup_test(self, base_mock_form, tmp_upload_path, tmp_download_path, test_data_loc):
        """Reset global aspects before a test, so that different tests won't interfere with each other,
        and save references to fixtures.
        """

        # Remove the global log file if one exists
        try:
            os.remove(GLOBAL_LOG_FILENAME)
        except FileNotFoundError:
            pass

        # Clear any existing loggers so new ones will be created fresh
        logging.Logger.manager.loggerDict.clear()

        # Save test data location
        self.test_data_loc = test_data_loc

        # Save tmp directories
        self.base_mock_form = base_mock_form
        self.tmp_upload_path = tmp_upload_path
        self.tmp_download_path = tmp_download_path

    def get_input_info(self, filename: str, **kwargs,):
        """Sets up a mock form for input and gets various variables we'll want to use for checks on output

        Parameters
        ----------
        filename : str
            The name of the file to use as input for the test
        """

        self.mock_form = deepcopy(self.base_mock_form)

        for key, value in kwargs.items():
            if key in self.mock_form:
                self.mock_form[key] = value
            else:
                raise RuntimeError(f"Invalid key {key} provided for form")

        # Save some variables from input we'll be using throughout this test
        self.source_filename = os.path.join(self.test_data_loc, filename)
        self.files = get_file_storage(self.source_filename)
        self.filename = self.files[FILE_TO_UPLOAD_KEY].filename
        self.filename_base = os.path.splitext(filename)[0]
        self.to_format = self.mock_form["to"]

    def run_converter(self, expect_code=None, **kwargs):
        """_summary_
        """

        self.test_converter = FileConverter(files=self.files,
                                            form=self.mock_form,
                                            file_to_convert=FILE_TO_UPLOAD_KEY,
                                            upload_dir=self.tmp_upload_path,
                                            download_dir=self.tmp_download_path,
                                            **kwargs)

        # Check that the input file now exists where we expect it to
        ex_input_filename = os.path.join(self.tmp_upload_path, self.filename)
        assert os.path.exists(ex_input_filename)

        if expect_code is None:
            # If we don't expect an error, just try running the converter
            self.test_converter.run()
        else:
            with pytest.raises(FileConverterAbortException) as esc_info:
                self.test_converter.run()
            assert esc_info.value.status_code == expect_code

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

        ex_input_filename = os.path.join(self.tmp_upload_path, self.files[FILE_TO_UPLOAD_KEY].filename)

        ex_output_filename_base = os.path.splitext(self.files[FILE_TO_UPLOAD_KEY].filename)[0]
        ex_output_filename = os.path.join(self.tmp_download_path, f"{ex_output_filename_base}.{self.to_format}")

        for check_condition, filename in ((input_exist, ex_input_filename),
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

        self.global_log_filename = GLOBAL_LOG_FILENAME
        self.local_log_filename = os.path.join(self.tmp_download_path,
                                               f"{self.filename}-{self.filename_base}.{self.to_format}.{LOCAL_LOG_EXT}")
        self.output_log_filename = os.path.join(self.tmp_download_path,
                                                f"{self.filename_base}.{OUTPUT_LOG_EXT}")

        self.global_log_text: str | None = None
        self.local_log_text: str | None = None
        self.output_log_text: str | None = None

        for log_type in ("global", "local", "output"):
            setattr(self, f"{log_type}_log_text", open(getattr(self, f"{log_type}_log_filename")).read())

    def test_mmcif_to_pdb(self):
        """Run a test of the converter on a straightforward `.mmcif` to `.pdb` conversion
        """

        self.get_input_info(filename="1NE6.mmcif")

        self.run_converter()

        # Check that the input file has been deleted and the output file exists where we expect it to
        self.check_file_status(input_exist=False, output_exist=True)

        # Check that the logs are as we expect
        self.get_logs()

        # Check that the global log file is empty
        assert len(self.global_log_text) == 0

        # Check that the local log and output logs contain expected information

        for log_text in (self.local_log_text, self.output_log_text):
            assert re.compile(r"File name:\s+"+self.filename_base).search(log_text)
            assert "Open Babel Warning" in log_text
            assert "Failed to kekulize aromatic bonds" in log_text

            # Check that we only have the timestamp in the local log, not the output log
            timestamp_re = re.compile(DATETIME_RE_RAW)
            if log_text is self.local_log_text:
                assert timestamp_re.search(log_text)
            else:
                assert not timestamp_re.search(log_text)

    def test_exceed_output_file_size(self):
        """Run a test of the converter to ensure it reports an error properly if the output file size is too large
        """

        self.get_input_info(filename="1NE6.mmcif")

        self.run_converter(expect_code=STATUS_CODE_SIZE,
                           max_file_size=0)

        # Check that the input and output files have properly been deleted
        self.check_file_status(input_exist=False, output_exist=False)

        # Check that the logs are as we expect
        self.get_logs()

        # Check that all logs contain the expected error
        for log_type in ("global", "local", "output"):
            log_text = getattr(self, f"{log_type}_log_text")
            assert "Output file exceeds maximum size" in log_text, ("Did not find expected error message in "
                                                                    f"{log_type} log at " +
                                                                    getattr(self, f"{log_type}_log_filename"))

    def test_invalid_converter(self):
        """Run a test of the converter to ensure it reports an error properly if an invalid converter is requested
        """

        self.get_input_info(filename="1NE6.mmcif",
                            converter="INVALID")

        self.run_converter(expect_code=STATUS_CODE_BAD_METHOD)

        # Check that the input and output files have properly been deleted
        self.check_file_status(input_exist=False, output_exist=False)

        # Check that the logs are as we expect
        self.get_logs()

        # Check that all logs contain the expected error
        for log_type in ("global", "local", "output"):
            log_text = getattr(self, f"{log_type}_log_text")
            assert "ERROR: Unknown converter" in log_text, ("Did not find expected error message in "
                                                            f"{log_type} log at " +
                                                            getattr(self, f"{log_type}_log_filename"))

    def test_xyz_to_inchi(self):
        """Run a test of the converter on a straightforward `.xyz` to `.inchi` conversion
        """

        self.get_input_info(filename="quartz.xyz",
                            to="inchi")

        # "from" is a reserved word so we can't set it as a kwarg in the function call above
        self.mock_form["from"] = "xyz"

        self.run_converter()

        # Check that the input file has been deleted and the output file exists where we expect it to
        self.check_file_status(input_exist=False, output_exist=True)

    def test_atomsk(self):
        """Run a test of the Atomsk converter on a straightforward `.pdb` to `.cif` conversion
        """

        self.get_input_info(filename="hemoglobin.pdb",
                            to="cif",
                            converter=CONVERTER_ATO)

        # "from" is a reserved word so we can't set it as a kwarg in the function call above
        self.mock_form["from"] = "pdb"

        self.run_converter()

        # Check that the input file has been deleted and the output file exists where we expect it to
        self.check_file_status(input_exist=False, output_exist=True)

    def test_xyz_to_inchi_err(self):
        """Run a test of the converter on an `.xyz` to `.inchi` conversion we expect to fail
        """

        self.get_input_info(filename="quartz_err.xyz",
                            to="inchi")

        # "from" is a reserved word so we can't set it as a kwarg in the function call above
        self.mock_form["from"] = "xyz"

        self.run_converter(expect_code=STATUS_CODE_GENERAL)

        # Check that the input and output files have properly been deleted
        self.check_file_status(input_exist=False, output_exist=False)
