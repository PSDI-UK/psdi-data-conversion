"""@file psdi-data-conversion/tests/converter_test.py

Created 2024-12-15 by Bryan Gillis.

Unit tests of the converter class
"""

from copy import deepcopy
import logging
import math
import os
import re
import shutil
import pytest

from psdi_data_conversion import constants as const
from psdi_data_conversion.converter import L_REGISTERED_CONVERTERS, run_converter
from psdi_data_conversion.converters.atomsk import CONVERTER_ATO
from psdi_data_conversion.converters.base import FileConverterAbortException
from psdi_data_conversion.converters.c2x import CONVERTER_C2X
from psdi_data_conversion.converters.openbabel import CONVERTER_OB, OBFileConverter
from psdi_data_conversion.file_io import split_archive_ext
from psdi_data_conversion.main import FileConverterInputException


@pytest.fixture()
def base_mock_data():
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
def base_input():
    """A fixture providing a dict of default input values for the converter which aren't provided in the `data`
    """
    return {'from_format': None,
            'to_format': 'pdb',
            'out_format': 'pdb'}


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
    def setup_test(self, base_mock_data, base_input, tmp_upload_path, tmp_download_path, test_data_loc):
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

        # Save test data location
        self.test_data_loc = test_data_loc

        # Save fixtures
        self.base_mock_data = base_mock_data
        self.base_input = base_input
        self.tmp_upload_path = tmp_upload_path
        self.tmp_download_path = tmp_download_path

    def get_input_info(self, filename: str, **kwargs,):
        """Sets up a mock `data` for input and gets various variables we'll want to use for checks on output

        Parameters
        ----------
        filename : str
            The name of the file to use as input for the test
        """

        self.mock_data = deepcopy(self.base_mock_data)
        self.input = deepcopy(self.base_input)

        for key, value in kwargs.items():
            if key in self.mock_data:
                self.mock_data[key] = value
            elif key in self.input:
                self.input[key] = value
            else:
                raise RuntimeError(f"Invalid key {key} provided for data or input")

        # Save some variables from input we'll be using throughout this test
        self.local_filename = filename
        self.source_filename = os.path.join(self.test_data_loc, filename)
        self.filename_base = split_archive_ext(filename)[0]
        self.to_format = self.input["to_format"]
        self.out_format = self.input["out_format"]

    def get_converter_kwargs(self, **kwargs):
        """Get the keyword arguments to be passed to a FileConverter for testing
        """
        standard_kwargs = {"filename": self.source_filename,
                           "data": self.mock_data,
                           "upload_dir": self.tmp_upload_path,
                           "download_dir": self.tmp_download_path}
        standard_kwargs.update(self.input)
        standard_kwargs.update(kwargs)
        return standard_kwargs

    def run_converter(self, expect_exception=None, expect_code=None, **kwargs):
        """Runs a test on a file converter and checks that it returns successfully or else fails with an expected error
        code.
        """

        converter_kwargs = self.get_converter_kwargs(**kwargs)

        if expect_exception is None and expect_code is None:
            # If we don't expect an error, just try running the converter
            run_converter(**converter_kwargs)
        elif expect_exception is not None:
            with pytest.raises(expect_exception) as esc_info:
                run_converter(**converter_kwargs)
        else:
            with pytest.raises(FileConverterAbortException) as esc_info:
                run_converter(**converter_kwargs)
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

    def test_mmcif_to_pdb(self):
        """Run a test of the converter on a straightforward `.mmcif` to `.pdb` conversion
        """

        self.get_input_info(filename="1NE6.mmcif")

        self.run_converter()

        # Check that the input file has not been deleted and the output file exists where we expect it to
        self.check_file_status(input_exist=True, output_exist=True)

        # Check that the logs are as we expect
        self.get_logs()

        # Check that the global log file is empty
        assert len(self.global_log_text) == 0

        # Check that the output log contains expected information
        assert re.compile(r"File name:\s+"+self.filename_base).search(self.output_log_text)
        assert "Open Babel Warning" in self.output_log_text
        assert "Failed to kekulize aromatic bonds" in self.output_log_text

        timestamp_re = re.compile(const.DATETIME_RE_RAW)
        assert timestamp_re.search(self.output_log_text)

    def test_exceed_output_file_size(self):
        """Run a test of the converter to ensure it reports an error properly if the output file size is too large
        """

        self.get_input_info(filename="1NE6.mmcif")

        self.run_converter(expect_code=const.STATUS_CODE_SIZE,
                           max_file_size=0.0001)

        # Check that the input file remains but no output file has been created
        self.check_file_status(input_exist=True, output_exist=False)

        # Check that the logs are as we expect
        self.get_logs()

        # Check that all logs contain the expected error
        for log_type in ("global", "output"):
            log_text = getattr(self, f"{log_type}_log_text")
            assert "Output file exceeds maximum size" in log_text, ("Did not find expected error message in "
                                                                    f"{log_type} log at " +
                                                                    getattr(self, f"{log_type}_log_filename"))

        # Now check that 0 properly works to indicate unlimited size
        self.get_input_info(filename="1NE6.mmcif")
        self.run_converter(max_file_size=0)
        self.check_file_status(input_exist=True, output_exist=True)

        # Check that it stops partway through when running on an archive. The tarball here is smaller than the size
        # limit, but the unpacked files are larger
        self.get_input_info(filename="caffeine-smi.tar.gz")
        self.run_converter(expect_code=const.STATUS_CODE_SIZE,
                           max_file_size=0.0005)
        self.check_file_status(input_exist=True, output_exist=False)
        self.get_logs()
        assert "Output file exceeds maximum size" in self.output_log_text

    def test_invalid_converter(self):
        """Run a test of the converter to ensure it reports an error properly if an invalid converter is requested
        """

        self.get_input_info(filename="1NE6.mmcif")

        self.run_converter(expect_exception=FileConverterInputException,
                           name="INVALID")

        # Check that the input file remains but no output file has been created
        self.check_file_status(input_exist=True, output_exist=False)

    def test_xyz_to_inchi(self):
        """Run a test of the converter on an `.xyz` to `.inchi` conversion which we expect to have warnings about data
        loss and extrapolation
        """

        self.get_input_info(filename="quartz.xyz",
                            to_format="inchi",
                            out_format="inchi")

        self.run_converter()

        # Check that the input file has not been deleted and the output file exists where we expect it to
        self.check_file_status(input_exist=True, output_exist=True)

        self.get_logs()

        # XYZ has the 2D and 3D properties while inchi doesn't, while inchi has the Connections property while xyz
        # doesn't, so check these exist in the logs as a warning
        assert "WARNING" in self.output_log_text
        assert const.QUAL_NOTE_OUT_MISSING.format(const.QUAL_2D_LABEL) in self.output_log_text
        assert const.QUAL_NOTE_OUT_MISSING.format(const.QUAL_3D_LABEL) in self.output_log_text
        assert const.QUAL_NOTE_IN_MISSING.format(const.QUAL_CONN_LABEL) in self.output_log_text

    def test_cleanup(self):
        """Test that input files are deleted if requested
        """

        self.get_input_info(filename="hemoglobin.pdb",
                            to_format="cif",
                            out_format="cif")

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

    def test_atomsk(self):
        """Run tests of Atomsk on some conversions we expect to succeed without issue
        """

        # pdb to cif

        self.get_input_info(filename="hemoglobin.pdb",
                            to_format="cif",
                            out_format="cif")

        self.run_converter(name=CONVERTER_ATO)

        # Check that the input file has not been deleted and the output file exists where we expect it to
        self.check_file_status(input_exist=True, output_exist=True)

        # cif to xyz

        self.get_input_info(filename="nacl.cif",
                            to_format="xyz",
                            out_format="xyz")
        self.run_converter(name=CONVERTER_ATO)
        self.check_file_status(input_exist=True, output_exist=True)

    def test_c2x(self):
        """Run tests of c2x on some conversions we expect to succeed without issue
        """

        self.get_input_info(filename="hemoglobin.pdb",
                            to_format="cif",
                            out_format="cif")

        self.run_converter(name=CONVERTER_C2X)

        # Check that the input file has not been deleted and the output file exists where we expect it to
        self.check_file_status(input_exist=True, output_exist=True)

        # cif to cml

        self.get_input_info(filename="nacl.cif",
                            to_format="cml",
                            out_format="cml")
        self.run_converter(name=CONVERTER_C2X)
        self.check_file_status(input_exist=True, output_exist=True)

    def test_xyz_to_inchi_err(self):
        """Run a test of the converter on an `.xyz` to `.inchi` conversion we expect to fail
        """

        self.get_input_info(filename="quartz_err.xyz",
                            to_format="inchi",
                            out_format="inchi")

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
