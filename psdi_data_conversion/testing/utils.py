"""
# utils.py

This module defines general classes and methods used for unit tests.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from collections.abc import Callable, Iterable
import os
from tempfile import TemporaryDirectory
from typing import Any

import pytest

from psdi_data_conversion.constants import CONVERTER_DEFAULT, GLOBAL_LOG_FILENAME, OUTPUT_LOG_EXT
from psdi_data_conversion.converter import run_converter
from psdi_data_conversion.testing.constants import INPUT_TEST_DATA_LOC


@dataclass
class ConversionTestInfo:
    """Information about a tested conversion."""

    test_spec: SingleConversionTestSpec
    """The specification of the test conversion which was run to produce this"""

    input_dir: str
    """The directory used to store input data for the test"""

    output_dir: str
    """The directory used to create output data in for the test"""

    success: bool = True
    """Whether or not the conversion was successful"""

    @property
    def qualified_in_filename(self):
        """Get the fully-qualified name of the input file"""
        return os.path.join(self.input_dir, self.test_spec.filename)

    @property
    def qualified_out_filename(self):
        """Get the fully-qualified name of the output file"""
        return os.path.join(self.output_dir, self.test_spec.out_filename)

    @property
    def qualified_log_filename(self):
        """Get the fully-qualified name of the log file"""
        return os.path.join(self.output_dir, self.test_spec.log_filename)

    @property
    def qualified_global_log_filename(self):
        """Get the fully-qualified name of the log file"""
        return self.test_spec.global_log_filename


@dataclass
class LibraryConversionTestInfo(ConversionTestInfo):
    """Information about a tested conversion, specifically for when it was tested through a call to the library"""

    exc_info: pytest.ExceptionInfo | None = None
    """If the test conversion raised an exception, that exception's info, otherwise None"""


@dataclass
class CLAConversionTestInfo(ConversionTestInfo):
    """Information about a tested conversion, specifically for when it was tested through a the command-line
    application (CLA)
    """

    captured_stdout: str | None = None
    """If the test was run through the CLA, any output to stdout, otherwise None"""

    captured_stderr: str | None = None
    """If the test was run through the CLA, any output to stderr, otherwise None"""


@dataclass
class GUIConversionTestInfo(ConversionTestInfo):
    """Information about a tested conversion, specifically for when it was tested through the GUI (the local version of
    the web app)
    """


@dataclass
class ConversionTestSpec:
    """Class providing a specification for a test file conversion.

    All attributes of this class can be provided either as a single value or a list of values. In the case that a list
    is provided for one or more attributes, the lists must all be the same length, and they will be iterated through
    (as if using zip on the multiple lists) to test each element in turn.
    """

    filename: str | Iterable[str]
    """The name of the input file, relative to the input test data location, or a list thereof"""

    to_format: str | Iterable[str]
    """The format to test converting the input file to, or a list thereof"""

    name: str | Iterable[str] = CONVERTER_DEFAULT
    """The name of the converter to be used for the test, or a list thereof"""

    conversion_kwargs: dict[str, Any] | Iterable[dict[str, Any]] = field(default_factory=dict)
    """Any keyword arguments to be provided to the call to `run_converter`, aside from those listed above, or a list
    thereof"""

    expect_success: bool | Iterable[bool] = True
    """Whether or not to expect the test to succeed"""

    post_conversion_callback: (Callable[[ConversionTestInfo], str] |
                               Iterable[Callable[[ConversionTestInfo], str]] | None) = None
    """Function to be called after the conversion is performed to check in detail whether results are as expected. It
    should take as its only argument a `ConversionTestInfo` and return a string. The string should be empty if the check
    is passed and should explain the failure otherwise."""

    def __post_init__(self):
        """Regularize the lengths of all attribute lists, in case some were provided as single values and others as
        lists, and set up initial values
        """

        # To ease maintainability, we get the list of this class's attributes automatically from its __dict__, excluding
        # any which start with an underscore
        self._l_attr_names: list[str] = [attr_name for attr_name in self.__dict__ if not attr_name.startswith("_")]

        l_single_val_attrs = []
        self._len: int = 1

        # Check if each attribute of this class is provided as a list, and if any are, make sure that all lists are
        # the same length
        for attr_name in self._l_attr_names:
            val = getattr(self, attr_name)

            val_len = 1

            # Check first if the attr is a str or a dict, which are iterable, but are single-values for the purpose
            # of values here
            if isinstance(val, (str, dict)):
                l_single_val_attrs.append(attr_name)
            else:
                # It's not a str or a dict, so test if we can get the length of it, which indicates it is iterable
                try:
                    val_len = len(val)
                    # If it's a single value in a list, unpack it for now
                    if val_len == 1:
                        # Pylint for some reason things `Any` objects aren't subscriptable, but here we know it is
                        val: Iterable[Any]
                        setattr(self, attr_name, val[0])
                except TypeError:
                    l_single_val_attrs.append(attr_name)

            # Check if there are any conflicts with some lists being provided as different lengths
            if (self._len > 1) and (val_len > 1) and (val_len != self._len):
                raise ValueError("All lists of values which are set as attributes for a `ConversionTestSpec` must be "
                                 "the same length.")
            self._len = val_len

        # At this point, self._len will be either 1 if all attrs are single values, or the length of the lists for attrs
        # that aren't. To keep everything regularised, we make everything a list of this length
        for attr_name in self._l_attr_names:
            if attr_name in l_single_val_attrs:
                setattr(self, attr_name, [getattr(self, attr_name)]*self._len)

    def __len__(self):
        """Get the length from the member - valid only after `__post_init__` has been called"""
        return self._len

    def __iter__(self):
        """Allow to iterate over the class, getting a `SingleConversionTestSpec` for each value
        """
        l_l_attr_vals = zip(*[getattr(self, attr_name) for attr_name in self._l_attr_names])
        for l_attr_vals in l_l_attr_vals:
            yield SingleConversionTestSpec(**dict(zip(self._l_attr_names, l_attr_vals)))


@dataclass
class SingleConversionTestSpec:
    """Class providing a specification for a single test file conversion, produced by iterating over a
    `ConversionTestSpec` object
    """

    filename: str
    """The name of the input file, relative to the input test data location"""

    to_format: str
    """The format to test converting the input file to"""

    name: str | Iterable[str] = CONVERTER_DEFAULT
    """The name of the converter to be used for the test"""

    conversion_kwargs: dict[str, Any] = field(default_factory=dict)
    """Any keyword arguments to be provided to the call to `run_converter`, aside from those listed above and
    `upload_dir` and `download_dir` (for which temporary directories are used)"""

    expect_success: bool = True
    """Whether or not to expect the test to succeed"""

    post_conversion_callback: (Callable[[ConversionTestInfo], str] | None) = None
    """Function to be called after the conversion is performed to check in detail whether results are as expected. It
    should take as its only argument a `ConversionTestInfo` and return a string. The string should be empty if the check
    is passed and should explain the failure otherwise."""

    @property
    def out_filename(self) -> str:
        """The unqualified name of the output file which should have been created by the conversion."""
        return f"{os.path.splitext(self.filename)[0]}.{self.to_format}"

    @property
    def log_filename(self) -> str:
        """The unqualified name of the log file which should have been created by the conversion."""
        return f"{os.path.splitext(self.filename)[0]}{OUTPUT_LOG_EXT}"

    @property
    def global_log_filename(self) -> str:
        """The unqualified name of the global log file which stores info on all conversions."""
        return GLOBAL_LOG_FILENAME


def run_test_conversion_with_library(test_spec: ConversionTestSpec):
    """Runs a test conversion or series thereof through a call to the python library's `run_converter` function.

    Parameters
    ----------
    test_spec : ConversionTestSpec
        The specification for the test or series of tests to be run
    """
    # Make temporary directories for the input and output files to be stored in
    with TemporaryDirectory("_input") as input_dir, TemporaryDirectory("_output") as output_dir:
        # Iterate over the test spec to run each individual test it defines
        for single_test_spec in test_spec:
            _run_single_test_conversion_with_library(test_spec=single_test_spec,
                                                     input_dir=input_dir,
                                                     output_dir=output_dir)


def _run_single_test_conversion_with_library(test_spec: SingleConversionTestSpec,
                                             input_dir: str,
                                             output_dir: str):
    """Runs a single test conversion through a call to the python library's `run_converter` function.

    Parameters
    ----------
    test_spec : _SingleConversionTestSpec
        The specification for the test to be run
    input_dir : str
        A directory which can be used to store input data
    output_dir : str
        A directory which can be used to create output data
    """

    # Symlink the input file to the input directory
    qualified_in_filename = os.path.join(input_dir, test_spec.filename)
    os.symlink(os.path.join(INPUT_TEST_DATA_LOC, test_spec.filename),
               qualified_in_filename)

    exc_info: pytest.ExceptionInfo | None = None
    if test_spec.expect_success:
        run_converter(filename=qualified_in_filename,
                      to_format=test_spec.to_format,
                      name=test_spec.name,
                      upload_dir=input_dir,
                      download_dir=output_dir,
                      **test_spec.conversion_kwargs)
        success = True
    else:
        with pytest.raises(Exception) as exc_info:
            run_converter(filename=qualified_in_filename,
                          to_format=test_spec.to_format,
                          name=test_spec.name,
                          upload_dir=input_dir,
                          download_dir=output_dir,
                          **test_spec.conversion_kwargs)
        success = False

    # Compile output info for the test and call the callback function if one is provided
    if test_spec.post_conversion_callback:
        test_info = LibraryConversionTestInfo(test_spec=test_spec,
                                              input_dir=input_dir,
                                              output_dir=output_dir,
                                              success=success,
                                              exc_info=exc_info)
        callback_msg = test_spec.post_conversion_callback(test_info)
        assert not callback_msg, callback_msg
