"""
# utils.py

This module defines general classes and methods used for unit tests.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from collections.abc import Callable, Iterable
import os
from typing import Any

from psdi_data_conversion.constants import CONVERTER_DEFAULT, LOG_EXT


@dataclass
class ConversionTestInfo:
    """Information about a tested conversion."""

    test_spec: ConversionTestSpec
    """The specification of the test conversion which was run to produce this"""

    to_format: str
    """The format the input file was tested to be converted to"""

    name: str
    """The name of the converter to used for the test"""

    success: bool = True
    """Whether or not the conversion was successful"""

    conversion_kwargs: dict[str | Any] = field(default_factory=dict)
    """Any keyword arguments provided to the call to `run_converter`, aside from those listed above"""

    def __post_init__(self):
        """Set some properties from the test spec.
        """
        self.filename = self.test_spec.filename

    @property
    def out_filename(self) -> str:
        """The unqualified name of the output file which should have been created by the conversion."""
        return f"{os.path.splitext(self.filename)[0]}.{self.to_format}"

    @property
    def log_filename(self) -> str:
        """The unqualified name of the log file which should have been created by the conversion."""
        return f"{os.path.splitext(self.filename)[0]}.{LOG_EXT}"


@dataclass
class LibraryConversionTestInfo(ConversionTestInfo):
    """Information about a tested conversion, specifically for when it was tested through a call to the library"""

    caught_exception: Exception | None = None
    """If the test conversion raised an exception, that exception, otherwise None"""


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
                setattr(self, attr_name, [attr_name]*self._len)

    def __len__(self):
        """Get the length from the member - valid only after `__post_init__` has been called"""
        return self._len

    def __iter__(self):
        """Allow to iterate over the class, getting a `SingleConversionTestSpec` for each value
        """
        l_l_attr_vals = zip(*[getattr(self, attr_name) for attr_name in self._l_attr_names])
        while True:
            l_attr_vals = next(l_l_attr_vals)
            yield _SingleConversionTestSpec(**dict(zip(self._l_attr_names, l_attr_vals)))


@dataclass
class _SingleConversionTestSpec:
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
    """Any keyword arguments to be provided to the call to `run_converter`, aside from those listed above"""

    expect_success: bool = True
    """Whether or not to expect the test to succeed"""

    post_conversion_callback: (Callable[[ConversionTestInfo], str] | None) = None
    """Function to be called after the conversion is performed to check in detail whether results are as expected. It
    should take as its only argument a `ConversionTestInfo` and return a string. The string should be empty if the check
    is passed and should explain the failure otherwise."""


def run_test_conversion_with_library(test_spec: ConversionTestSpec):
    """Runs a test conversion or series thereof through a call to the python library's `run_converter` function.

    Parameters
    ----------
    test_spec : ConversionTestSpec
        The specification for the test or series of tests to be run
    """

    # Iterate over the test spec to run each individual test it defines
    for single_test_spec in test_spec:
        _run_single_test_conversion_with_library(single_test_spec)


def _run_single_test_conversion_with_library(single_test_spec: _SingleConversionTestSpec):
    """Runs a single test conversion through a call to the python library's `run_converter` function.

    Parameters
    ----------
    single_test_spec : _SingleConversionTestSpec
        The specification for the test to be run
    """
    pass
