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

    caught_exception: Exception | None = None
    """If the test conversion raised an exception, that exception, otherwise None"""

    captured_stdout: str | None = None
    """If the test was run through the CLA, any output to stdout, otherwise None"""

    captured_stderr: str | None = None
    """If the test was run through the CLA, any output to stderr, otherwise None"""

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
class ConversionTestSpec:
    """Class providing a specification for a test file conversion.
    """

    filename: str | Iterable[str]
    """The name of the input file, relative to the input test data location, or a list thereof"""

    to_format: str | Iterable[str]
    """The format to test converting the input file to, or a list thereof"""

    name: str | Iterable[str] = CONVERTER_DEFAULT
    """The name of the converter to be used for the test, or a list thereof"""

    conversion_kwargs: dict[str, Any] | list[dict[str, Any]] = field(default_factory=dict)
    """Any keyword arguments to be provided to the call to `run_converter`, aside from those listed above, or a list
    thereof"""

    test_all_combinations: bool = False
    """When more than one of `to_format`, `name`, and `conversion_kwargs` is provided as a list, the default behaviour
    (when this is set to False) is to zip through the lists, e.g.:
    
    ```
    to_format = ["pdb", "mmcif"]
    name = ["Open Babel", "c2x"]
    ```
    
    will run one test of converting to "pdb" with the "Open Babel" converter, and one test of converting to "mmcif" with
    the "c2x" converter.

    When this is set to True, instead all possibilities between the lists will be tested, so the above will run four
    test conversions:
    - Converting to "pdb" with "Open Babel"
    - Converting to "pdb" with "c2x"
    - Converting to "mmcif" with "Open Babel"
    - Converting to "mmcif" with "c2x"
    """

    expect_success: bool = True
    """Whether or not to expect the test to succeed"""

    post_conversion_callback: Callable[[ConversionTestInfo], str] | None = None
    """Function to be called after the conversion is performed to check in detail whether results are as expected. It
    should take as its only argument a `ConversionTestInfo` and return a string. The string should be empty if the check
    is passed and should explain the failure otherwise."""
