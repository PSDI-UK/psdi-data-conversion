"""
# conversion_test_spec.py

This module defines a class to provide a specification for a test conversion to be run.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from collections.abc import Callable, Iterable
from typing import Any

from psdi_data_conversion.constants import CONVERTER_DEFAULT

@dataclass
class ConversionTestSpec:
    """Class providing a specification for a test file conversion.
    """

    filename: str
    """The name of the input file, relative to the input test data location"""

    to_format: str | Iterable[str]
    """The format to test converting the input file to, or a list thereof"""

    name: str | Iterable[str] = CONVERTER_DEFAULT
    """The name of the converter to be used for the test"""

    conversion_kwargs: dict[str|Any] | list[dict[str|Any]] = field(default_factory=dict)
    """Any keyword arguments to be provided to the call to `run_converter`, aside from those listed above."""

    expect_success: bool = True
    """Whether or not to expect the test to succeed."""

    post_conversion_callback: Callable[[ConversionTestSpec], str] | None = None
    """Function to be called after the conversion is performed to check in detail whether results are as expected. It
    should take as its first argument a `ConversionTestSpec` (this object) and return a string. The string should be
    empty if the check is passed and should explain the failure otherwise."""