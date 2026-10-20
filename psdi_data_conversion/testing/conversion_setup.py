"""psdi_data_conversion/testing/conversion_setup.py
===============

Created 2026-10-20 by Bryan Gillis.

This module provides functions and callable classes which can be used as setup callbacks for test specs.
"""

import os
from collections.abc import Callable, Iterable
from typing import Concatenate

from psdi_data_conversion.database import get_format_info
from psdi_data_conversion.testing.utils import SingleConversionTestSpec
from psdi_data_conversion.utils import replace_extension

JUNK_FILE_CONTENT = "TEST_JUNK_CONTENT"
"""A text string added to junk files created, which can be checked for to see if output files contain this string"""


def write_junk_file(dir: str, filename: str):
    """Write out a junk file at a desired location"""
    open(os.path.join(dir, filename), "w").write(JUNK_FILE_CONTENT)


class MultiSetup:
    """Callable class which stores a list of other setup tasks to call, and runs them all in sequence. All tasks will
    be run, with no try blocks being applied to stop exceptions from propagating upwards as soon as they happen
    """

    def __init__(self, *l_setups: Callable[Concatenate[SingleConversionTestSpec, ...], None]):

        self.l_setups: Iterable[Callable[Concatenate[SingleConversionTestSpec, ...], None]] = l_setups
        """A list of setup tasks to be called in turn"""

    def __call__(self, test_spec: SingleConversionTestSpec, **kwargs) -> None:
        """When called and passed a test spec, run each stored task on the test spec in turn"""

        for setup in self.l_setups:
            setup(test_spec, **kwargs)


class BlockOutput:
    """Callable class which will block the intended output filename of a conversion with a junk file before the
    conversion is performed
    """

    def __call__(self, test_spec: SingleConversionTestSpec, output_dir, **_) -> None:
        # Get the expected name of the output file
        if test_spec.ex_out_filename:
            ex_out_filename = test_spec.ex_out_filename
        else:
            ex_out_filename = replace_extension(test_spec.filename, get_format_info(test_spec.to_format).name)

        write_junk_file(output_dir, ex_out_filename)


class BlockOutputExt:
    """Callable class which will block the files in the output directory which have the output filename with its
    extension replaced by any in the provided list
    """

    def __init__(self, l_ext: str | Iterable[str]):
        if isinstance(l_ext, str):
            self.l_ext = [l_ext]
        else:
            self.l_ext = l_ext

    def __call__(self, test_spec: SingleConversionTestSpec, output_dir, **_) -> None:
        for ext in self.l_ext:
            write_junk_file(output_dir, replace_extension(test_spec.filename, ext))
