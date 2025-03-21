"""
# conversion_callbacks.py

This module provides functions and callable classes which can be used as callbacks for to check the results of a
conversion test, run with the functions and classes defined in the `utils.py` module.
"""

from collections.abc import Callable, Iterable
from dataclasses import dataclass
from genericpath import isfile
import os

from psdi_data_conversion.testing.utils import ConversionTestInfo


@dataclass
class MultiCallback:
    """Callable class which stores a list of other callbacks to call on the input, and runs them all in sequence. All
    callbacks will be run, even if some fail, and the results will be joined to an output string.
    """

    l_callbacks: Iterable[Callable[[ConversionTestInfo]]]
    """A list of callbacks to be called in turn"""

    def __call__(self, test_info: ConversionTestInfo) -> str:
        """When called and passed test info, run each stored callback on the test info in turn and join the results
        """

        l_results = [callback(test_info) for callback in self.l_callbacks]

        # Only join non-empty results so we don't end up with excess newlines
        res = "\n".join([x for x in l_results if x]).strip()

        return res


@dataclass
class CheckOutputStatus:
    """Callable class which checks the status of standard output files from a conversion"""

    expect_output_exists: bool = True
    """Whether to expect that the output file of the conversion exists (with non-zero size) or not, by default True"""

    expect_log_exists: bool = True
    """Whether to expect that the log exists (with non-zero size) or not, by default True"""

    def __call__(self, test_info: ConversionTestInfo) -> str:
        """Perform the check on output file and log status"""

        l_errors: list[str] = []

        # Check the status of the output file
        qualified_out_filename = test_info.qualified_out_filename
        if self.expect_output_exists:
            if not os.path.isfile(qualified_out_filename):
                l_errors.append(f"ERROR: Expected output file from conversion '{qualified_out_filename}' does not "
                                "exist")
            elif os.path.getsize(qualified_out_filename) == 0:
                l_errors.append(f"ERROR: Expected output file from conversion '{qualified_out_filename}' exists but "
                                "is unexpectedly empty")
        elif os.path.isfile(qualified_out_filename):
            l_errors.append(f"ERROR: Output file from conversion '{qualified_out_filename}' exists, but was expected "
                            "to not exist")

        qualified_log_filename = test_info.qualified_log_filename
        if self.expect_log_exists:
            if not os.path.isfile(qualified_log_filename):
                l_errors.append(f"ERROR: Expected log file from conversion '{qualified_log_filename}' does not "
                                "exist")
            elif os.path.getsize(qualified_log_filename) == 0:
                l_errors.append(f"ERROR: Expected log file from conversion '{qualified_log_filename}' exists but "
                                "is unexpectedly empty")
        elif os.path.isfile(qualified_log_filename):
            l_errors.append(f"ERROR: Log file from conversion '{qualified_log_filename}' exists, but was expected "
                            "to not exist")

        # Join any errors for output
        res = "\n".join(l_errors)
        return res
