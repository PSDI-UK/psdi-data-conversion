"""
# conversion_callbacks.py

This module provides functions and callable classes which can be used as callbacks for to check the results of a
conversion test, run with the functions and classes defined in the `utils.py` module.
"""

import abc
import os
import re
from collections.abc import Callable, Iterable
from dataclasses import dataclass, field
from tempfile import TemporaryDirectory

from psdi_data_conversion.constants import DATETIME_RE_RAW
from psdi_data_conversion.file_io import unpack_zip_or_tar
from psdi_data_conversion.log_utility import string_with_placeholders_matches
from psdi_data_conversion.testing.constants import OUTPUT_TEST_DATA_LOC_IN_PROJECT
from psdi_data_conversion.testing.utils import ConversionTestInfo, check_file_match


class MultiCallback:
    """Callable class which stores a list of other callbacks to call on the input, and runs them all in sequence. All
    callbacks will be run, even if some fail, and the results will be joined to an output string.
    """

    def __init__(self, *l_callbacks: Callable[[ConversionTestInfo], str]):

        self.l_callbacks: Iterable[Callable[[ConversionTestInfo], str]] = l_callbacks
        """A list of callbacks to be called in turn"""

    def __call__(self, test_info: ConversionTestInfo) -> str:
        """When called and passed test info, run each stored callback on the test info in turn and join the results
        """

        l_results = [callback(test_info) for callback in self.l_callbacks]

        # Only join non-empty results so we don't end up with excess newlines
        res = "\n".join([x for x in l_results if x]).strip()

        return res


@dataclass
class CheckFileStatus:
    """Callable class which checks the presence or absence of standard input and output files for a conversion"""

    expect_input_exists: bool | None = True
    """Whether to expect that the input file of the conversion exists (with non-zero size) or not. If None, will
    not check either way."""

    expect_output_exists: bool | None = True
    """Whether to expect that the output file of the conversion exists (with non-zero size) or not. If None, will
    not check either way."""

    expect_log_exists: bool | None = True
    """Whether to expect that the log exists (with non-zero size) or not. If None, will not check either way."""

    expect_global_log_exists: bool | None = None
    """Whether to expect that the global log exists (zero-size allowed) or not. If None, will not check either way."""

    def __call__(self, test_info: ConversionTestInfo) -> str:
        """Perform the check on output file and log status"""

        l_errors: list[str] = []

        # Check the status of the input file
        qualified_in_filename = test_info.qualified_in_filename
        if self.expect_input_exists:
            if not os.path.isfile(qualified_in_filename):
                l_errors.append(f"ERROR: Expected input file for conversion '{qualified_in_filename}' does not "
                                "exist")
            elif os.path.getsize(qualified_in_filename) == 0:
                l_errors.append(f"ERROR: Expected input file for conversion '{qualified_in_filename}' exists but "
                                "is unexpectedly empty")
        elif self.expect_output_exists is False and os.path.isfile(qualified_in_filename):
            l_errors.append(f"ERROR: Input file from conversion '{qualified_in_filename}' exists, but was expected "
                            "to not exist")

        # Check the status of the output file
        qualified_out_filename = test_info.qualified_out_filename
        if self.expect_output_exists:
            if not os.path.isfile(qualified_out_filename):
                l_errors.append(f"ERROR: Expected output file from conversion '{qualified_out_filename}' does not "
                                "exist")
            elif os.path.getsize(qualified_out_filename) == 0:
                l_errors.append(f"ERROR: Expected output file from conversion '{qualified_out_filename}' exists but "
                                "is unexpectedly empty")
        elif self.expect_output_exists is False and os.path.isfile(qualified_out_filename):
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
        elif self.expect_log_exists is False and os.path.isfile(qualified_log_filename):
            l_errors.append(f"ERROR: Log file from conversion '{qualified_log_filename}' exists, but was expected "
                            "to not exist")

        qualified_global_log_filename = test_info.qualified_global_log_filename
        if self.expect_global_log_exists:
            if not os.path.isfile(qualified_global_log_filename):
                l_errors.append(f"ERROR: Expected global log file from conversion '{qualified_global_log_filename}' "
                                "does not exist")
        elif self.expect_global_log_exists is False and os.path.isfile(qualified_global_log_filename):
            l_errors.append(f"ERROR: Global log file from conversion '{qualified_global_log_filename}' exists, but was "
                            "expected to not exist")

        # Join any errors for output
        res = "\n".join(l_errors)
        return res


@dataclass
class CheckArchiveContents:
    """Callable class which checks that an archive file created as the result of a conversion contains files with the
    expected names and all of the expected type"""

    l_filename_bases: Iterable[str]
    """List of unqualified filenames without extensions, representing the files that should be found in the archive"""

    to_format: str
    """Format (extension) that all files in the archive should have"""

    def __call__(self, test_info: ConversionTestInfo):
        """Run the check on archive contents"""

        # First, check that the archive file exists
        qualified_out_filename = test_info.qualified_out_filename
        if not os.path.isfile(qualified_out_filename):
            raise FileNotFoundError(f"ERROR: Expected output file from conversion '{qualified_out_filename}' does not "
                                    "exist")

        l_errors = []

        # Use a temporary directory to unpack the archive, then check for each file in it
        with TemporaryDirectory() as extract_dir:
            unpack_zip_or_tar(qualified_out_filename, extract_dir=extract_dir)

            for filename_base in self.l_filename_bases:
                filename = f"{filename_base}.{self.to_format}"
                if not os.path.isfile(os.path.join(extract_dir, filename)):
                    l_errors.append(f"ERROR: Expected file '{filename}' was not found in archive "
                                    f"{qualified_out_filename}")

        return "\n".join(l_errors)


@dataclass
class CheckTextContents(abc.ABC):
    """Callable class which checks the contents of some text (e.g. stdout or a log file) from a conversion"""

    l_strings_to_find: str | Iterable[str] = field(default_factory=list)
    """One or more strings which must be found in the text. These may optionally include formatting placeholders,
    e.g. "The filename is: {file}", in which case any text in the place of the placeholder will be considered valid for
    a match.
    """

    l_strings_to_exclude: str | Iterable[str] = field(default_factory=list)
    """One or more strings which must NOT be found in the text. These may optionally include formatting
    placeholders, e.g. "The filename is: {file}", in which case any text in the place of the placeholder will be
    considered valid for a match.
    """

    l_regex_to_find: str | Iterable[str] = field(default_factory=list)
    """One or more uncompiled regular expressions which must be matched somewhere in the text"""

    l_regex_to_exclude: str | Iterable[str] = field(default_factory=list)
    """One or more uncompiled regular expressions which must NOT be matched anywhere in the text"""

    def __post_init__(self):
        """If any input was provided as just a single string, coerce it to a list"""
        if isinstance(self.l_strings_to_find, str):
            self.l_strings_to_find = [self.l_strings_to_find]
        if isinstance(self.l_strings_to_exclude, str):
            self.l_strings_to_exclude = [self.l_strings_to_exclude]
        if isinstance(self.l_regex_to_find, str):
            self.l_regex_to_find = [self.l_regex_to_find]
        if isinstance(self.l_regex_to_exclude, str):
            self.l_regex_to_exclude = [self.l_regex_to_exclude]

    def get_default_strings_to_find(self, test_info: ConversionTestInfo) -> list[str]:
        """Get a default list of strings to find in the text; can be overridden by child classes to implement
        defaults"""
        return []

    def get_default_strings_to_exclude(self, test_info: ConversionTestInfo) -> list[str]:
        """Get a default list of strings to NOT find in the text; can be overridden by child classes to implement
        defaults"""
        return []

    def get_default_regex_to_find(self, test_info: ConversionTestInfo) -> list[str]:
        """Get a default list of uncompiled regular expressions to match in the text; can be overridden by child
        classes to implement defaults"""
        return []

    def get_default_regex_to_exclude(self, test_info: ConversionTestInfo) -> list[str]:
        """Get a default list of uncompiled regular expressions to NOT match in the text; can be overridden by
        child classes to implement defaults"""
        return []

    @abc.abstractmethod
    def _get_text(self, test_info: ConversionTestInfo) -> str:
        """Abstract method which must be overridden to get the text to be checked"""
        pass

    @abc.abstractmethod
    def _get_text_source_label(self, test_info: ConversionTestInfo) -> str:
        """Abstract method which must be overridden to label the source of the text to be checked"""
        pass

    def __call__(self, test_info: ConversionTestInfo) -> str:
        """Perform the check on text contents"""

        test_text = self._get_text(test_info)

        l_errors: list[str] = []

        # Add the default checks to the lists of strings/regexes to check
        l_strings_to_find: list[str] = list(self.l_strings_to_find) + self.get_default_strings_to_find(test_info)
        l_strings_to_exclude: list[str] = (list(self.l_strings_to_exclude) +
                                           self.get_default_strings_to_exclude(test_info))
        l_regex_to_find: list[str] = list(self.l_regex_to_find) + self.get_default_regex_to_find(test_info)
        l_regex_to_exclude: list[str] = list(self.l_regex_to_exclude) + self.get_default_regex_to_exclude(test_info)

        # Check that all expected strings are present
        for string_to_find in l_strings_to_find:
            if not string_with_placeholders_matches(string_to_find, test_text):
                l_errors.append(f"ERROR: String \"{string_to_find}\" was expected in "
                                f"{self._get_text_source_label(test_info)} but was not found. Text:\n {test_text}")

        # Check that all excluded strings are not present
        for l_strings_to_exclude in l_strings_to_exclude:
            if string_with_placeholders_matches(l_strings_to_exclude, test_text):
                l_errors.append(f"ERROR: String \"{l_strings_to_exclude}\" was not expected in "
                                f"{self._get_text_source_label(test_info)} but was found. Text:\n {test_text}")

        # Check that all expected regexes are present
        for regex_to_find in l_regex_to_find:
            compiled_regex = re.compile(regex_to_find)
            if not compiled_regex.search(test_text):
                l_errors.append(f"ERROR: Regex /{regex_to_find}/ was expected in "
                                f"{self._get_text_source_label(test_info)} but was not found. Text:\n {test_text}")

        # Check that all excluded regexes are not present
        for regex_to_exclude in l_regex_to_exclude:
            compiled_regex = re.compile(regex_to_exclude)
            if compiled_regex.search(test_text):
                l_errors.append(f"ERROR: Regex /{regex_to_exclude}/ was not expected in "
                                f"{self._get_text_source_label(test_info)} but was found. Text:\n {test_text}")

        # Join any errors for output
        res = "\n".join(l_errors)
        return res


@dataclass
class CheckLogContents(CheckTextContents):
    """Implementation of `CheckTextContents` which checks the contents of the output log file of a conversion"""

    def _get_text(self, test_info: ConversionTestInfo) -> str:
        """Get the text from the output log file"""

        # First, check that the log exists
        qualified_log_filename = test_info.qualified_log_filename
        if not os.path.isfile(qualified_log_filename):
            raise FileNotFoundError(f"ERROR: Expected log file from conversion '{qualified_log_filename}' does not "
                                    "exist")

        return open(qualified_log_filename, "r").read()

    def _get_text_source_label(self, test_info: ConversionTestInfo) -> str:
        """Label as coming from the appropriate log file"""
        return f"log file '{test_info.qualified_log_filename}'"


@dataclass
class CheckLogContentsSuccess(CheckLogContents):
    """Specialized callback to check log contents for a successful conversion"""

    def get_default_strings_to_exclude(self, test_info: ConversionTestInfo) -> Iterable[str]:
        """Exclude strings which indicate something likely went wrong"""
        return ["ERROR", "exception", "Exception"]

    def get_default_regex_to_find(self, test_info: ConversionTestInfo) -> Iterable[str]:
        """Check for the filename and date and time in the log"""
        return [r"File name:\s*"+os.path.splitext(test_info.test_spec.filename)[0],
                DATETIME_RE_RAW]


@dataclass
class CheckStdoutContents(CheckTextContents):
    """Implementation of `CheckTextContents` which checks the output to stdout"""

    def _get_text(self, test_info: ConversionTestInfo) -> str:
        """Get the text from the captured stdout"""

        return test_info.captured_stdout

    def _get_text_source_label(self, test_info: ConversionTestInfo) -> str:
        """Label as coming from stdout"""
        return "stdout"


@dataclass
class CheckStderrContents(CheckTextContents):
    """Implementation of `CheckTextContents` which checks the output to stderr"""

    def _get_text(self, test_info: ConversionTestInfo) -> str:
        """Get the text from the captured stderr"""

        return test_info.captured_stderr

    def _get_text_source_label(self, test_info: ConversionTestInfo) -> str:
        """Label as coming from stderr"""
        return "stderr"


@dataclass
class CheckException:
    """Callable class which checks an exception raised for its type, status code, and message. Tests will only be
    run on tests with the python library, as that's the only route that provides exceptions."""

    ex_type: type[Exception] | None = None
    """The expected type of the raised exception (subclasses of it will also be allowed)"""

    ex_message: str | None = None
    """A string pattern (with placeholders allowed) which must be found in the exceptions message"""

    ex_status_code: int | None = None
    """The required status code of the exception - will not pass if this is provided and the exception doesn't come
    with a status code"""

    def __call__(self, test_info: ConversionTestInfo) -> str:
        """Perform the check on the exception"""

        # Skip check on CLA, since this won't catch any exceptions
        if test_info.run_type == "cla":
            return ""

        # Confirm that an exception was indeed raised
        exc_info = test_info.exc_info
        if not exc_info:
            return f"ERROR: No exception was raised - expected an exception of type {self.ex_type}"

        l_errors: list[str] = []

        # Check the exception type
        if self.ex_type and not issubclass(exc_info.type, self.ex_type):
            l_errors.append(f"ERROR: Raised exception is of type '{exc_info.type}', but expected '{self.ex_type}'. "
                            f"Raised exception's message was: {str(exc_info)}")

        exc = exc_info.value

        # Check the exception message, if applicable
        if self.ex_message:
            if len(exc.args) == 0:
                l_errors.append(f"ERROR: Expected string \"{self.ex_message}\" in exception of type '{type(exc)}''s "
                                "message, but exception does not have any message")
            elif not string_with_placeholders_matches(self.ex_message, str(exc.args[0])):
                l_errors.append(f"ERROR: Expected string \"{self.ex_message}\" not found in exception message: "
                                f"\"{exc.args[0]}\"")

        # Check the exception status code, if applicable
        if self.ex_status_code:
            if not hasattr(exc, "status_code"):
                l_errors.append(f"ERROR: Expected status code {self.ex_status_code} for exception of type {type(exc)}, "
                                "but this exception type does not have a status code attribute")
            elif self.ex_status_code != exc.status_code:
                l_errors.append(f"ERROR: Expected status code {self.ex_status_code} does not match that of {exc}: "
                                f"{exc.status_code}")

        # Join any errors for output
        res = "\n".join(l_errors)
        return res


@dataclass
class MatchOutputFile:
    """Callable class which checks that the output file of a conversion numerically matches an expected output file"""

    ex_output_filename: str
    """The name of the file which the output of the conversion should match, relative to the test_data/output directory
    """

    def __call__(self, test_info: ConversionTestInfo):
        """Run the check comparing the two files"""

        qualified_ex_output_filename = os.path.join(OUTPUT_TEST_DATA_LOC_IN_PROJECT, self.ex_output_filename)

        return check_file_match(test_info.qualified_out_filename, qualified_ex_output_filename)
