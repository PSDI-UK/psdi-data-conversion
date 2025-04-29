"""
# gui.py

Utilities to aid in testing of the GUI
"""

import os
import shutil
import time
from dataclasses import dataclass
from tempfile import TemporaryDirectory

import pytest
from selenium.common.exceptions import TimeoutException
from selenium.webdriver.common.alert import Alert
from selenium.webdriver.common.by import By
from selenium.webdriver.firefox.webdriver import WebDriver
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.select import Select
from selenium.webdriver.support.ui import WebDriverWait

from psdi_data_conversion.constants import STATUS_CODE_GENERAL
from psdi_data_conversion.converters.base import (FileConverterAbortException, FileConverterException,
                                                  FileConverterInputException)
from psdi_data_conversion.converters.openbabel import (COORD_GEN_KEY, COORD_GEN_QUAL_KEY, DEFAULT_COORD_GEN_QUAL,
                                                       L_ALLOWED_COORD_GEN_QUALS, L_ALLOWED_COORD_GENS)
from psdi_data_conversion.database import get_format_info
from psdi_data_conversion.file_io import split_archive_ext
from psdi_data_conversion.testing.constants import DEFAULT_ORIGIN
from psdi_data_conversion.testing.utils import (ConversionTestInfo, ConversionTestSpec, SingleConversionTestSpec,
                                                get_input_test_data_loc)

# Standard timeout at 10 seconds
TIMEOUT = 10


def wait_for_element(root: WebDriver | EC.WebElement, xpath: str, by=By.XPATH):
    """Shortcut for boilerplate to wait until a web element is visible"""
    WebDriverWait(root, TIMEOUT).until(EC.element_to_be_clickable((by, xpath)))


def wait_and_find_element(root: WebDriver | EC.WebElement, xpath: str, by=By.XPATH) -> EC.WebElement:
    """Finds a web element, after first waiting to ensure it's visible"""
    wait_for_element(root, xpath, by=by)
    return root.find_element(by, xpath)


def wait_for_cover_hidden(root: WebDriver):
    """Wait until the page cover is removed"""
    WebDriverWait(root, TIMEOUT).until(EC.invisibility_of_element((By.XPATH, "//div[@id='cover']")))


@dataclass
class GuiTestSpecRunner():
    """Class which provides an interface to run test conversions through the GUI
    """

    driver: WebDriver
    """The WebDriver to be used for testing"""

    origin: str = DEFAULT_ORIGIN
    """The address of the homepage of the testing server"""

    def run(self, test_spec: ConversionTestSpec):
        """Run the test conversions outlined in a test spec"""

        self._test_spec = test_spec

        # Make temporary directories for the input and output files to be stored in
        with TemporaryDirectory("_input") as input_dir, TemporaryDirectory("_output") as output_dir:

            # Iterate over the test spec to run each individual test it defines
            for single_test_spec in test_spec:
                if single_test_spec.skip:
                    print(f"Skipping single test spec {single_test_spec}")
                    continue

                print(f"Running single test spec: {single_test_spec}")

                GuiSingleTestSpecRunner(parent=self,
                                        input_dir=input_dir,
                                        output_dir=output_dir,
                                        single_test_spec=single_test_spec).run()

                print(f"Success for test spec: {single_test_spec}")


class GuiSingleTestSpecRunner:
    """Class which handles running an individual test conversion
    """

    def __init__(self,
                 parent: GuiTestSpecRunner,
                 input_dir: str,
                 output_dir: str,
                 single_test_spec: SingleConversionTestSpec):
        """

        Parameters
        ----------
        parent : GuiTestSpecRunner
            The GuiTestSpecRunner which created this and is running it
        input_dir : str
            The temporary directory to be used for input data
        output_dir : str
            The temporary directory to be used for output data
        single_test_spec : SingleConversionTestSpec
            The test spec that is currently being tested
        """

        self.input_dir: str = input_dir
        self.output_dir: str = output_dir
        self.single_test_spec: SingleConversionTestSpec = single_test_spec

        # Inherit data from the parent class

        self.driver: WebDriver = parent.driver
        """The WebDriver to be used for testing"""

        self.origin: str = parent.origin
        """The address of the homepage of the testing server"""

        # Interpret information from the test spec that we'll need for testing

        # Get just the local filename
        self._filename = os.path.split(self.single_test_spec.filename)[1]

        # Default options for conversion
        self._base_filename, ext = split_archive_ext(self._filename)
        self._strict = True
        self._from_flags: str | None = None
        self._to_flags: str | None = None
        self._from_options: str | None = None
        self._to_options: str | None = None
        self._coord_gen = None
        self._coord_gen_qual = None

        # Get the from_format from the extension if not provided
        from_format = single_test_spec.from_format
        if not from_format:
            from_format = ext

        # Get the format info for each format, which we'll use to get the name and note of each
        self._from_format_info = get_format_info(from_format, which=0)
        self._to_format_info = get_format_info(single_test_spec.to_format, which=0)

        # For each argument in the conversion kwargs, interpret it as the appropriate option for this conversion,
        # overriding defaults set above
        for key, val in self.single_test_spec.conversion_kwargs.items():
            if key == "log_mode":
                raise ValueError(f"The conversion kwarg {key} is not valid with conversions through the GUI")
            elif key == "delete_input":
                raise ValueError(f"The conversion kwarg {key} is not valid with conversions through the GUI")
            elif key == "strict":
                self._strict = val
            elif key == "max_file_size":
                raise ValueError(f"The conversion kwarg {key} is not valid with conversions through the GUI")
            elif key == "data":
                for subkey, subval in val.items():
                    if subkey == "from_flags":
                        self._from_flags = subval
                    elif subkey == "to_flags":
                        self._to_flags = subval
                    elif subkey == "from_options":
                        self._from_options = subval
                    elif subkey == "to_options":
                        self._to_options = subval
                    elif subkey == COORD_GEN_KEY:
                        self._coord_gen = subval
                        if COORD_GEN_QUAL_KEY in val:
                            self._coord_gen_qual = val[COORD_GEN_QUAL_KEY]
                        else:
                            self._coord_gen_qual = DEFAULT_COORD_GEN_QUAL
                    elif subkey == COORD_GEN_QUAL_KEY:
                        # Handled alongside COORD_GEN_KEY above
                        pass
                    else:
                        pytest.fail(f"The key 'data[\"{subkey}\"]' was passed to `conversion_kwargs` but could not be "
                                    "interpreted")
            else:
                pytest.fail(f"The key '{key}' was passed to `conversion_kwargs` but could not be interpreted")

    def run(self):
        """Run the conversion outlined in the test spec"""

        exc_info: pytest.ExceptionInfo | None = None
        if self.single_test_spec.expect_success:
            try:
                self._run_conversion()
                success = False
            except Exception:
                print(f"Unexpected exception raised for single test spec {self.single_test_spec}")
                raise
        else:
            with pytest.raises(FileConverterException) as exc_info:
                self._run_conversion()
            success = False

        # Compile output info for the test and call the callback function if one is provided
        if self.single_test_spec.callback:
            test_info = ConversionTestInfo(run_type="gui",
                                           test_spec=self.single_test_spec,
                                           input_dir=self.input_dir,
                                           output_dir=self.output_dir,
                                           success=success,
                                           exc_info=exc_info)
            callback_msg = self.single_test_spec.callback(test_info)
            if callback_msg:
                pytest.fail(callback_msg)

    def _run_conversion(self):
        """Run a conversion through the GUI
        """

        self._set_up_files()

        self._select_formats_and_converter()

        self._set_conversion_settings()

        self._provide_input_file()

        self._request_conversion()

        self._move_output()

    def _set_up_files(self):
        """Set up the filenames we expect and initialize them - delete any leftover files and symlink the input file
        to the desired location
        """
        # Set up the expected filenames
        source_input_file = os.path.realpath(os.path.join(get_input_test_data_loc(), self.single_test_spec.filename))
        self._input_file = os.path.join(self.input_dir, self.single_test_spec.filename)

        self._log_file = os.path.realpath(os.path.join(os.path.expanduser("~/Downloads"),
                                                       self.single_test_spec.log_filename))

        self._output_file = os.path.realpath(os.path.join(os.path.expanduser("~/Downloads"),
                                                          self.single_test_spec.out_filename))

        # Clean up any leftover files
        if (os.path.isfile(self._input_file)):
            os.unlink(self._input_file)
        if (os.path.isfile(self._log_file)):
            os.remove(self._log_file)
        if (os.path.isfile(self._output_file)):
            os.remove(self._output_file)

        # Symlink the input file to the desired location
        os.symlink(source_input_file, self._input_file)

    def _select_formats_and_converter(self):
        """Handle the tasks on the format and converter selection page when running a test:

        1. Load the main page (waiting for it to fully load)
        2. Select the input and output formats
        3. Select the converter
        4. Click the "Yes" button to confirm and go to the convert page
        """

        # Get the homepage and wait for the cover to be removed
        self.driver.get(f"{self.origin}/")
        wait_for_cover_hidden(self.driver)

        wait_for_element(self.driver, "//select[@id='fromList']/option")

        # Select from_format from the 'from' list.
        full_from_format = f"{self._from_format_info.name}: {self._from_format_info.note}"
        self.driver.find_element(
            By.XPATH, f"//select[@id='fromList']/option[starts-with(.,'{full_from_format}')]").click()

        # Select to_format from the 'to' list.
        full_to_format = f"{self._to_format_info.name}: {self._to_format_info.note}"
        self.driver.find_element(
            By.XPATH, f"//select[@id='toList']/option[starts-with(.,'{full_to_format}')]").click()

        # Select converter from the available conversion options list.
        self.driver.find_element(
            By.XPATH, f"//select[@id='success']/option[contains(.,'{self.single_test_spec.converter_name}')]").click()

        # Click on the "Yes" button to accept the converter and go to the conversion page, and wait for the cover to be
        # removed there
        self.driver.find_element(By.XPATH, "//input[@id='yesButton']").click()
        wait_for_cover_hidden(self.driver)

    def _set_conversion_settings(self):
        """Set settings on the convert page appropriately for the desired conversion
        """
        # Request non-strict filename checking if desired
        if not self._strict:
            wait_and_find_element(self.driver, "//input[@id='extCheck']").click()

        # Request the log file
        wait_and_find_element(self.driver, "//input[@id='requestLog']").click()

        # Set appropriate format and converter settings for this conversion
        self._select_format_flags()
        self._set_format_options()
        self._apply_radio_settings()

    def _provide_input_file(self):
        """Provide the input file for the conversion, checking if any alert is raised in response
        """
        # Select the input file
        wait_and_find_element(self.driver, "//input[@id='fileToUpload']").send_keys(str(self._input_file))

        # An alert may be present here, which we check for using a try block
        try:
            WebDriverWait(self.driver, 0.2).until(EC.alert_is_present())
            alert = Alert(self.driver)
            alert_text = alert.text
            alert.dismiss()
            raise FileConverterInputException(alert_text)
        except TimeoutException:
            pass

    def _select_format_flags(self):
        """Select desired format flags. The options in the select box only have a text attribute, so we need to find
        the one that starts with each flag - since we don't have too many, iterating over all possible combinations is
        the easiest way
        """
        for (l_flags, select_id) in ((self._from_flags, "inFlags"),
                                     (self._to_flags, "outFlags")):
            if not l_flags:
                continue
            flags_select = Select(wait_and_find_element(self.driver, f"//select[@id='{select_id}']"))
            for flag in l_flags:
                for option in flags_select.options:
                    if option.text.startswith(f"{flag}:"):
                        flags_select.select_by_visible_text(option.text)
                        break
                else:
                    raise ValueError(f"Flag {flag} was not found in {select_id} selection box for conversion from "
                                     f"{self._from_format_info.name} to {self._to_format_info.name} with "
                                     f"converter {self.single_test_spec.converter_name}")

    def _set_format_options(self):
        """Set desired format options
        """
        for (options_string, table_id) in ((self._from_options, "in_argFlags"),
                                           (self._to_options, "out_argFlags")):
            if not options_string:
                continue

            # Split each option into words, of which the first letter of each is the key and the remainder is the value
            l_options = options_string.split()

            # Get the rows in the options table
            options_table = wait_and_find_element(self.driver, f"//table[@id='{table_id}']")
            l_rows = options_table.find_elements(By.XPATH, "./tr")

            # Look for and set each option
            for option in l_options:
                for row in l_rows:
                    l_items = row.find_elements(By.XPATH, "./td")
                    label = l_items[1]
                    if not label.text.startswith(option[0]):
                        continue

                    # Select the option by clicking the box at the first element in the row to make the input appear
                    l_items[0].click()

                    # Input the option in the input box that appears in the third position in the row
                    input_box = wait_and_find_element(l_items[2], "./input")
                    input_box.send_keys(option[1:])

                    break

                else:
                    raise ValueError(f"Option {option} was not found in {table_id} options table for conversion from "
                                     f"{self._from_format_info.name} to {self._to_format_info.name} with "
                                     f"converter {self.single_test_spec.converter_name}")

    def _apply_radio_settings(self):
        """Apply any radio-button settings desired for this conversion by clicking the appropriate radio buttons
        """

        for setting, name, l_allowed in ((self._coord_gen, "coord_gen", L_ALLOWED_COORD_GENS),
                                         (self._coord_gen_qual, "coord_gen_qual", L_ALLOWED_COORD_GEN_QUALS)):
            if not setting:
                continue

            if setting not in l_allowed:
                raise ValueError(f"Invalid {name} value supplied: {setting}. Allowed values are: " +
                                 str(l_allowed))

            setting_radio = wait_and_find_element(self.driver, f"//input[@value='{setting}']")
            setting_radio.click()

    def _request_conversion(self):
        """Request the conversion, handle the alert box that appears, and wait for the files to be downloaded
        """
        # Click on the "Convert" button.
        wait_and_find_element(self.driver, "//input[@id='uploadButton']").click()

        # Handle alert box.
        WebDriverWait(self.driver, TIMEOUT).until(EC.alert_is_present())
        alert = Alert(self.driver)
        alert_text = alert.text
        alert.dismiss()

        if alert_text.startswith("ERROR:"):
            # Raise an appropriate exception type depending on if it's a recognised input issue or not
            if "unexpected exception" in alert_text:
                raise FileConverterAbortException(STATUS_CODE_GENERAL, alert_text)
            raise FileConverterInputException(alert_text)

        # Wait until the log file exists, since it's downloaded second
        time_elapsed = 0
        while not os.path.isfile(self._log_file):
            time.sleep(1)
            time_elapsed += 1
            if time_elapsed > TIMEOUT:
                pytest.fail(f"Download of {self._output_file} and {self._log_file} timed out")

        time.sleep(1)

    def _move_output(self):
        """Move the created output files out of the default Downloads directory and into the desired output files
        directory.
        """
        # Check for the presence of the output file
        if not os.path.isfile(self._output_file):
            raise FileConverterAbortException("ERROR: No output file was produced. Log contents:\n" +
                                              open(self._log_file, "r").read())

        # Move the output file and log file to the expected locations
        for qual_filename in self._output_file, self._log_file:
            self._base_filename = os.path.split(qual_filename)[1]
            target_filename = os.path.join(self.output_dir, self._base_filename)
            if os.path.isfile(target_filename):
                os.remove(target_filename)
            if os.path.isfile(qual_filename):
                shutil.move(qual_filename, target_filename)
