"""
# gui.py

Utilities to aid in testing of the GUI
"""


import os
import shutil
from tempfile import TemporaryDirectory

import time
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
from psdi_data_conversion.file_io import split_archive_ext
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


def run_test_conversion_with_gui(test_spec: ConversionTestSpec,
                                 driver: WebDriver,
                                 origin: str):
    """Runs a test conversion or series thereof through the GUI. Note that this requires the server to be started before
    this is called.

    Parameters
    ----------
    test_spec : ConversionTestSpec
        The specification for the test or series of tests to be run
    driver : WebDriver
        The WebDriver to be used for testing
    origin : str
        The address of the homepage of the testing server
    """
    # Make temporary directories for the input and output files to be stored in
    with TemporaryDirectory("_input") as input_dir, TemporaryDirectory("_output") as output_dir:
        # Iterate over the test spec to run each individual test it defines
        for single_test_spec in test_spec:
            if single_test_spec.skip:
                print(f"Skipping single test spec {single_test_spec}")
                continue
            print(f"Running single test spec: {single_test_spec}")
            _run_single_test_conversion_with_gui(test_spec=single_test_spec,
                                                 input_dir=input_dir,
                                                 output_dir=output_dir,
                                                 driver=driver,
                                                 origin=origin)
            print(f"Success for test spec: {single_test_spec}")


def _run_single_test_conversion_with_gui(test_spec: SingleConversionTestSpec,
                                         input_dir: str,
                                         output_dir: str,
                                         driver: WebDriver,
                                         origin: str):
    """Runs a single test conversion through the GUI.

    Parameters
    ----------
    test_spec : SingleConversionTestSpec
        The specification for the test to be run
    input_dir : str
        A directory which can be used to store input data before uploading
    output_dir : str
        A directory which can be used to create output data after downloading
    driver : WebDriver
        The WebDriver to be used for testing
    origin : str
        The address of the homepage of the testing server
    """

    exc_info: pytest.ExceptionInfo | None = None
    if test_spec.expect_success:
        try:
            run_converter_through_gui(test_spec=test_spec,
                                      input_dir=input_dir,
                                      output_dir=output_dir,
                                      driver=driver,
                                      origin=origin,
                                      **test_spec.conversion_kwargs)
            success = True
        except Exception:
            print(f"Unexpected exception raised for single test spec {test_spec}")
            raise
    else:
        with pytest.raises(FileConverterException) as exc_info:
            run_converter_through_gui(test_spec=test_spec,
                                      input_dir=input_dir,
                                      output_dir=output_dir,
                                      driver=driver,
                                      origin=origin,
                                      **test_spec.conversion_kwargs)
        success = False

    # Compile output info for the test and call the callback function if one is provided
    if test_spec.callback:
        test_info = ConversionTestInfo(run_type="gui",
                                       test_spec=test_spec,
                                       input_dir=input_dir,
                                       output_dir=output_dir,
                                       success=success,
                                       exc_info=exc_info)
        callback_msg = test_spec.callback(test_info)
        assert not callback_msg, callback_msg


def run_converter_through_gui(test_spec: SingleConversionTestSpec,
                              input_dir: str,
                              output_dir: str,
                              driver: WebDriver,
                              origin: str,
                              **conversion_kwargs):
    """_summary_

    Parameters
    ----------
    test_spec : SingleConversionTestSpec
        The specification for the test to be run
    input_dir : str
        The directory which contains the input file
    output_dir : str
        The directory which contains the output file
    driver : WebDriver
        The WebDriver to be used for testing
    origin : str
        The address of the homepage of the testing server
    """

    # Get just the local filename
    filename = os.path.split(test_spec.filename)[1]

    # Default options for conversion
    base_filename, from_format = split_archive_ext(filename)
    strict = True
    from_flags: str | None = None
    to_flags: str | None = None
    from_options: str | None = None
    to_options: str | None = None
    coord_gen = None
    coord_gen_qual = None

    # For each argument in the conversion kwargs, interpret it as the appropriate option for this conversion, overriding
    # defaults set above
    for key, val in conversion_kwargs.items():
        if key == "from_format":
            from_format = val
        elif key == "log_mode":
            raise ValueError(f"The conversion kwarg {key} is not valid with conversions through the GUI")
        elif key == "delete_input":
            raise ValueError(f"The conversion kwarg {key} is not valid with conversions through the GUI")
        elif key == "strict":
            strict = val
        elif key == "max_file_size":
            raise ValueError(f"The conversion kwarg {key} is not valid with conversions through the GUI")
        elif key == "data":
            for subkey, subval in val.items():
                if subkey == "from_flags":
                    from_flags = subval
                elif subkey == "to_flags":
                    to_flags = subval
                elif subkey == "from_options":
                    from_options = subval
                elif subkey == "to_options":
                    to_options = subval
                elif subkey == COORD_GEN_KEY:
                    coord_gen = subval
                    if COORD_GEN_QUAL_KEY in val:
                        coord_gen_qual = val[COORD_GEN_QUAL_KEY]
                    else:
                        coord_gen_qual = DEFAULT_COORD_GEN_QUAL
                elif subkey == COORD_GEN_QUAL_KEY:
                    # Handled alongside COORD_GEN_KEY above
                    pass
                else:
                    pytest.fail(f"The key 'data[\"{subkey}\"]' was passed to `conversion_kwargs` but could not be "
                                "interpreted")
        else:
            pytest.fail(f"The key '{key}' was passed to `conversion_kwargs` but could not be interpreted")

    # Cleanup of arguments
    if from_format.startswith("."):
        from_format = from_format[1:]

    # Set up the input file where we expect it to be
    source_input_file = os.path.realpath(os.path.join(get_input_test_data_loc(), test_spec.filename))
    input_file = os.path.join(input_dir, test_spec.filename)
    if (os.path.isfile(input_file)):
        os.unlink(input_file)
    os.symlink(source_input_file, input_file)

    # Remove test files from Downloads directory if they exist.

    log_file = os.path.realpath(os.path.join(os.path.expanduser("~/Downloads"), test_spec.log_filename))
    if (os.path.isfile(log_file)):
        os.remove(log_file)

    output_file = os.path.realpath(os.path.join(os.path.expanduser("~/Downloads"), test_spec.out_filename))
    if (os.path.isfile(output_file)):
        os.remove(output_file)

    # Get the homepage
    driver.get(f"{origin}/")

    wait_for_element(driver, "//select[@id='fromList']/option")

    # Select from_format from the 'from' list.
    driver.find_element(By.XPATH, f"//select[@id='fromList']/option[starts-with(.,'{from_format}:')]").click()

    # Select to_format from the 'to' list.
    driver.find_element(By.XPATH, f"//select[@id='toList']/option[starts-with(.,'{test_spec.to_format}:')]").click()

    # Select converter from the available conversion options list.
    driver.find_element(By.XPATH, f"//select[@id='success']/option[contains(.,'{test_spec.converter_name}')]").click()

    # Click on the "Yes" button to accept the converter and go to the conversion page
    driver.find_element(By.XPATH, "//input[@id='yesButton']").click()

    # Request non-strict filename checking if desired
    if not strict:
        wait_and_find_element(driver, "//input[@id='extCheck']").click()

    # Select the input file
    wait_and_find_element(driver, "//input[@id='fileToUpload']").send_keys(str(input_file))

    # An alert may be present here, which we check for using a try block
    try:
        WebDriverWait(driver, 0.2).until(EC.alert_is_present())
        alert = Alert(driver)
        alert_text = alert.text
        alert.dismiss()
        raise FileConverterInputException(alert_text)
    except TimeoutException:
        pass

    # Request the log file
    wait_and_find_element(driver, "//input[@id='requestLog']").click()

    # Select appropriate format args. The args only have a text attribute, so we need to find the one that starts with
    # each flag - since we don't have too many, iterating over all possible combinations is the easiest way
    for (l_flags, select_id) in ((from_flags, "inFlags"),
                                 (to_flags, "outFlags")):
        if not l_flags:
            continue
        flags_select = Select(wait_and_find_element(driver, f"//select[@id='{select_id}']"))
        for flag in l_flags:
            found = False
            for option in flags_select.options:
                if option.text.startswith(f"{flag}:"):
                    flags_select.select_by_visible_text(option.text)
                    found = True
                    break
            if not found:
                raise ValueError(f"Flag {flag} was not found in {select_id} selection box for conversion from "
                                 f"{from_format} to {test_spec.to_format} with converter {test_spec.converter_name}")

    for (options_string, table_id) in ((from_options, "in_argFlags"),
                                       (to_options, "out_argFlags")):
        if not options_string:
            continue

        # Split each option into words, of which the first letter of each is the key and the remainder is the value
        l_options = options_string.split()

        # Get the rows in the options table
        options_table = wait_and_find_element(driver, f"//table[@id='{table_id}']")
        l_rows = options_table.find_elements(By.XPATH, "./tr")

        # Look for and set each option
        for option in l_options:
            found = False
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

                found = True
                break

            if not found:
                raise ValueError(f"Option {option} was not found in {table_id} options table for conversion from "
                                 f"{from_format} to {test_spec.to_format} with converter {test_spec.converter_name}")

    # If radio-button settings are supplied, apply them now
    for setting, name, l_allowed in ((coord_gen, "coord_gen", L_ALLOWED_COORD_GENS),
                                     (coord_gen_qual, "coord_gen_qual", L_ALLOWED_COORD_GEN_QUALS)):
        if not setting:
            continue

        if setting not in l_allowed:
            raise ValueError(f"Invalid {name} value supplied: {setting}. Allowed values are: " +
                             str(l_allowed))

        setting_radio = wait_and_find_element(driver, f"//input[@value='{setting}']")
        setting_radio.click()

    # Click on the "Convert" button.
    wait_and_find_element(driver, "//input[@id='uploadButton']").click()

    # Handle alert box.
    WebDriverWait(driver, TIMEOUT).until(EC.alert_is_present())
    alert = Alert(driver)
    alert_text = alert.text
    alert.dismiss()

    if alert_text.startswith("ERROR:"):
        # Raise an appropriate exception type depending on if it's a recognised input issue or not
        if "unexpected exception" in alert_text:
            raise FileConverterAbortException(STATUS_CODE_GENERAL, alert_text)
        raise FileConverterInputException(alert_text)

    # Wait until the log file exists, since it's downloaded second
    time_elapsed = 0
    while not os.path.isfile(log_file):
        time.sleep(1)
        time_elapsed += 1
        if time_elapsed > TIMEOUT:
            pytest.fail(f"Download of {output_file} and {log_file} timed out")

    time.sleep(1)

    if not os.path.isfile(output_file):
        raise FileConverterAbortException("ERROR: No output file was produced. Log contents:\n" +
                                          open(log_file, "r").read())

    # Move the output file and log file to the expected locations
    for qual_filename in output_file, log_file:
        base_filename = os.path.split(qual_filename)[1]
        target_filename = os.path.join(output_dir, base_filename)
        if os.path.isfile(target_filename):
            os.remove(target_filename)
        if os.path.isfile(qual_filename):
            shutil.move(qual_filename, target_filename)
