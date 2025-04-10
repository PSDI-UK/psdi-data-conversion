"""
# gui.py

Utilities to aid in testing of the GUI
"""


from dataclasses import dataclass
import os
import shutil
from tempfile import TemporaryDirectory

from pathlib import Path
import time
import pytest
from selenium.webdriver.common.alert import Alert
from selenium.webdriver.common.by import By
from selenium.webdriver.firefox.webdriver import WebDriver
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.ui import WebDriverWait

from psdi_data_conversion.converters.openbabel import (COORD_GEN_KEY, COORD_GEN_QUAL_KEY, DEFAULT_COORD_GEN,
                                                       DEFAULT_COORD_GEN_QUAL)
from psdi_data_conversion.file_io import split_archive_ext
from psdi_data_conversion.testing.utils import (ConversionTestInfo, ConversionTestSpec, SingleConversionTestSpec,
                                                get_input_test_data_loc)

# Standard timeout at 10 seconds
TIMEOUT = 10


def wait_for_element(driver: WebDriver, xpath: str, by=By.XPATH):
    """Shortcut for boilerplate to wait until a web element is visible"""
    WebDriverWait(driver, TIMEOUT).until(EC.element_to_be_clickable((by, xpath)))


def wait_and_find_element(driver: WebDriver, xpath: str, by=By.XPATH) -> EC.WebElement:
    """Finds a web element, after first waiting to ensure it's visible"""
    wait_for_element(driver, xpath, by=by)
    return driver.find_element(by, xpath)


@dataclass
class GuiConversionTestInfo(ConversionTestInfo):
    """Information about a tested conversion, specifically for when it was tested through the GUI (the local version of
    the web app)
    """


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
    test_spec : _SingleConversionTestSpec
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

    # Symlink the input file to the input directory
    qualified_in_filename = os.path.realpath(os.path.join(input_dir, test_spec.filename))
    try:
        os.symlink(os.path.join(get_input_test_data_loc(), test_spec.filename),
                   qualified_in_filename)
    except FileExistsError:
        pass

    try:
        success = run_converter_through_gui(filename=qualified_in_filename,
                                            to_format=test_spec.to_format,
                                            name=test_spec.converter_name,
                                            input_dir=input_dir,
                                            output_dir=output_dir,
                                            log_file=os.path.join(output_dir, test_spec.log_filename),
                                            driver=driver,
                                            origin=origin,
                                            **test_spec.conversion_kwargs)
    except Exception:
        print(f"Unexpected exception raised for single test spec {test_spec}")
        raise

    # Compile output info for the test and call the callback function if one is provided
    if test_spec.callback:
        test_info = GuiConversionTestInfo(test_spec=test_spec,
                                          input_dir=input_dir,
                                          output_dir=output_dir,
                                          success=success)
        callback_msg = test_spec.callback(test_info)
        assert not callback_msg, callback_msg


def run_converter_through_gui(filename: str,
                              to_format: str,
                              name: str,
                              input_dir: str,
                              output_dir: str,
                              log_file: str,
                              driver: WebDriver,
                              origin: str,
                              **conversion_kwargs):
    """_summary_

    Parameters
    ----------
    filename : str
        The (unqualified) name of the input file to be converted
    to_format : str
        The format to convert the input file to
    name : str
        The name of the converter to use
    input_dir : str
        The directory which contains the input file
    output_dir : str
        The directory which contains the output file
    log_file : str
        The desired name of the log file
    driver : WebDriver
        The WebDriver to be used for testing
    origin : str
        The address of the homepage of the testing server
    """

    # Get just the local filename
    filename = os.path.split(filename)[1]

    # Default options for conversion
    base_filename, from_format = split_archive_ext(filename)
    strict = True
    from_flags = None
    to_flags = None
    from_options = None
    to_options = None
    coord_gen = DEFAULT_COORD_GEN
    coord_gen_qual = DEFAULT_COORD_GEN_QUAL

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

    input_file = os.path.realpath(os.path.join(input_dir, filename))
    output_file = Path.home().joinpath("Downloads", f"{base_filename}.{to_format}")
    log_file = Path.home().joinpath("Downloads", f"{base_filename}.log.txt")

    # Remove test files from Downloads directory if they exist.

    if (Path.is_file(log_file)):
        Path.unlink(log_file)

    if (Path.is_file(output_file)):
        Path.unlink(output_file)

    # Get the homepage
    driver.get(f"{origin}/")

    wait_for_element(driver, "//select[@id='fromList']/option")

    # Select cdxml from the 'from' list.
    driver.find_element(By.XPATH, f"//select[@id='fromList']/option[contains(.,'{from_format}:')]").click()

    # Select InChI from the 'to' list.
    driver.find_element(By.XPATH, f"//select[@id='toList']/option[contains(.,'{to_format}:')]").click()

    # Select Open Babel from the available conversion options list.
    driver.find_element(By.XPATH, f"//select[@id='success']/option[contains(.,'{name}')]").click()

    # Click on the "Yes" button to accept the converter and go to the conversion page
    driver.find_element(By.XPATH, "//input[@id='yesButton']").click()

    # Select the input file.
    wait_and_find_element(driver, "//input[@id='fileToUpload']").send_keys(str(input_file))

    # Request the log file
    wait_and_find_element(driver, "//input[@id='requestLog']").click()

    # Request non-strict filename checking if desired
    if not strict:
        wait_and_find_element(driver, "//input[@id='extCheck']").click()

    # Click on the "Convert" button.
    wait_and_find_element(driver, "//input[@id='uploadButton']").click()

    # Handle alert box.
    WebDriverWait(driver, TIMEOUT).until(EC.alert_is_present())
    Alert(driver).dismiss()

    # Wait until the log file exists - on failure, the output file won't exist, but the log file always will
    time_elapsed = 0
    while not Path.is_file(log_file):
        time.sleep(1)
        time_elapsed += 1
        if time_elapsed > TIMEOUT:
            pytest.fail(f"Download of {output_file} and {log_file} timed out")

    time.sleep(1)

    if not Path.is_file(output_file):
        success = False
    else:
        success = True

    # Move the output file and log file to the expected locations
    shutil.move(output_file, output_dir)
    shutil.move(log_file, output_dir)

    return success
