#!/usr/bin/env python

# Selenium test script for PSDI Data Conversion Service.

import os
import time
from multiprocessing import Process

from pathlib import Path
import pytest

from psdi_data_conversion.testing.gui import (TIMEOUT, run_test_conversion_with_gui, wait_and_find_element,
                                              wait_for_element)
from psdi_data_conversion.testing.conversion_test_specs import l_gui_test_specs

# Skip all tests in this module if required packages for GUI testing aren't installed
try:
    from selenium.common.exceptions import NoSuchElementException
    from selenium import webdriver
    from selenium.webdriver.common.alert import Alert
    from selenium.webdriver.common.by import By
    from selenium.webdriver import FirefoxOptions
    from selenium.webdriver.firefox.service import Service as FirefoxService
    from selenium.webdriver.firefox.webdriver import WebDriver
    from selenium.webdriver.support import expected_conditions as EC
    from selenium.webdriver.support.ui import WebDriverWait
    from webdriver_manager.firefox import GeckoDriverManager

    from psdi_data_conversion.app import start_app

except ImportError:
    # We put the importorskip commands here rather than above so that standard imports can be used by static analysis
    # tools where possible, and the importorskip is used here so pytest will stop processing immediately if things can't
    # be imported - pytest.mark.skip won't do that
    pytest.importorskip("Flask")
    pytest.importorskip("selenium")
    pytest.importorskip("webdriver_manager.firefox")


import psdi_data_conversion
from psdi_data_conversion.testing.utils import get_input_test_data_loc

origin = os.environ.get("ORIGIN", "http://127.0.0.1:5000")


@pytest.fixture(scope="module", autouse=True)
def common_setup():
    """Autouse fixture which starts the app before tests and stops it afterwards"""

    server = Process(target=start_app)
    server.start()

    # Change to the root dir of the project for running the tests, in case this was invoked elsewhere
    old_cwd = os.getcwd()
    os.chdir(os.path.join(psdi_data_conversion.__path__[0], ".."))

    yield

    server.terminate()
    server.join()

    # Change back to the previous directory
    os.chdir(old_cwd)


@pytest.fixture(scope="module")
def driver():
    """Get a headless Firefox web driver"""

    driver_path = os.environ.get("DRIVER")

    if not driver_path:
        driver_path = GeckoDriverManager().install()
        print(f"Gecko driver installed to {driver_path}")

    opts = FirefoxOptions()
    opts.add_argument("--headless")
    ff_driver = webdriver.Firefox(service=FirefoxService(driver_path),
                                  options=opts)
    yield ff_driver
    ff_driver.quit()


@pytest.fixture(scope="module", autouse=True)
def setup(driver: WebDriver):
    """Run common tasks for each test"""

    driver.get(f"{origin}/")


def test_initial_frontpage(driver: WebDriver):

    # Check that the front page contains the header "Data Conversion Service".

    element = wait_and_find_element(driver, "//header//h5")
    assert element.text == "Data Conversion Service"

    # Check that the 'from' and 'to' lists contains "abinit" and "acesin" respectively.

    wait_for_element(driver, "//select[@id='fromList']/option")
    driver.find_element(By.XPATH, "//select[@id='fromList']/option[contains(.,'abinit: ABINIT output')]")

    wait_for_element(driver, "//select[@id='toList']/option")
    driver.find_element(By.XPATH, "//select[@id='toList']/option[contains(.,'acesin: ACES input')]")

    # Check that the available conversions list is empty.

    with pytest.raises(NoSuchElementException):
        driver.find_element(By.XPATH, "//select[@id='success']/option")


def test_cdxml_to_inchi_conversion(driver: WebDriver):

    test_file = "standard_test"

    input_file = os.path.realpath(os.path.join(get_input_test_data_loc(), f"{test_file}.cdxml"))
    output_file = Path.home().joinpath("Downloads", f"{test_file}.inchi")
    log_file = Path.home().joinpath("Downloads", f"{test_file}.log.txt")

    # Remove test files from Downloads directory if they exist.

    if (Path.is_file(log_file)):
        Path.unlink(log_file)

    if (Path.is_file(output_file)):
        Path.unlink(output_file)

    wait_for_element(driver, "//select[@id='fromList']/option")

    # Select cdxml from the 'from' list.
    driver.find_element(By.XPATH, "//select[@id='fromList']/option[contains(.,'cdxml: ChemDraw CDXML')]").click()

    # Select InChI from the 'to' list.
    driver.find_element(By.XPATH, "//select[@id='toList']/option[contains(.,'inchi: InChI')]").click()

    # Select Open Babel from the available conversion options list.
    driver.find_element(By.XPATH, "//select[@id='success']/option[contains(.,'Open Babel')]").click()

    # Click on the "Yes" button to accept the converter and go to the conversion page
    driver.find_element(By.XPATH, "//input[@id='yesButton']").click()

    # Select the input file.
    wait_and_find_element(driver, "//input[@id='fileToUpload']").send_keys(str(input_file))

    # Request the log file
    wait_and_find_element(driver, "//input[@id='requestLog']").click()

    # Click on the "Convert" button.
    wait_and_find_element(driver, "//input[@id='uploadButton']").click()

    # Handle alert box.
    WebDriverWait(driver, TIMEOUT).until(EC.alert_is_present())
    Alert(driver).dismiss()

    # Wait until files exist.

    time_elapsed = 0
    while (not Path.is_file(log_file)) or (not Path.is_file(output_file)):
        time.sleep(1)
        time_elapsed += 1
        if time_elapsed > TIMEOUT:
            pytest.fail(f"Download of {output_file} and {log_file} timed out")

    time.sleep(1)

    # Verify that the InChI file is correct.
    assert output_file.read_text().strip() == "InChI=1S/C12NO/c1-12(2)6-7-13-11-5-4-9(14-3)8-10(11)12"


@pytest.mark.parametrize("test_spec", l_gui_test_specs,
                         ids=lambda x: x.name)
def test_conversions(driver, test_spec):
    """Run all conversion tests in the defined list of test specifications
    """
    run_test_conversion_with_gui(test_spec=test_spec,
                                 driver=driver)
