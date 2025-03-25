#!/usr/bin/env python

# Selenium test script for PSDI Data Conversion Service.

import os
import time

from pathlib import Path
import pytest
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

env_driver = os.environ.get("DRIVER")
origin = os.environ.get("ORIGIN")

if (env_driver is None):
    driver_path = GeckoDriverManager().install()
else:
    driver_path = env_driver

if (origin is None):
    print("ORIGIN environment variable must be set.")
    exit(1)

# Standard timeout at 10 seconds
TIMEOUT = 10


def wait_for_element(driver: WebDriver, xpath: str, by=By.XPATH):
    """Shortcut for boilerplate to wait until a web element is visible"""
    WebDriverWait(driver, TIMEOUT).until(EC.visibility_of_element_located((by, xpath)))


def wait_and_find_element(driver: WebDriver, xpath: str, by=By.XPATH) -> EC.WebElement:
    """Finds a web element, after first waiting to ensure it's visible"""
    wait_for_element(driver, xpath, by=by)
    return driver.find_element(by, xpath)


@pytest.fixture(scope="module")
def driver():
    """Get a headless Firefox web driver"""
    opts = FirefoxOptions()
    opts.add_argument("--headless")
    ff_driver = webdriver.Firefox(service=FirefoxService(driver_path),
                                  options=opts)
    yield ff_driver
    ff_driver.quit()


def test_initial_frontpage(driver: WebDriver):

    driver.get(f"{origin}/")

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

    input_file = Path.cwd().joinpath("files", f"{test_file}.cdxml")
    output_file = Path.home().joinpath("Downloads", f"{test_file}.inchi")
    log_file = Path.home().joinpath("Downloads", f"{test_file}.log.txt")

    # Remove test files from Downloads directory if they exist.

    if (Path.is_file(log_file)):
        Path.unlink(log_file)

    if (Path.is_file(output_file)):
        Path.unlink(output_file)

    driver.get(f"{origin}/")

    wait_for_element(driver, "//select[@id='fromList']/option")

    # Select cdxml from the 'from' list.
    driver.find_element(By.XPATH, "//select[@id='fromList']/option[contains(.,'cdxml: ChemDraw CDXML')]").click()

    # Select InChI from the 'to' list.
    driver.find_element(By.XPATH, "//select[@id='toList']/option[contains(.,'inchi: InChI')]").click()

    # Select Open Babel from the available conversion options list.
    driver.find_element(By.XPATH, "//select[@id='success']/option[contains(.,'Open Babel')]").click()

    # Click on the "Yes" button.
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
            assert False, f"Download of {output_file} and {log_file} timed out"

    time.sleep(1)

    # Verify that the InChI file is correct.

    assert output_file.read_text().strip() == "InChI=1S/C12NO/c1-12(2)6-7-13-11-5-4-9(14-3)8-10(11)12"
