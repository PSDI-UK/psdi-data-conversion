#!/usr/bin/env python

# Selenium test script for PSDI Data Conversion Service.

import os
from multiprocessing import Process

import pytest

from psdi_data_conversion.testing.constants import DEFAULT_ORIGIN
from psdi_data_conversion.testing.gui import (GuiTestSpecRunner, wait_and_find_element, wait_for_cover_hidden,
                                              wait_for_element)
from psdi_data_conversion.testing.conversion_test_specs import l_gui_test_specs

# Skip all tests in this module if required packages for GUI testing aren't installed
try:
    from selenium.common.exceptions import NoSuchElementException
    from selenium import webdriver
    from selenium.webdriver.common.by import By
    from selenium.webdriver import FirefoxOptions
    from selenium.webdriver.firefox.service import Service as FirefoxService
    from selenium.webdriver.firefox.webdriver import WebDriver
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

origin = os.environ.get("ORIGIN", DEFAULT_ORIGIN)


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


def test_initial_frontpage(driver: WebDriver):

    # Load the home page and wait for the page cover to be removed
    driver.get(f"{origin}/")
    wait_for_cover_hidden(driver)

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


@pytest.mark.parametrize("test_spec", l_gui_test_specs,
                         ids=lambda x: x.name)
def test_conversions(driver, test_spec):
    """Run all conversion tests in the defined list of test specifications
    """
    GuiTestSpecRunner(driver=driver, origin=origin).run(test_spec)
