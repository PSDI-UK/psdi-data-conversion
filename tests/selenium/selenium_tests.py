#!/usr/bin/env python

# Selenium test script for PSDI Data Conversion Service.

import os
import time
import unittest

from pathlib import Path
from selenium.common.exceptions import NoSuchElementException
from selenium import webdriver
from selenium.webdriver.common.alert import Alert
from selenium.webdriver.common.by import By
from selenium.webdriver.firefox.service import Service as FirefoxService
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

# print(f"origin: {origin}")
# print(f"driver: {driver_path}")

driver = webdriver.Firefox(service=FirefoxService(driver_path))


class TestBasicOperations(unittest.TestCase):

    def test_initial_frontpage(self):

        driver.get(f"{origin}/")

        # Check that the front page contains the header "Data Conversion Service".

        WebDriverWait(driver, 10).until(EC.visibility_of_element_located((By.XPATH, "//header")))
        element = driver.find_element(By.XPATH, "//header//h5")
        self.assertEqual(element.text, "Data Conversion Service")

        # Check that the 'from' and 'to' lists contains "abinit" and "acesin" respectively.

        WebDriverWait(driver, 10).until(EC.visibility_of_element_located(
            (By.XPATH, "//select[@id='fromList']/option")))
        driver.find_element(By.XPATH, "//select[@id='fromList']/option[contains(.,'abinit: ABINIT output')]")

        WebDriverWait(driver, 10).until(EC.visibility_of_element_located((By.XPATH, "//select[@id='toList']/option")))
        driver.find_element(By.XPATH, "//select[@id='toList']/option[contains(.,'acesin: ACES input')]")

        # Check that the available conversions list is empty.

        with self.assertRaises(NoSuchElementException):
            driver.find_element(By.XPATH, "//select[@id='success']/option")

    def test_cdxml_to_inchi_conversion(self):

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

        WebDriverWait(driver, 10).until(EC.visibility_of_element_located(
            (By.XPATH, "//select[@id='fromList']/option")))

        # Select cdxml from the 'from' list.

        driver.find_element(By.XPATH, "//select[@id='fromList']/option[contains(.,'cdxml: ChemDraw CDXML')]").click()

        # Select InChI from the 'to' list.

        driver.find_element(By.XPATH, "//select[@id='toList']/option[contains(.,'inchi: InChI')]").click()

        # Select Open Babel from the available conversion options list.

        driver.find_element(By.XPATH, "//select[@id='success']/option[contains(.,'Open Babel')]").click()

        # Click on the "Yes" button.

        driver.find_element(By.XPATH, "//input[@id='yesButton']").click()

        # Select the input file.

        WebDriverWait(driver, 10).until(EC.visibility_of_element_located((By.XPATH, "//input[@id='fileToUpload']")))
        driver.find_element(By.XPATH, "//input[@id='fileToUpload']").send_keys(str(input_file))

        # Click on the "Convert" button.

        WebDriverWait(driver, 10).until(EC.visibility_of_element_located((By.XPATH, "//input[@id='uploadButton']")))
        driver.find_element(By.XPATH, "//input[@id='uploadButton']").click()

        # Handle alert box.

        WebDriverWait(driver, 10).until(EC.alert_is_present())
        Alert(driver).dismiss()

        # Wait until files exist.

        while (not Path.is_file(log_file)) or (not Path.is_file(output_file)):
            time.sleep(1)

        time.sleep(1)

        # Verify that the InChI file is correct.

        self.assertEqual(output_file.read_text().strip(), "InChI=1S/C12NO/c1-12(2)6-7-13-11-5-4-9(14-3)8-10(11)12")


if __name__ == '__main__':
    unittest.main()
