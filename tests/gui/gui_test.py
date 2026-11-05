#!/usr/bin/env python

# Selenium test script for PSDI Data Conversion Service.

import os
from multiprocessing import Process

import pytest
import requests
import time
import json

from http import HTTPStatus

from psdi_data_conversion.testing.constants import DEFAULT_ORIGIN
from psdi_data_conversion.testing.conversion_test_specs import l_gui_test_specs
from psdi_data_conversion.testing.gui import GuiTestSpecRunner, wait_for_cover_hidden, wait_for_element

# Skip all tests in this module if required packages for GUI testing aren't installed
try:
    from selenium import webdriver
    from selenium.common.exceptions import NoSuchElementException
    from selenium.webdriver import FirefoxOptions
    from selenium.webdriver.common.by import By
    from selenium.webdriver.firefox.service import Service as FirefoxService
    from selenium.webdriver.firefox.webdriver import WebDriver
    from webdriver_manager.firefox import GeckoDriverManager

    from psdi_data_conversion.gui.setup import start_app

except ImportError:
    # We put the importorskip commands here rather than above so that standard imports can be used by static analysis
    # tools where possible, and the importorskip is used here so pytest will stop processing immediately if things can't
    # be imported - pytest.mark.skip won't do that
    pytest.importorskip("Flask")
    pytest.importorskip("selenium")
    pytest.importorskip("webdriver_manager.firefox")


import psdi_data_conversion

origin = os.environ.get("ORIGIN", DEFAULT_ORIGIN)
api_origin = f"{origin}/api/v1"

@pytest.fixture(scope="module", autouse=True)
def common_setup():
    """Autouse fixture which starts the app before tests and stops it afterwards"""

    # If the origin is set to something else, don't start the local server here
    if origin != DEFAULT_ORIGIN:
        yield
        return

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

    element = wait_for_element(driver, "//header//h5")
    assert element.text == "Data Conversion Service"

    # Check that the 'from' and 'to' lists contains "abinit" and "acesin" respectively.

    wait_for_element(driver, "//select[@id='fromList']/option")
    driver.find_element(By.XPATH, "//select[@id='fromList']/option[contains(.,'abinit: ABINIT output')]")

    wait_for_element(driver, "//select[@id='toList']/option")
    driver.find_element(By.XPATH, "//select[@id='toList']/option[contains(.,'acesin: ACES input')]")

    # Check that the available conversions list is empty.

    with pytest.raises(NoSuchElementException):
        driver.find_element(By.XPATH, "//select[@id='success']/option")


@pytest.fixture(scope="module")
def gui_test_spect_runner(driver):
    return GuiTestSpecRunner(driver=driver, origin=origin)


@pytest.mark.parametrize("test_spec", l_gui_test_specs,
                         ids=lambda x: x.name)
def test_conversions(gui_test_spect_runner, test_spec, subtests):
    """Run all conversion tests in the defined list of test specifications
    """
    gui_test_spect_runner.run(test_spec, subtests)


# API tests

openbabel = "f2d7f6c0-04c3-4677-b42e-1e7181d3618d"
cdxml = "c868b7b8-6d84-40de-899f-cf6071c50982"
inchi = "d2763887-718b-4819-a321-54b821379f07"
unknown = "00000000-0000-4000-8000-000000000000"

def wait_for_server(url: str, timeout: float = 10.0) -> None:
    deadline = time.monotonic() + timeout

    while time.monotonic() < deadline:
        try:
            response = requests.get(url, timeout=2)
            if response.ok:
                return
        except requests.RequestException:
            pass  # The server may not be listening yet.

        time.sleep(0.2)

    raise TimeoutError(f"Server did not become ready at {url} within {timeout}s")


def test_api_get_all_converters():

    wait_for_server(f"{origin}/")

    response = requests.get(f"{api_origin}/converters", headers={ "Accept": "application/json" })

    assert response.status_code == 200
    assert "application/json" in response.headers["Content-Type"]

    data = response.json()

    # The result should be a list.
    assert isinstance(data, list)

    # The result should not be an empty list.
    assert data

    # Each result should have a name entry.
    assert all("name" in item for item in data)

    # Each result should have a id entry.
    assert all("id" in item for item in data)

    # Each result should have a description entry.
    assert all("description" in item for item in data)

    # Each result should have a url entry.
    assert all("url" in item for item in data)


def test_api_get_converter():

    wait_for_server(f"{origin}/")

    response = requests.get(f"{api_origin}/converters/{openbabel}",
        headers={ "Accept": "application/json" })

    assert response.status_code == 200
    assert "application/json" in response.headers["Content-Type"]

    data = response.json()

    # The result should be a dict.
    assert isinstance(data, dict)

    # The result should not be empty list.
    assert data

    assert data["id"] == openbabel
    assert data["name"] == "openbabel"
    assert data["description"] == "Comprehensive converter"
    assert data["url"] == "https://open-babel.readthedocs.io/en/latest/"


def test_api_get_unknown_converter():

    wait_for_server(f"{origin}/")

    response = requests.get(f"{api_origin}/converters/{unknown}",
        headers={ "Accept": "application/json" })

    assert response.status_code == 404

    data = response.json()

    assert data["error"] == f"Converter {unknown} was not found"


def test_api_get_converter_formats():

    wait_for_server(f"{origin}/")

    response = requests.get(f"{api_origin}/converters/{openbabel}/formats",
        headers={ "Accept": "application/json" })

    assert response.status_code == 200
    assert "application/json" in response.headers["Content-Type"]

    data = response.json()

    # The result should be a dict with inputFormats and outputFormats lists.
    assert isinstance(data, dict)
    assert isinstance(data["inputFormats"], list)
    assert isinstance(data["outputFormats"], list)

    # The result should not be empty list.
    assert data

    # Each input format should have an id.
    assert all("id" in item for item in data["inputFormats"])

    # Each input format should have a label.
    assert all("label" in item for item in data["inputFormats"])

    # Each output format should have an id.
    assert all("id" in item for item in data["outputFormats"])

    # Each output format should have a label.
    assert all("label" in item for item in data["outputFormats"])


def test_api_get_unknown_converter_formats():

    wait_for_server(f"{origin}/")

    response = requests.get(f"{api_origin}/converters/{unknown}/formats",
        headers={ "Accept": "application/json" })

    assert response.status_code == 404

    data = response.json()

    assert data["error"] == f"Converter {unknown} was not found"


def test_api_get_conversion_info():

    wait_for_server(f"{origin}/")

    response = requests.get(f"{api_origin}/converters/{openbabel}/conversions/{cdxml}/{inchi}",
        headers={ "Accept": "application/json" })

    assert response.status_code == 200
    assert "application/json" in response.headers["Content-Type"]

    data = response.json()

    # The result should be a dict with inputFormats and outputFormats lists.
    assert isinstance(data, dict)
    assert isinstance(data["quality"], str)
    assert isinstance(data["from"], dict)
    assert isinstance(data["to"], dict)


    assert data["quality"] == "very good"

    assert data["from"]["compositionSupported"] == True
    assert data["from"]["connectionsSupported"] == True
    assert data["from"]["twoDimSupported"] == True
    assert data["from"]["threeDimSupported"] == True

    assert data["to"]["compositionSupported"] == True
    assert data["to"]["connectionsSupported"] == True
    assert data["to"]["twoDimSupported"] == False
    assert data["to"]["threeDimSupported"] == False


def test_api_get_conversion_info_unknown_converter():

    wait_for_server(f"{origin}/")

    response = requests.get(f"{api_origin}/converters/{unknown}/conversions/{cdxml}/{inchi}",
        headers={ "Accept": "application/json" })

    assert response.status_code == 404

    data = response.json()

    assert data["error"] == f"Converter {unknown} was not found"


def test_api_get_conversion_info_unknown_from_format():

    wait_for_server(f"{origin}/")

    response = requests.get(f"{api_origin}/converters/{openbabel}/conversions/{unknown}/{inchi}",
        headers={ "Accept": "application/json" })

    assert response.status_code == 404

    data = response.json()

    assert data["error"] == f"Format {unknown} was not found"


def test_api_get_conversion_info_unknown_to_format():

    wait_for_server(f"{origin}/")

    response = requests.get(f"{api_origin}/converters/{openbabel}/conversions/{cdxml}/{unknown}",
        headers={ "Accept": "application/json" })

    assert response.status_code == 404

    data = response.json()

    assert data["error"] == f"Format {unknown} was not found"


def test_api_get_cdxml_format():

    wait_for_server(f"{origin}/")

    response = requests.get(f"{api_origin}/formats/{cdxml}",
        headers={ "Accept": "application/json" })

    assert response.status_code == 200
    assert "application/json" in response.headers["Content-Type"]

    data = response.json()

    # The result should be a dict.
    assert isinstance(data, dict)

    assert data["id"] == cdxml
    assert data["name"] == "cdxml"
    assert data["twoDim"] == True
    assert data["threeDim"] == True


def test_api_get_unknown_format():

    wait_for_server(f"{origin}/")

    response = requests.get(f"{api_origin}/formats/{unknown}",
        headers={ "Accept": "application/json" })

    assert response.status_code == 404

    data = response.json()

    assert data["error"] == f"Format {unknown} was not found"


def test_api_convert_success():

    wait_for_server(f"{origin}/")

    with open("test_data/standard_test.cdxml", "rb") as input_file:
        response = requests.post(
            f"{api_origin}/converters/{openbabel}/conversions/{cdxml}/{inchi}/convert",
            headers={
                "Accept": "application/octet-stream",
            },
            data={
            },
            files={
                "args": (
                    None,
                    json.dumps({}),
                    "application/json",
                ),
                "input_file": (
                    "standard_test.cdxml",
                    input_file,
                    "application/octet-stream",
                ),
            },
        )

    assert response.status_code == HTTPStatus.OK
    assert response.headers["Content-Type"].startswith(
        "application/octet-stream"
    )

    assert response.content == (
        b"InChI=1S/C12H19NO/c1-12(2)6-7-13-11-5-4-9(14-3)8-10(11)12/"
        b"h6-7,9,13H,4-5,8H2,1-3H3\n"
    )


def test_api_convert_success_output_flag():

    wait_for_server(f"{origin}/")

    with open("test_data/standard_test.cdxml", "rb") as input_file:
        response = requests.post(
            f"{api_origin}/converters/{openbabel}/conversions/{cdxml}/{inchi}/convert",
            headers={
                "Accept": "application/octet-stream",
            },
            data={
            },
            files={
                "args": (
                    None,
                    json.dumps({"t": True}),
                    "application/json",
                ),
                "input_file": (
                    "standard_test.cdxml",
                    input_file,
                    "application/octet-stream",
                ),
            },
        )

    assert response.status_code == HTTPStatus.OK
    assert response.headers["Content-Type"].startswith(
        "application/octet-stream"
    )

    assert response.content == (
        b"InChI=1S/C12H19NO/c1-12(2)6-7-13-11-5-4-9(14-3)8-10(11)12/"
        b"h6-7,9,13H,4-5,8H2,1-3H3 5\n"
    )


def test_api_convert_success_output_option():

    wait_for_server(f"{origin}/")

    with open("test_data/standard_test.cdxml", "rb") as input_file:
        response = requests.post(
            f"{api_origin}/converters/{openbabel}/conversions/{cdxml}/{inchi}/convert",
            headers={
                "Accept": "application/octet-stream",
            },
            data={
            },
            files={
                "args": (
                    None,
                    json.dumps({"T": "/formula"}),
                    "application/json",
                ),
                "input_file": (
                    "standard_test.cdxml",
                    input_file,
                    "application/octet-stream",
                ),
            },
        )

    assert response.status_code == HTTPStatus.OK
    assert response.headers["Content-Type"].startswith(
        "application/octet-stream"
    )

    assert response.content == (
        b"InChI=1S/C12H19NO\n"
    )


def test_api_convert_unknown_converter():

    wait_for_server(f"{origin}/")

    with open("test_data/standard_test.cdxml", "rb") as input_file:
        response = requests.post(
            f"{api_origin}/converters/{unknown}/conversions/{cdxml}/{inchi}/convert",
            headers={
                "Accept": "application/octet-stream",
            },
            data={
            },
            files={
                "args": (
                    None,
                    json.dumps({"t": True}),
                    "application/json",
                ),
                "input_file": (
                    "standard_test.cdxml",
                    input_file,
                    "application/octet-stream",
                ),
            },
        )

    assert response.status_code == HTTPStatus.NOT_FOUND

    data = response.json()

    assert data["error"] == f"Converter {unknown} was not found"


def test_api_convert_unknown_from_format():

    wait_for_server(f"{origin}/")

    with open("test_data/standard_test.cdxml", "rb") as input_file:
        response = requests.post(
            f"{api_origin}/converters/{openbabel}/conversions/{unknown}/{inchi}/convert",
            headers={
                "Accept": "application/octet-stream",
            },
            data={
            },
            files={
                "args": (
                    None,
                    json.dumps({"t": True}),
                    "application/json",
                ),
                "input_file": (
                    "standard_test.cdxml",
                    input_file,
                    "application/octet-stream",
                ),
            },
        )

    assert response.status_code == HTTPStatus.NOT_FOUND

    data = response.json()

    assert data["error"] == f"Format {unknown} was not found"


def test_api_convert_unknown_to_format():

    wait_for_server(f"{origin}/")

    with open("test_data/standard_test.cdxml", "rb") as input_file:
        response = requests.post(
            f"{api_origin}/converters/{openbabel}/conversions/{cdxml}/{unknown}/convert",
            headers={
                "Accept": "application/octet-stream",
            },
            data={
            },
            files={
                "args": (
                    None,
                    json.dumps({"t": True}),
                    "application/json",
                ),
                "input_file": (
                    "standard_test.cdxml",
                    input_file,
                    "application/octet-stream",
                ),
            },
        )

    assert response.status_code == HTTPStatus.NOT_FOUND

    data = response.json()

    assert data["error"] == f"Format {unknown} was not found"
