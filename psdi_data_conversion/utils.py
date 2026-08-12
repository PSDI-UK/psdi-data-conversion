"""
# utils.py

Miscellaneous utility functions used by this project
"""


import json
import sys
import textwrap
from functools import lru_cache
from importlib.metadata import Distribution
from pathlib import Path

from psdi_data_conversion.constants import TERM_WIDTH
from psdi_data_conversion.file_io import get_package_path

# Common types
JsonDict = dict[str, None | int | str | bool | dict | list]
JsonMainDict = dict[str, None | int | str | bool | JsonDict | list[JsonDict]]


class TextColors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


def get_wrapped_str(s: str, **kwargs):
    """Get a string wrapped to the terminal width"""
    return textwrap.fill(s, width=TERM_WIDTH, **kwargs)


def print_wrap(s: str, newline=False, err=False, **kwargs):
    """Print a string wrapped to the terminal width
    """
    if err:
        file = sys.stderr
    else:
        file = sys.stdout
    for line in s.split("\n"):
        print(get_wrapped_str(line, **kwargs), file=file)
    if newline:
        print("")


def regularize_name(name: str):
    """Regularizes a name for comparisons, making it lowercase and stripping spaces

    Parameters
    ----------
    name : str
        The name, e.g. "Open Babel"

    Returns
    -------
    str
        The regularized name, e.g. "openbabel"
    """
    return name.lower().replace(" ", "")


def in_editable_mode():
    """Checks if the `psdi_data_conversion` module is installed in editable mode
    """

    direct_url = Distribution.from_name("psdi_data_conversion").read_text("direct_url.json")
    is_editable: bool = json.loads(direct_url).get("dir_info", {}).get("editable", False)

    return is_editable


def confirm_editable_mode():
    """Checks if the `psdi_data_conversion` module is installed in editable mode, and exits the program if not
    """
    if not in_editable_mode():
        print_wrap(f"{TextColors.FAIL}ERROR:{TextColors.ENDC} To run this script, the package must be installed in "
                   f"editable mode. Please reinstall with:\n")
        print(f"{TextColors.WARNING}pip install --editable .{TextColors.ENDC}\n")
        print_wrap("and re-run this script.")
        exit(1)


@lru_cache(maxsize=1)
def get_project_path() -> Path:
    """Gets the absolute path to where the project is on disk, using the package path to find it and checking that it
    contains the expected files
    """

    project_path = (get_package_path() / "..").resolve()

    # Check that the project path contains the expected test_data folder
    if not (project_path / "pyproject.toml").is_file():
        raise FileNotFoundError(f"Project path was expected to be '{project_path}', but this does not contain the "
                                f"expected file 'pyproject.toml'")

    return project_path
