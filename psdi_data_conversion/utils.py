"""
# utils.py

Miscellaneous utility functions used by this project
"""


import json
import re
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
    """ANSI escape codes that can be used to color text printed to the terminal. E.g. to give text the header color,
    you could do `print(f"{TextColors.MAGENTA}Header text{TextColors.OFF}")`
    """

    # Text color codes

    RED = "\033[91m"
    """Start coloring red"""

    DARKRED = ERROR = FAIL = "\033[31m"
    """Start coloring dark red"""

    GREEN = SUCCESS = "\033[92m"
    """Start coloring green"""

    DARKGREEN = "\033[32m"
    """Start coloring dark green"""

    YELLOW = CODE = "\033[93m"
    """Start coloring yellow"""

    DARKYELLOW = WARNING = "\033[33m"
    """Start coloring dark yellow"""

    BLUE = ID = "\033[94m"
    """Start coloring blue"""

    DARKBLUE = "\033[34m"
    """Start coloring dark blue"""

    MAGENTA = "\033[95m"
    """Start coloring magenta"""

    DARKMAGENTA = "\033[35m"
    """Start coloring dark magenta"""

    CYAN = PATH = "\033[96m"
    """Start coloring cyan"""

    DARKCYAN = MESSAGE = "\033[36m"
    """Start coloring dark cyan"""

    # Text formatting codes

    BOLD = "\033[1m"
    """Start formatting bold - NOT compatible with coloring"""

    DIM = "\033[2m"
    """Start formatting dim (opposite of bold) - NOT compatible with coloring"""

    UNDERLINE = "\033[4m"
    """Start underlining - compatible with coloring"""

    # Combined codes

    HEADER = "\033[95m\033[4m"
    """Start header section - magenta underlined"""

    # Other codes

    OFF = "\033[0m"
    """End all coloring and formatting"""

    @classmethod
    def _get_codes(cls):
        return [x for x in dir(cls) if not x.startswith("_") and x.upper() == x]

    def display(self):
        """Displays all color codes"""
        l_codes_and_vals = [(x, getattr(self, x)) for x in self._get_codes()]
        l_codes_and_vals.sort(key=lambda x: int(x[1].replace("\033[", "").replace("m", "")) if x[1] else x[0])
        for code, val in l_codes_and_vals:
            print(f"{val}{code}{self.OFF}")

    def enable(self):
        """Enable colors"""
        l_codes = self._get_codes()
        for code in l_codes:
            setattr(self, code, getattr(type(self), code))

    def disable(self):
        """Disable colors"""
        l_codes = self._get_codes()
        for code in l_codes:
            setattr(self, code, "")


tc = TextColors()


def disable_colors():
    """Globally disable color formatting in output text"""
    tc.disable()


def enable_colors():
    """Globally (re)enable color formatting in output text"""
    tc.enable()


CONTROL_CODE_RE = re.compile("\033\\[\\d+?m")


def strip_control_codes(s: str):
    """Strip all control codes from a string"""
    return CONTROL_CODE_RE.sub("", s)


def displaylen(s: str):
    """Get the length of a string as it would be displayed in the terminal - this is, stripping out control codes"""
    return len(strip_control_codes(s))


def get_wrapped_str(s: str, color: str | None = None, **kwargs):
    """Get a string wrapped to the terminal width"""
    if color:
        s_colored = color+s+TextColors.OFF
    else:
        s_colored = s
    return textwrap.fill(s_colored, width=TERM_WIDTH, **kwargs)


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
        print_wrap(f"{TextColors.RED}ERROR:{TextColors.OFF} To run this script, the package must be installed in "
                   f"editable mode. Please reinstall with:\n")
        print_wrap("pip install --editable .\n", color=TextColors.WARNING)
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
