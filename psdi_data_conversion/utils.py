"""
# utils.py

Miscellaneous utility functions used by this project
"""


import json
import sys
import textwrap
from importlib.metadata import Distribution

from psdi_data_conversion.constants import TERM_WIDTH


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
