"""@file psdi_data_conversion/dist.py

Created 2025-02-25 by Bryan Gillis.

Functions and utilities related to handling multiple user OSes and distributions
"""

import os
import psdi_data_conversion
import sys

# Labels for each platform (which we use for the folder in this project), and the head of the name each platform will
# have in `sys.platform`

LINUX_LABEL = "linux"
LINUX_NAME_HEAD = "linux"

WINDOWS_LABEL = "windows"
WINDOWS_NAME_HEAD = "win"

MAC_LABEL = "mac"
MAC_NAME_HEAD = "darwin"

D_DIST_NAME_HEADS = {LINUX_LABEL: LINUX_NAME_HEAD,
                     WINDOWS_LABEL: WINDOWS_NAME_HEAD,
                     MAC_LABEL: MAC_NAME_HEAD, }

# Determine the dist when this module is first imported
_dist: str | None = None
for label, name_head in D_DIST_NAME_HEADS.items():
    if sys.platform.startswith(name_head):
        _dist = label
        break
DIST = _dist

# Determine the fully-qualified binary directory when this module is first imported
BIN_DIR: str = os.path.join(psdi_data_conversion.__path__[0], "bin")


def get_bin_path(bin_name: str) -> str | None:
    """Gets the path to an appropriate binary for the user's OS/distribution, if one exists

    Parameters
    ----------
    bin_name : str
        The name of the binary relative to the ``psdi_data_conversion/bin/$DIST`` directory

    Returns
    -------
    str | None
        If an appropriate binary exists for the user's distribution, a fully-qualified path to it. Otherwise, None
    """

    # If DIST is None, then the user's OS/distribution is unsupported
    if not DIST:
        return None

    bin_path = os.path.join(BIN_DIR, DIST, bin_name)

    # Check if the binary exists in the path for the user's OS/distribution
    if not os.path.isfile(bin_path):
        return None

    return bin_path


def bin_exists(bin_name: str) -> bool:
    """Gets whether or not a binary of the given name exists for the user's distribution
    """

    return get_bin_path(bin_name) is not None
