"""@file tests/dist_test.py

Created 2025-02-25 by Bryan Gillis.

Tests of the `dist` module, for determining the user's platform and appropriate binaries
"""

import os
import sys
from unittest.mock import patch
from psdi_data_conversion import dist


def test_get_dist():
    """Test that the dist is determined correctly for each platform
    """
    # Test each known platform
    for label, platform in dist.D_DIST_NAME_HEADS.items():
        with patch.object(sys, 'platform', platform):
            assert dist.get_dist() == label

    # Test an unknown platform
    with patch.object(sys, 'platform', "unknown"):
        assert dist.get_dist() is None


def test_get_bin():
    """Test that binaries can be found correctly for each platform
    """

    # Confirm the base binary directory is correct
    assert os.path.isdir(dist.BIN_DIR)

    # Test that binaries are found only when we expect them to be
    for bin_name, platform, should_exist in (("atomsk", dist.LINUX_NAME_HEAD, True),
                                             ("atomsk", dist.WINDOWS_NAME_HEAD, False),
                                             ("atomsk", dist.MAC_NAME_HEAD, True),
                                             ("c2x", dist.LINUX_NAME_HEAD, True),
                                             ("c2x", dist.WINDOWS_NAME_HEAD, False),
                                             ("c2x", dist.MAC_NAME_HEAD, True),
                                             ("bad_bin", dist.LINUX_NAME_HEAD, False),
                                             ("c2x", "unknown", False),):
        with patch.object(sys, 'platform', platform):
            bin_path = dist.get_bin_path(bin_name)
            if should_exist:
                assert os.path.isfile(bin_path), (bin_name, platform)
                assert dist.bin_exists(bin_name)
            else:
                assert bin_path is None or "psdi_data_conversion" not in bin_path, (bin_name, platform)
                if bin_path is None:
                    assert not dist.bin_exists(bin_name)
                else:
                    assert dist.bin_exists(bin_name)
