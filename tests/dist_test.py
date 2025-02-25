"""@file tests/dist_test.py

Created 2025-02-25 by Bryan Gillis.

Tests of the `dist` module, for determining the user's platform and appropriate binaries
"""

import sys
from unittest.mock import patch
from psdi_data_conversion import dist


def test_get_dist():
    """Test that the dist is determined correctly for each platform
    """
    # Test each known platform
    for label, platform in dist.D_DIST_NAME_HEADS.items():
        with patch.object(sys, 'platform', platform):
            assert dist._get_dist() == label

    # Test an unknown platform
    with patch.object(sys, 'platform', "unknown"):
        assert dist._get_dist() is None
