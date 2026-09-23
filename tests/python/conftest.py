"""tests/python/conftest.py
===============

Created 2026-09-22 by Bryan Gillis.

Setup for all tests in this folder
"""

import pytest

# Compatibility patch for subtests in case a version of pytest prior to 9.0 is used, substituting a dummy fixture. It
# won't be quite as informative if tests fail, but will still fail properly if anything fails and succeed if everything
# succeeds
try:
    from pytest import Subtests
    _ = Subtests
except ImportError:

    from psdi_data_conversion.compatibility import DummySubtests

    @pytest.fixture
    def subtests():
        return DummySubtests()
