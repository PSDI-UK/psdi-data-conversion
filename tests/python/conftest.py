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

# VSCode workaround so that if any subtests fail, the whole test is reported as failed in VSCode, from
# https://github.com/microsoft/vscode-python/issues/25824#issuecomment-4460161647
from collections.abc import Generator

_calling = set[str]()


@pytest.hookimpl(wrapper=True)
def pytest_runtest_call(item: pytest.Item) -> Generator[None]:
    _calling.add(item.nodeid)
    yield
    _calling.remove(item.nodeid)


@pytest.hookimpl(wrapper=True)
def pytest_runtest_makereport(item: pytest.Item,
                              call: pytest.CallInfo) -> Generator[None, pytest.TestReport, pytest.TestReport]:
    report = yield
    if call.when == "call" and item.nodeid in _calling:
        report.nodeid += "(subtest)"
    return report
