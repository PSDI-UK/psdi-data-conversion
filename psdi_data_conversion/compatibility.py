"""psdi_data_conversion/compatibility.py
===============

Created 2026-09-22 by Bryan Gillis.

Different versions of methods etc. to help with compatibility with different versions of Python. Official tools are
used where possible, with fallback to implementations here if not
"""

from contextlib import AbstractContextManager, contextmanager

__all__ = ["batched", "DummySubtests"]

try:
    from itertools import batched
except ImportError:

    from collections.abc import Generator, Iterable

    def batched(x: Iterable, n: int) -> Generator[tuple]:
        """Fallback implementation of `itertools.batched`, which was introduced in Python 3.12"""
        ix = iter(x)
        while True:
            # We can't quite use list comprehension here to fully match the behaviour of `itertools.batched``, as that
            # will return a partial batch at the end if the iterable isn't evenly divisible by `n`, while list
            # comprehension won't
            l_out = []
            for _ in range(n):
                try:
                    l_out.append(next(ix))
                except StopIteration:
                    pass
            if not l_out:
                return
            yield tuple(l_out)


@contextmanager
def _dummy_context_manager():
    yield
    return


class DummySubtests:
    """Fallback for Pytest's `subtests` fixtrue, which just sets up a dummy context manager which will result
    in any subtest failing immediately failing the enclosing test"""

    def __init__(*args, **kwargs):
        return

    def test(self, *args, **kwargs) -> AbstractContextManager:
        return _dummy_context_manager()
