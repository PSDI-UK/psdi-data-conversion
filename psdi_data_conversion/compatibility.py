"""psdi_data_conversion/compatibility.py
===============

Created 2026-09-22 by Bryan Gillis.

Different versions of methods etc. to help with compatibility with different versions of Python. Official tools are
used where possible, with fallback to implementations here if not
"""

from collections.abc import Generator, Iterable

__all__ = ["batched"]

try:
    from itertools import batched
except ImportError:

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
