"""@file psdi_data_conversion/security.py

Created 2025-02-24 by Bryan Gillis.

Functions related to security
"""

import re

# Short list of safe allowed characters:
# \w: All letters and digits
# \s: All whitespace characters
# .: Period
# \-: Hyphen
# :: Colon
# +: Plus symbol
# *: Asterisk
# =: Equals sign
# $: Dollar sign
# /: Forward-slash
# \\: Backslash
_SAFE_CHARS = r"[\w\s.\-:+*=$/\\]"
SAFE_CHAR_RE = re.compile(_SAFE_CHARS)
SAFE_STRING_RE = re.compile(f"{_SAFE_CHARS}*")


def char_is_safe(s: str):
    """Checks whether a character is in the set of predefined safe characters.

    Will return False if the string does not contain exactly one character.
    """
    return bool(SAFE_CHAR_RE.fullmatch(s))


def string_is_safe(s: str):
    """Checks whether a string contains only characters in the predefined list of safe characters.
    """
    return bool(SAFE_STRING_RE.fullmatch(s))
