"""@file tests/security_test.py

Created 2025-02-24 by Bryan Gillis.

Tests of security-related functions
"""

from psdi_data_conversion.security import char_is_safe, string_is_safe


def test_char_is_safe():
    """Tests of the `char_is_safe` function.
    """
    # Test some strings we expect to pass
    assert char_is_safe("a")
    assert char_is_safe("0")
    assert char_is_safe("W")
    assert char_is_safe("0")
    assert char_is_safe(".")
    assert char_is_safe("-")
    assert char_is_safe("+")
    assert char_is_safe("*")
    assert char_is_safe("=")
    assert char_is_safe("$")
    assert char_is_safe(" ")
    assert char_is_safe("\t")
    assert char_is_safe("\n")
    assert char_is_safe("/")
    assert char_is_safe("\\")

    # Test some strings we expect to fail
    assert not char_is_safe("}")
    assert not char_is_safe("\"")
    assert not char_is_safe("'")
    assert not char_is_safe(";")
    assert not char_is_safe("&")
    assert not char_is_safe("|")
    assert not char_is_safe("`")
    assert not char_is_safe("")
    assert not char_is_safe("aa")


def test_string_is_safe():
    """Tests of the `string_is_safe` function.
    """
    # Test some strings we expect to pass
    assert string_is_safe("alphabet_and_D16IT5")
    assert string_is_safe("like_a_file.name")
    assert string_is_safe("-0.5+15*2/3=0")
    assert string_is_safe("\tmy\nmultiline string")
    assert string_is_safe("path/to/my/file")
    assert string_is_safe("C:\\windows\\style\\path")

    # Test some strings we expect to fail
    assert not string_is_safe("file.txt; hack_the_mainframe")
    assert not string_is_safe("} ")
    assert not string_is_safe("this & that")
    assert not string_is_safe("now | never")
    assert not string_is_safe("`backticked command`")
