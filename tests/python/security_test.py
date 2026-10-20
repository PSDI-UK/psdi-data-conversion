"""@file tests/security_test.py

Created 2025-02-24 by Bryan Gillis.

Tests of security-related functions
"""

import os
import shlex
import sys
from unittest.mock import patch

import pytest
from psdi_data_conversion.main import main
from psdi_data_conversion.security import char_is_safe, string_is_safe
from psdi_data_conversion.testing.utils import get_input_test_data_loc


def run_with_arg_string(s):
    """Runs the convert script with the provided argument string
    """
    l_args = shlex.split("test " + s)
    with patch.object(sys, 'argv', l_args):
        main()


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


def test_format_arg_security(tmp_path_factory, capsys):
    """Test that format flags are processed correctly and results in the right conversions
    """
    input_dir = tmp_path_factory.mktemp("input")
    output_dir = tmp_path_factory.mktemp("output")

    test_filename_base = "caffeine"

    from_format = "inchi"
    to_format = "smi"

    input_filename = f"{test_filename_base}.{from_format}"

    # Symlink the input file from the test_data directory to the input directory
    os.symlink(os.path.join(get_input_test_data_loc(), input_filename),
               os.path.join(input_dir, input_filename))

    basic_arg_string = f"{input_filename} -t {to_format} -i {input_dir} -o {output_dir}"

    # Run for each set of flags that might trigger a security error
    for from_flags, to_flags, from_options, to_options in (("'; stop'", None, None, None),
                                                           (None, "'`command`'", None, None),
                                                           (None, None, "'this & that'", None),
                                                           (None, None, None, "(nope)"),
                                                           (None, None, None, "'now || never'"),
                                                           (None, None, None, "'} a b'"),):

        arg_string = basic_arg_string
        if from_flags is not None:
            arg_string += f" --from-flags {from_flags}"
        if to_flags is not None:
            arg_string += f" --to-flags {to_flags}"
        if from_options is not None:
            arg_string += f" --from-options {from_options}"
        if to_options is not None:
            arg_string += f" --to-options {to_options}"

        # Run the conversion
        with pytest.raises(SystemExit):
            run_with_arg_string(arg_string)

        # Make sure the expected error is output, and not the full traceback
        captured = capsys.readouterr()
        assert "ERROR" in captured.err, arg_string
        assert "security" in captured.err, arg_string
        assert "Traceback" not in captured.out, arg_string
        assert "Traceback" not in captured.err, arg_string
