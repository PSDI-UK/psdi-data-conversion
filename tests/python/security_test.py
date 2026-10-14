"""@file tests/security_test.py

Created 2025-02-24 by Bryan Gillis.

Tests of security-related functions
"""

import os
import shlex
import sys
from unittest.mock import patch

import pytest

from psdi_data_conversion.constants import CONVERTER_OB
from psdi_data_conversion.main import main
from psdi_data_conversion.security import char_is_safe, string_is_safe
from psdi_data_conversion.testing.utils import get_input_test_data_loc


def run_with_arg_string(s):
    """Runs the convert script with the provided argument string
    """
    l_args = shlex.split("test " + s)
    with patch.object(sys, 'argv', l_args):
        main()


@pytest.mark.parametrize("c", ("a", "0", "W", "0", ".", "-", "+", "*", "=", "$", " ", "\t", "\n", "/", "\\"))
def test_char_is_safe_true(c):
    """Test the `char_is_safe` function returns True for various safe characters"""
    assert char_is_safe(c)


@pytest.mark.parametrize("c", ("}", "\"", "'", ";", "&", "|", "`", "", "aa"))
def test_char_is_safe_false(c):
    """Test the `char_is_safe` function returns False for various unsafe characters"""
    assert not char_is_safe(c)


@pytest.mark.parametrize("s", ("alphabet_and_D16IT5", "like_a_file.name", "-0.5+15*2/3=0", "\tmy\nmultiline string",
                               "path/to/my/file", "C:\\windows\\style\\path"))
def test_string_is_safe_true(s):
    """Test the `string_is_safe` function returns True for various safe strings"""
    assert string_is_safe(s)


@pytest.mark.parametrize("s", ("file.txt; hack_the_mainframe", "} ", "this & that", "now | never",
                               "`backticked command`"))
def test_string_is_safe_false(s):
    """Test the `string_is_safe` function returns False for various unsafe strings"""
    assert not string_is_safe(s)


@pytest.mark.parametrize("from_flags, to_flags, from_options, to_options", (("'; stop'", None, None, None),
                                                                            (None, "'`command`'", None, None),
                                                                            (None, None, "'this & that'", None),
                                                                            (None, None, None, "(nope)"),
                                                                            (None, None, None, "'now || never'"),
                                                                            (None, None, None, "'} a b'"),))
def test_format_arg_security(from_flags, to_flags, from_options, to_options, tmp_path_factory, capsys, subtests):
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

    arg_string = f"{input_filename} -t {to_format} -i {input_dir} -o {output_dir} -w {CONVERTER_OB}"

    if from_flags is not None:
        arg_string += f" --from-flags {from_flags}"
    if to_flags is not None:
        arg_string += f" --to-flags {to_flags}"
    if from_options is not None:
        arg_string += f" --from-options {from_options}"
    if to_options is not None:
        arg_string += f" --to-options {to_options}"

    # Run the conversion
    with subtests.test("system exit raised"), pytest.raises(SystemExit):
        run_with_arg_string(arg_string)

    # Make sure the expected error is output, and not the full traceback
    captured = capsys.readouterr()
    with subtests.test("security error output to stderr"):
        assert "ERROR" in captured.err
        assert "security" in captured.err
    with subtests.test("no traceback in stdout"):
        assert "Traceback" not in captured.out
    with subtests.test("no traceback in stderr"):
        assert "Traceback" not in captured.err
