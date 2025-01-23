"""@file psdi-data-conversion/psdi_data_conversion/converter.py

Created 2024-12-10 by Bryan Gillis.

Class and functions to perform file conversion
"""

from psdi_data_conversion.converters.base import FileConverter


def get_converter(**converter_kwargs) -> FileConverter:
    """_summary_
    """
    return FileConverter(**converter_kwargs)


def run_converter(**converter_kwargs):
    """Shortcut to create and run a FileConverter in one step
    """

    return get_converter(**converter_kwargs).run()
