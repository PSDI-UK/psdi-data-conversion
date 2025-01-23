"""@file psdi-data-conversion/psdi_data_conversion/converter.py

Created 2024-12-10 by Bryan Gillis.

Class and functions to perform file conversion
"""

from psdi_data_conversion import constants as const
from psdi_data_conversion.converters.base import FileConverter, FileConverterInputException


def get_converter(converter=const.CONVERTER_DEFAULT, **converter_kwargs) -> FileConverter:
    """Get a FileConverter of the proper subclass for the requested converter type

    Parameters
    ----------
    converter : str
        The desired converter type, by default 'Open Babel'

    Returns
    -------
    FileConverter
        A subclassed FileConverter for the desired converter type

    Raises
    ------
    FileConverterInputException
        If the converter isn't recognized
    """
    if converter not in const.L_ALLOWED_CONVERTERS:
        raise FileConverterInputException(f"Converter {converter} not recognized. Allowed converters are: " +
                                          f"{const.L_ALLOWED_CONVERTERS}")
    return FileConverter(converter=converter, **converter_kwargs)


def run_converter(**converter_kwargs):
    """Shortcut to create and run a FileConverter in one step
    """

    return get_converter(**converter_kwargs).run()
