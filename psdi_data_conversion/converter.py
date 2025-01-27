"""@file psdi-data-conversion/psdi_data_conversion/converter.py

Created 2024-12-10 by Bryan Gillis.

Class and functions to perform file conversion
"""

import os
import importlib
import sys
import traceback
from typing import NamedTuple
from psdi_data_conversion import constants as const
from psdi_data_conversion.converters import base

import glob

# Find all modules for specific converters
l_converter_modules = glob.glob(os.path.dirname(base.__file__) + "/*.py")

try:

    class NameAndClass(NamedTuple):
        name: str
        converter_class: type[base.FileConverter]

    def get_converter_name_and_class(module_path: str) -> NameAndClass | None:

        module_name = os.path.splitext(os.path.basename(module_path))[0]

        # Skip the base module and the package __init__
        if module_name in ("base", "__init__"):
            return None

        package_name = "psdi_data_conversion.converters"
        module = importlib.import_module(f".{module_name}", package=package_name)

        # Check that the module defines a converter
        if not hasattr(module, "converter") or not issubclass(module.converter, base.FileConverter):
            print(f"ERROR: Module `{module_name}` in package `{package_name}` fails to define a converter to the "
                  "variable `converter` which is a subclass of `FileConverter`.", file=sys.stderr)
            return None

        converter_class = module.converter
        name = converter_class.name

        return NameAndClass(name, converter_class)

    # Get a list of all converter names and FileConverter subclasses
    l_converter_names_and_classes = [get_converter_name_and_class(module_name) for
                                     module_name in l_converter_modules]
    # Remove the None entry from the list, which corresponds to the 'base' module
    l_converter_names_and_classes = [x for x in l_converter_names_and_classes if x is not None]

    # Make constant dict and list of registered converters
    D_REGISTERED_CONVERTERS = dict(l_converter_names_and_classes)
    L_REGISTERED_CONVERTERS = [x for x in D_REGISTERED_CONVERTERS.keys()]

except Exception:
    print(f"ERROR: Failed to register converters. Exception was: \n{traceback.format_exc()}", file=sys.stderr)
    D_REGISTERED_CONVERTERS = {}
    L_REGISTERED_CONVERTERS = []


def get_converter(name=const.CONVERTER_DEFAULT, **converter_kwargs) -> base.FileConverter:
    """Get a FileConverter of the proper subclass for the requested converter type

    Parameters
    ----------
    name : str
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
    if name not in L_REGISTERED_CONVERTERS:
        raise base.FileConverterInputException(f"Converter {name} not recognized. Allowed converters are: " +
                                               f"{L_REGISTERED_CONVERTERS}")
    converter_class = D_REGISTERED_CONVERTERS[name]

    return converter_class(name=name, **converter_kwargs)


def run_converter(**converter_kwargs):
    """Shortcut to create and run a FileConverter in one step
    """

    return get_converter(**converter_kwargs).run()
