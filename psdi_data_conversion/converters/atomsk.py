"""@file psdi_data_conversion/converters/atomsk.py

Created 2025-01-23 by Bryan Gillis.

Atomsk FileConverter
"""

from psdi_data_conversion.converters.base import ScriptFileConverter

CONVERTER_ATO = 'Atomsk'


class AtoFileConverter(ScriptFileConverter):
    """File Converter specialized to use Atomsk for conversions.

    This converter does not yet support any additional configuration options provided at class init to the `data` kwarg.
    """

    name = CONVERTER_ATO
    script = "atomsk.sh"
    required_bin = "atomsk"
    info = ("Atomsk binaries compiled for 64-bit Linux and MacOS systems are distributed with this package. It may be "
            "registered on other systems by compiling it locally and adding the compiled 'atomsk' binary (with this "
            "exact name - rename it or make a symbolic link to it if necessary) to your $PATH.\n"
            "\n"
            "Atomsk is licensed under GPLv3, the full text of which may be found at "
            "https://www.gnu.org/licenses/gpl-3.0.en.html. Its binaries are redistributed here under the terms of this "
            "license, and any further redistribution must also follow these terms. Its corresponding source code "
            "may be found at https://github.com/pierrehirel/atomsk/")


# Assign this converter to the `converter` variable - this lets the psdi_data_conversion.converter module detect and
# register it, making it available for use by the command-line script, python library, and web app
converter = AtoFileConverter
