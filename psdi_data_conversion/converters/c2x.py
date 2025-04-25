"""@file psdi_data_conversion/converters/c2x.py

Created 2025-01-23 by Bryan Gillis.

c2x FileConverter
"""

from psdi_data_conversion.converters.base import ScriptFileConverter

CONVERTER_C2X = 'c2x'


class C2xFileConverter(ScriptFileConverter):
    """File Converter specialized to use c2x for conversions.

    This converter does not yet support any additional configuration options provided at class init to the `data` kwarg.
    """

    name = CONVERTER_C2X
    script = "c2x.sh"
    required_bin = "c2x"
    info = ("c2x binaries compiled for 64-bit Linux and MacOS systems are distributed with this package. It may be "
            "registered on other systems by compiling it locally and adding the compiled 'c2x' binary (with this "
            "exact name - rename it or make a symbolic link to it if necessary) to your $PATH.\n"
            "\n"
            "c2x is licensed under GPLv3, the full text of which may be found at "
            "https://www.gnu.org/licenses/gpl-3.0.en.html. Its binaries are redistributed here under the terms of this "
            "license, and any further redistribution must also follow these terms. Its corresponding source code "
            "may be downloaded from https://www.c2x.org.uk/downloads/")
    supports_ambiguous_extensions = True

    def _get_script_args(self):
        """Override the standard script arguments so we can set the different format names expected by c2x
        """
        l_script_args = super()._get_script_args()

        # Update the output format to c2x style
        l_script_args[0] = "--" + self.to_format_info.c2x_format

        # TODO - check if the input file has an extension which will be accepted by c2x for its format, and handle if
        # not

        return l_script_args


# Assign this converter to the `converter` variable - this lets the psdi_data_conversion.converter module detect and
# register it, making it available for use by the command-line script, python library, and web app
converter = C2xFileConverter
