"""@file psdi_data_conversion/converters/c2x.py

Created 2025-01-23 by Bryan Gillis.

C2X FileConverter
"""

from psdi_data_conversion.converters.base import ScriptFileConverter

CONVERTER_C2X = 'c2x'


class C2xFileConverter(ScriptFileConverter):
    """File Converter specialized to use C2X for conversions
    """

    name = CONVERTER_C2X
    script = "c2x.sh"


# Assign this converter to the `converter` variable - this lets the psdi_data_conversion.converter module detect and
# register it, making it available for use by the CLI and web app
converter = C2xFileConverter
