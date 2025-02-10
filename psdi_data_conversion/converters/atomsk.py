"""@file psdi_data_conversion/converters/atomsk.py

Created 2025-01-23 by Bryan Gillis.

Atomsk FileConverter
"""

from psdi_data_conversion.converters.base import ScriptFileConverter

CONVERTER_ATO = 'Atomsk'


class AtoFileConverter(ScriptFileConverter):
    """File Converter specialized to use Atomsk for conversions
    """

    name = CONVERTER_ATO
    script = "atomsk.sh"
    info = "See documentation for Atomsk at https://atomsk.univ-lille.fr/"


# Assign this converter to the `converter` variable - this lets the psdi_data_conversion.converter module detect and
# register it, making it available for use by the CLI and web app
converter = AtoFileConverter
