"""@file psdi_data_conversion/converters/atomsk.py

Created 2025-01-23 by Bryan Gillis.

Atomsk FileConverter
"""

from psdi_data_conversion.constants import CONVERTER_ATO
from psdi_data_conversion.converters.base import ScriptFileConverter


class AtoFileConverter(ScriptFileConverter):
    """File Converter specialized to use Atomsk for conversions
    """

    converter = CONVERTER_ATO
    script = "atomsk.sh"


converter = AtoFileConverter
