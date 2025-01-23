"""@file psdi_data_conversion/converters/atomsk.py

Created 2025-01-23 by Bryan Gillis.

Atomsk FileConverter
"""

from psdi_data_conversion.converters.base import ScriptFileConverter


class AtoFileConverter(ScriptFileConverter):
    """File Converter specialized to use Atomsk for conversions
    """

    script = "atomsk.sh"
