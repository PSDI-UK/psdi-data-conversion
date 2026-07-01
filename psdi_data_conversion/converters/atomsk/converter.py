"""@file psdi_data_conversion/converters/script_template/converter.py

Atomsk file converter
"""

from psdi_data_conversion.converters.base import FileConverterMeta, ScriptFileConverter


class AtomskFileConverter(ScriptFileConverter):
    """File converter specialised to use Atomsk for conversions"""

    meta: FileConverterMeta = FileConverterMeta.load(__file__)

    script = ""
    required_bin = ""

    has_in_format_flags_or_options = False
    has_out_format_flags_or_options = False

    allowed_flags = ()
    allowed_options = ()

    def _create_message(self) -> str:
        """Optional override: Create a log of options passed to the converter"""
        return ""


# Assign this converter to the `converter` variable - this lets the psdi_data_conversion.converter module detect and
# register it, making it available for use by the CLI and web app
converter = AtomskFileConverter
