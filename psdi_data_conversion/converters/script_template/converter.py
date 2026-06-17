"""@file psdi_data_conversion/converters/script_template/converter.py

ScriptTemplate file converter
"""

from psdi_data_conversion.converters.base import FileConverterMeta, ScriptFileConverter


class ScriptTemplateFileConverter(ScriptFileConverter):
    """File converter specialised to use ScriptTemplate for conversions"""

    meta = FileConverterMeta.load(__file__)

    script = ""
    required_bin = ""

    has_in_format_flags_or_options = False
    has_out_format_flags_or_options = False

    allowed_flags = ()
    allowed_options = ()

    def _convert(self):
        """Mandatory override: Perform the conversion"""
        pass

    def _create_message(self) -> str:
        """Optional override: Create a lot of options passed to the converter"""
        return ""


# Assign this converter to the `converter` variable - this lets the psdi_data_conversion.converter module detect and
# register it, making it available for use by the CLI and web app
converter = ScriptTemplateFileConverter
