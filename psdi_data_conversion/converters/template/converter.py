"""@file psdi_data_conversion/converters/template/converter.py

Template file converter
"""

from psdi_data_conversion.converters.base import FileConverter, FileConverterMeta


class TemplateFileConverter(FileConverter):
    """File converter specialised to use Template for conversions"""

    meta = FileConverterMeta.load(__file__)
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
converter = TemplateFileConverter
