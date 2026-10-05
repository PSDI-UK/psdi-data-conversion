"""@file test_data/test_converter_simple/converter.py

Test file converter, a simple case
"""

from psdi_data_conversion.converters.base import FileConverter, FileConverterMeta


class TestFileConverterSimple(FileConverter):
    """File converter specialised to use Example for conversions"""

    meta: FileConverterMeta = FileConverterMeta.load(__file__)
    has_in_format_flags_or_options = False
    has_out_format_flags_or_options = False

    allowed_flags = ()
    allowed_options = ()

    def _convert(self):
        """Mandatory override: Perform the conversion"""
        pass


# Assign this converter to the `converter` variable - this lets the psdi_data_conversion.converter module detect and
# register it, making it available for use by the CLI and web app
converter = TestFileConverterSimple
