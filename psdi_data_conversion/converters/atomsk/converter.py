"""@file psdi_data_conversion/converters/script_template/converter.py

Atomsk file converter
"""

import os

from psdi_data_conversion.converters.base import FileConverterMeta, ScriptFileConverter


class AtomskFileConverter(ScriptFileConverter):
    """File converter specialised to use Atomsk for conversions"""

    meta: FileConverterMeta = FileConverterMeta.load(__file__)

    script = "atomsk.sh"
    required_bin = "atomsk"

    has_in_format_flags_or_options = False
    has_out_format_flags_or_options = False

    allowed_flags = ()
    allowed_options = ()

    def _convert(self):
        """Atomsk has a bug with some formats which results in the file being created with the wrong extension, so
        we check for that here and correct it
        """
        super()._convert()

        if self.to_format_info.name in ["sxyz", "exyz"]:
            incorrect_out_filename = self.out_filename[:-4] + "xyz"
            os.rename(incorrect_out_filename, self.out_filename)


# Assign this converter to the `converter` variable - this lets the psdi_data_conversion.converter module detect and
# register it, making it available for use by the CLI and web app
converter = AtomskFileConverter
