"""@file psdi_data_conversion/converters/script_template/converter.py

Atomsk file converter
"""

import os

from psdi_data_conversion.converters.base import FileConverterMeta, ScriptFileConverter


def replace_extension(filename: str, new_ext: str):
    """Replace the extension of a filename with a new extension"""

    # Strip all extensions from the filename
    stripped_filename = filename
    last_stripped_filename = None
    while stripped_filename != last_stripped_filename:
        last_stripped_filename = stripped_filename
        stripped_filename, _ = os.path.splitext(stripped_filename)

    # Get the extension with a full-stop at the beginning
    if new_ext.startswith("."):
        ext = new_ext
    else:
        ext = "." + new_ext

    return stripped_filename + ext


class AtomskFileConverter(ScriptFileConverter):
    """File converter specialised to use Atomsk for conversions"""

    meta: FileConverterMeta = FileConverterMeta.load(__file__)

    script = "atomsk.sh"
    required_bin = "atomsk"

    has_in_format_flags_or_options = False
    has_out_format_flags_or_options = False

    allowed_flags = ()
    allowed_options = ()

    def _get_script_args(self):
        """Override the standard script arguments so we can set the different format names expected by atomsk
        """
        l_script_args = super()._get_script_args()

        # Update the formats to Atomsk style if relevant
        if self.from_format_info.name != self.from_format_name:
            l_script_args[1] = replace_extension(l_script_args[1], self.from_format_name)
        if self.to_format_info.name != self.to_format_name:
            l_script_args[2] = replace_extension(l_script_args[2], self.to_format_name)

        return l_script_args

    def _convert(self):
        """Atomsk has a bug with some formats which results in the file being created with the wrong extension, so
        we check for that here and correct it
        """

        super()._convert()

        if self.to_format_info.primary_name != self.to_format_name:
            incorrect_out_filename = replace_extension(self.out_filename, self.to_format_info.primary_name)
            os.rename(incorrect_out_filename, self.out_filename)


# Assign this converter to the `converter` variable - this lets the psdi_data_conversion.converter module detect and
# register it, making it available for use by the CLI and web app
converter = AtomskFileConverter
