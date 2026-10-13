"""@file psdi_data_conversion/converters/script_template/converter.py

Atomsk file converter
"""

import os

from psdi_data_conversion.converters.base import FileConverterMeta, ScriptFileConverter
from psdi_data_conversion.utils import replace_extension


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

        # Delete the output file if it already exists
        if os.path.exists(self.out_filename):
            os.remove(self.out_filename)

        # Check if we might need to delete an intermediate file too, and do so
        intermediate_filename: str = ""
        backup_filename: str = ""
        if self.to_format_info.primary_name != self.to_format_name:
            intermediate_filename = replace_extension(self.out_filename, self.to_format_info.primary_name)
            if os.path.exists(intermediate_filename):
                # Create a name for a backup file we'll keep of any intermediate file if one exists
                backup_filename = intermediate_filename + ".bak"
                while os.path.exists(backup_filename):
                    backup_filename = backup_filename + ".bak"
                os.rename(intermediate_filename, backup_filename)

        try:

            super()._convert()

            # If an intermediate filename was used, rename to the proper output filename now
            if intermediate_filename:
                os.rename(intermediate_filename, self.out_filename)

        finally:
            if intermediate_filename and backup_filename:
                os.rename(backup_filename, intermediate_filename)


# Assign this converter to the `converter` variable - this lets the psdi_data_conversion.converter module detect and
# register it, making it available for use by the CLI and web app
converter = AtomskFileConverter
