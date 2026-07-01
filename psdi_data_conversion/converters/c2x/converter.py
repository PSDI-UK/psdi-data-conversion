"""@file psdi_data_conversion/converters/script_template/converter.py

c2x file converter
"""

from psdi_data_conversion.converters.base import FileConverterMeta, ScriptFileConverter


class C2xFileConverter(ScriptFileConverter):
    """File converter specialised to use c2x for conversions"""

    meta: FileConverterMeta = FileConverterMeta.load(__file__)

    script = "c2x.sh"
    required_bin = "c2x"

    has_in_format_flags_or_options = False
    has_out_format_flags_or_options = False

    allowed_flags = ()
    allowed_options = ()

    def _get_script_args(self):
        """Override the standard script arguments so we can set the different format names expected by c2x
        """
        l_script_args = super()._get_script_args()

        # Update the output format to c2x style
        l_script_args[0] = "--" + self.to_format_info.c2x_format

        # TODO - check if the input file has an extension which will be accepted by c2x for its format, and handle if
        # not

        return l_script_args

    def _create_message(self) -> str:
        """Optional override: Create a log of options passed to the converter"""
        return ""


# Assign this converter to the `converter` variable - this lets the psdi_data_conversion.converter module detect and
# register it, making it available for use by the CLI and web app
converter = C2xFileConverter
