"""@file psdi_data_conversion/converters/example/converter.py

Example file converter
"""

from psdi_data_conversion.converters.base import FileConverter, FileConverterArgException, FileConverterMeta


def process_example_option(l_opts: list[str] | None) -> dict[str, str]:
    """Example method to process an option for this converter"""
    if l_opts is None:
        return {"opt": "None"}
    if len(l_opts) == 1:
        return {"opt": l_opts[0]}
    if len(l_opts) == 2:
        return {l_opts[0]: l_opts[1]}
    raise FileConverterArgException("Invalid number of arguments for option '--example-opt'")


class ExampleFileConverter(FileConverter):
    """File converter specialised to use Example for conversions"""

    meta = FileConverterMeta.load(__file__)
    has_in_format_flags_or_options = True
    has_out_format_flags_or_options = True
    database_key_prefix = "ob"

    allowed_flags = (("--example-flag-1",
                      "(Example converter only). Help text for example flag 1."),
                     ("--example-flag-2",
                     "(Example converter only). Help text for example flag 2."))
    allowed_options = (("--example-opt",
                        {"help": "(Example converter only). Help text for example option.",
                         "type": str,
                         "default": None,
                         "nargs": "+"},
                        process_example_option))

    def _convert(self):
        """Mandatory override: Perform the conversion"""
        pass

    def _create_message(self) -> str:
        """Optional override: Create a lot of options passed to the converter"""
        return ""


# Assign this converter to the `converter` variable - this lets the psdi_data_conversion.converter module detect and
# register it, making it available for use by the CLI and web app
converter = ExampleFileConverter
