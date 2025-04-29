"""@file psdi_data_conversion/converters/obenbabel.py

Created 2025-01-23 by Bryan Gillis.

Open Babel FileConverter
"""

from copy import deepcopy

import py
from openbabel import openbabel

from psdi_data_conversion.converters.base import FileConverter, FileConverterInputException
from psdi_data_conversion.security import SAFE_STRING_RE, string_is_safe

CONVERTER_OB = 'Open Babel'

# Constants related to command-line and library arguments unique to this converter
L_ALLOWED_COORD_GENS = ["Gen2D", "Gen3D", "neither"]
DEFAULT_COORD_GEN = "neither"
COORD_GEN_KEY = "coordinates"
L_ALLOWED_COORD_GEN_QUALS = ["fastest", "fast", "medium", "better", "best"]
DEFAULT_COORD_GEN_QUAL = "medium"
COORD_GEN_QUAL_KEY = "coordOption"

D_DEFAULT_OB_DATA = {COORD_GEN_KEY: DEFAULT_COORD_GEN,
                     COORD_GEN_QUAL_KEY: DEFAULT_COORD_GEN_QUAL}


def check_string_security(s: str):
    """Checks that a string is secure and raises an exception if it isn't.
    """
    if not string_is_safe(s):
        raise FileConverterInputException(f"Format option '{s}' does not pass security checks. It must pass the regex "
                                          f"/{SAFE_STRING_RE.pattern}/.", help=True)


def get_option_and_value(s: str):
    """Splits an option into the option character and value for it, checking for security
    """
    check_string_security(s)
    if len(s) == 0:
        return "", ""
    return s[0], s[1:]


def get_coord_gen(l_opts: list[str] | None) -> dict[str, str]:

    # Keyword arguments specific to OpenBabel conversion
    coord_gen: str
    if l_opts is None:
        coord_gen = DEFAULT_COORD_GEN
    else:
        coord_gen = l_opts[0]

    coord_gen_qual: str
    if l_opts is None or len(l_opts) == 1:
        coord_gen_qual = DEFAULT_COORD_GEN_QUAL
    else:
        coord_gen_qual = l_opts[1]

    # No more than two arguments supplied to --coord-gen
    if l_opts is not None and len(l_opts) > 2:
        raise FileConverterInputException("At most two arguments may be provided to --coord-gen, the mode and "
                                          "quality, e.g. '--coord-gen Gen3D best'", help=True)

    # Coordinate generation options are valid
    if coord_gen not in L_ALLOWED_COORD_GENS:
        raise FileConverterInputException(f"Coordinate generation type '{coord_gen}' not recognised. Allowed "
                                          f"types are: {L_ALLOWED_COORD_GENS}", help=True)
    if coord_gen_qual not in L_ALLOWED_COORD_GEN_QUALS:
        raise FileConverterInputException(f"Coordinate generation quality '{coord_gen_qual}' not recognised. "
                                          f"Allowed qualities are: {L_ALLOWED_COORD_GEN_QUALS}", help=True)

    return {COORD_GEN_KEY: coord_gen,
            COORD_GEN_QUAL_KEY: coord_gen_qual}


class OBFileConverter(FileConverter):
    """File Converter specialized to use Open Babel for conversions.

    This converter supports some additional configuration options which can be provided at class init or call to
    `run_converter()` through providing a dict to the `data` kwarg. The supported keys and values are:

    "from_flags": str
        String of concatenated one-letter flags for how to read the input file. To list the flags supported for a given
        input format, call ``psdi-data-convert -l -f <format> -w Open Babel`` at the command-line and look for the
        "Allowed input flags" section, if one exists, or alternatively call the library function
        ``psdi_data_conversion.database.get_in_format_args("Open Babel", <format>)`` from within Python code.

    "to_flags": str
        String of concatenated one-letter flags for how to write the output file. To list the flags supported for a
        given output format, call ``psdi-data-convert -l -f <format> -w Open Babel`` at the command-line and look for
        the "Allowed output flags" section, if one exists, or alternatively call the library function
        ``psdi_data_conversion.database.get_out_format_args("Open Babel", <format>)`` from within Python code.

    "from_options": str
        String of space-separated options for how to read the input file. Each option "word" in this string should start
        with the letter indicating which option is being used, followed by the value for that option. To list the
        options supported for a given input format, call ``psdi-data-convert -l -f <format> -w Open Babel`` at the
        command-line and look for the "Allowed input options" section, if one exists, or alternatively call the library
        function ``psdi_data_conversion.database.get_in_format_args("Open Babel", <format>)`` from within Python code.

    "to_options": str
        String of space-separated options for how to write the output file. Each option "word" in this string should
        start with the letter indicating which option is being used, followed by the value for that option. To list the
        options supported for a given output format, call ``psdi-data-convert -l -t <format> -w Open Babel`` at the
        command-line and look for the "Allowed output options" section, if one exists, or alternatively call the library
        function ``psdi_data_conversion.database.get_out_format_args("Open Babel", <format>)`` from within Python code.

    "coordinates": str
        One of "Gen2D", "Gen3D", or "neither", specifying how positional coordinates should be generated in the output
        file. Default "neither"

    "coordOption": str
        One of "fastest", "fast", "medium", "better", or "best", specifying the quality of the calculation of
        coordinates. Default "medium"

    Note that some other keys are supported for compatibility purposes, but these may be deprecated in the future.
    """

    name = CONVERTER_OB
    has_in_format_flags_or_options = True
    has_out_format_flags_or_options = True
    database_key_prefix = "ob"

    allowed_flags = ()
    allowed_options = (("--coord-gen",
                        {"help": "(Open Babel converter only). The mode to be used for Open Babel calculation of "
                         "atomic coordinates, and optionally the quality of the conversion. The mode should be one of "
                         "'Gen2D', 'Gen3D', or 'neither' (default 'neither'). The quality, if supplied, should be "
                         "one of 'fastest', 'fast', 'medium', 'better' or 'best' (default 'medium'). E.g. "
                         "'--coord-gen Gen2D' (quality defaults to 'medium'), '--coord-gen Gen3D best'",
                         "type": str,
                         "default": None,
                         "nargs": "+"},
                        get_coord_gen),)

    def _convert(self):

        self.logger.debug("Using OpenBabel's Python library to perform file conversion")

        # Apply default values to the data dict
        _data = deepcopy(D_DEFAULT_OB_DATA)
        _data.update(self.data)
        self.data = _data

        try:
            stdouterr_ob = py.io.StdCaptureFD(in_=False)

            ob_conversion = openbabel.OBConversion()
            ob_conversion.SetInAndOutFormats(self.from_format_info.name, self.to_format_info.name)

            # Retrieve 'from' and 'to' option flags and arguments
            from_flags = self.data.get("from_flags", "")
            to_flags = self.data.get("to_flags", "")
            from_arg_flags = self.data.get("from_arg_flags", "")
            to_arg_flags = self.data.get("to_arg_flags", "")
            from_args = self.data.get("from_args", "")
            to_args = self.data.get("to_args", "")

            # Add option flags and arguments as appropriate
            for char in from_flags:
                check_string_security(char)
                ob_conversion.AddOption(char, ob_conversion.INOPTIONS)

            for char in to_flags:
                check_string_security(char)
                ob_conversion.AddOption(char, ob_conversion.OUTOPTIONS)

            self.data["read_flags_args"] = []
            self.data["write_flags_args"] = []

            # Check if we were provided options by the command-line script/library or the web app, and handle them
            # appropriately
            if "from_options" in self.data:
                # From options were provided by the command-line script or library
                l_from_options = self.data["from_options"].split()
                for opt in l_from_options:
                    option, value = get_option_and_value(opt)
                    ob_conversion.AddOption(option, ob_conversion.INOPTIONS, value)
                self.logger.debug(f"Set Open Babel read flags arguments to: {self.data['from_options']}")
                # Store the options in the "read_flags_args" entry for the later logging
                self.data["read_flags_args"] = l_from_options
            else:
                # From options were provided by the command-line script or library
                for char in from_arg_flags:
                    index = from_args.find('£')
                    arg, from_args = from_args[0:index], from_args[index + 1:len(from_args)]
                    check_string_security(char), check_string_security(arg)
                    ob_conversion.AddOption(char, ob_conversion.INOPTIONS, arg)
                    self.data["read_flags_args"].append(char + "  " + arg)
                self.logger.debug(f"Set Open Babel read flags arguments to: {self.data['read_flags_args']}")

            if "to_options" in self.data:
                # From options were provided by the command-line script or library
                l_to_options = self.data["to_options"].split()
                for opt in l_to_options:
                    option, value = get_option_and_value(opt)
                    ob_conversion.AddOption(option, ob_conversion.OUTOPTIONS, value)
                self.logger.debug(f"Set Open Babel write flags arguments to: {self.data['to_options']}")
                # Store the options in the "write_flags_args" entry for the later logging
                self.data["write_flags_args"] = l_to_options
            else:
                # From options were provided by the command-line script or library
                for char in to_arg_flags:
                    index = to_args.find('£')
                    arg, to_args = to_args[0:index], to_args[index + 1:len(to_args)]
                    check_string_security(char), check_string_security(arg)
                    ob_conversion.AddOption(char, ob_conversion.OUTOPTIONS, arg)
                    self.data["write_flags_args"].append(char + "  " + arg)
                self.logger.debug(f"Set Open Babel write flags arguments to: {self.data['read_flags_args']}")

            # Read the file to be converted
            mol = openbabel.OBMol()
            ob_conversion.ReadFile(mol, self.in_filename)

            # Calculate atomic coordinates
            if self.data[COORD_GEN_KEY] == 'neither':
                self.data[COORD_GEN_QUAL_KEY] = 'N/A'
            else:
                # Retrieve coordinate calculation option (fastest, fast, medium, better, best)
                self.option = self.data[COORD_GEN_QUAL_KEY]

                gen = openbabel.OBOp.FindType(self.data[COORD_GEN_KEY])
                self.logger.debug(f"Performing Open Babel {self.data[COORD_GEN_KEY]} coordinate conversion with option "
                                  f"'{self.option}'")
                gen.Do(mol, self.data[COORD_GEN_QUAL_KEY])

            # Write the converted file
            ob_conversion.WriteFile(mol, self.out_filename)

            self.out, self.err = stdouterr_ob.reset()   # Grab stdout and stderr

        finally:
            # Reset stdout and stderr capture
            stdouterr_ob.done()

        if "Open Babel Error" in self.err:
            self._abort_from_err()

    def _create_message(self) -> str:
        """Overload method to create a log of options passed to the converter
        """

        message = ""

        label_length = 19

        for (label, key, multi) in (("Coord. gen.:", COORD_GEN_KEY, False),
                                    ("Coord. option:", "coord_option", False),
                                    ("Read options:", "from_flags", False),
                                    ("Write options:", "to_flags", False),
                                    ("Read opts + args:", "read_flags_args", True),
                                    ("Write opts + args:", "write_flags_args", True)):
            val = self.data.get(key)

            if not val:
                message += f"{label:<{label_length}}none\n"
                continue

            if multi:
                l_items = val
            else:
                l_items = (val,)

            for i, item in enumerate(l_items):
                if i == 0:
                    line_label = label
                else:
                    line_label = ""
                message += f"{line_label:<{label_length}}{item}\n"

        return message


# Assign this converter to the `converter` variable - this lets the psdi_data_conversion.converter module detect and
# register it, making it available for use by the CLI and web app
converter = OBFileConverter
