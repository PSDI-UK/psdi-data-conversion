"""@file psdi_data_conversion/converters/obenbabel.py

Created 2025-01-23 by Bryan Gillis.

Open Babel FileConverter
"""

from openbabel import openbabel
import py

from psdi_data_conversion.converters.base import FileConverter

CONVERTER_OB = 'Open Babel'


class OBFileConverter(FileConverter):
    """File Converter specialized to use Open Babel for conversions
    """

    name = CONVERTER_OB

    def _convert(self):
        stdouterr_ob = py.io.StdCaptureFD(in_=False)

        ob_conversion = openbabel.OBConversion()
        ob_conversion.SetInAndOutFormats(self.from_format, self.to_format)

        # Retrieve 'from' and 'to' option flags and arguments
        from_flags = self.data.get("from_flags", "")
        to_flags = self.data.get("to_flags", "")
        from_arg_flags = self.data.get("from_arg_flags", "")
        to_arg_flags = self.data.get("to_arg_flags", "")
        from_args = self.data.get("from_args", "")
        to_args = self.data.get("to_args", "")

        # Add option flags and arguments as appropriate
        for char in from_flags:
            ob_conversion.AddOption(char, ob_conversion.INOPTIONS)

        for char in to_flags:
            ob_conversion.AddOption(char, ob_conversion.OUTOPTIONS)

        self.data["read_flags_args"] = []
        self.data["write_flags_args"] = []

        for char in from_arg_flags:
            index = from_args.find('£')
            arg = from_args[0:index]
            from_args = from_args[index + 1:len(from_args)]
            self.data["read_flags_args"].append(char + "  " + arg)
            ob_conversion.AddOption(char, ob_conversion.INOPTIONS, arg)

        for char in to_arg_flags:
            index = to_args.find('£')
            arg = to_args[0:index]
            to_args = to_args[index + 1:len(to_args)]
            self.data["write_flags_args"].append(char + "  " + arg)
            ob_conversion.AddOption(char, ob_conversion.OUTOPTIONS, arg)

        # Read the file to be converted
        mol = openbabel.OBMol()
        ob_conversion.ReadFile(mol, self.in_filename)

        # Calculate atomic coordinates
        if self.data["coordinates"] == 'neither':
            self.data['coordOption'] = 'N/A'
        else:
            # Retrieve coordinate calculation option (fastest, fast, medium, better, best)
            self.option = self.data['coordOption']

            gen = openbabel.OBOp.FindType(self.data["coordinates"])
            gen.Do(mol, self.data['coordOption'])

        # Write the converted file
        ob_conversion.WriteFile(mol, self.out_filename)

        self.out, self.err = stdouterr_ob.reset()   # Grab stdout and stderr
        stdouterr_ob.done()

        if "Open Babel Error" in self.err:
            self._abort_from_err()

    def _create_message(self) -> str:
        """Overload method to create a log of options passed to the converter
        """

        message = super()._create_message()

        coordinates = self.data.get("coordinates", "none")
        if coordinates and coordinates != "neither":
            message += 'Coord. gen.:       ' + coordinates + '\n'

        coord_option = self.data.get("coordOption", "")
        if coord_option:
            message += 'Coord. option:     ' + coord_option + '\n'

        from_flags = self.data.get("coordOption", "")
        if from_flags:
            message += 'Read options:      ' + from_flags + '\n'

        to_flags = self.data.get("to_flags", "")
        if to_flags:
            message += 'Write options:     ' + to_flags + '\n'

        read_flags_args = self.data.get("read_flags_args", "")
        if read_flags_args:
            for pair in read_flags_args:
                message += 'Read opts + args:  ' + pair + '\n'

        write_flags_args = self.data.get("write_flags_args", "")
        if write_flags_args:
            for pair in write_flags_args:
                message += 'Write opts + args: ' + pair + '\n'

        return message


# Assign this converter to the `converter` variable - this lets the psdi_data_conversion.converter module detect and
# register it, making it available for use by the CLI and web app
converter = OBFileConverter
