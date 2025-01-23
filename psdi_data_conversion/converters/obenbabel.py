"""@file psdi_data_conversion/converters/obenbabel.py

Created 2025-01-23 by Bryan Gillis.

Open Babel FileConverter
"""

import openbabel
import os
import py

from psdi_data_conversion.converters.base import FileConverter


class OBFileConverter(FileConverter):
    """File Converter specialized to use Open Babel for conversions
    """

    def _convert(self):
        stdouterr_ob = py.io.StdCaptureFD(in_=False)

        ob_conversion = openbabel.OBConversion()
        ob_conversion.SetInAndOutFormats(self.from_format, self.to_format)

        # Retrieve 'from' and 'to' option flags and arguments
        self.from_flags = self.form['from_flags']
        self.to_flags = self.form['to_flags']
        from_arg_flags = self.form['from_arg_flags']
        to_arg_flags = self.form['to_arg_flags']
        from_args = self.form['from_args']
        to_args = self.form['to_args']

        # Add option flags and arguments as appropriate
        for char in self.from_flags:
            ob_conversion.AddOption(char, ob_conversion.INOPTIONS)

        for char in self.to_flags:
            ob_conversion.AddOption(char, ob_conversion.OUTOPTIONS)

        self.read_flags_args = []
        self.write_flags_args = []

        for char in from_arg_flags:
            index = from_args.find('£')
            arg = from_args[0:index]
            from_args = from_args[index + 1:len(from_args)]
            self.read_flags_args.append(char + "  " + arg)
            ob_conversion.AddOption(char, ob_conversion.INOPTIONS, arg)

        for char in to_arg_flags:
            index = to_args.find('£')
            arg = to_args[0:index]
            to_args = to_args[index + 1:len(to_args)]
            self.write_flags_args.append(char + "  " + arg)
            ob_conversion.AddOption(char, ob_conversion.OUTOPTIONS, arg)

        # Read the file to be converted
        mol = openbabel.OBMol()
        ob_conversion.ReadFile(mol, self.in_filename)

        # Retrieve coordinate calculation type (Gen2D, Gen3D, neither)
        self.calc_type = self.form['coordinates']

        self.option = 'N/A'

        # Calculate atomic coordinates
        if self.calc_type != 'neither':
            # Retrieve coordinate calculation option (fastest, fast, medium, better, best)
            self.option = self.form['coordOption']

            gen = openbabel.OBOp.FindType(self.calc_type)
            gen.Do(mol, self.option)

        # Write the converted file
        ob_conversion.WriteFile(mol, self.out_filename)

        self.out, self.err = stdouterr_ob.reset()   # Grab stdout and stderr
        stdouterr_ob.done()

        self.in_size, self.out_size = self._check_file_size()

        if self.file_to_convert != 'file':  # Website only (i.e., not command line option)
            if self.delete_input:
                os.remove(self.in_filename)
            self.from_format = self.form['from_full']
            self.to_format = self.form['to_full']
            self.quality = self.form['success']
        else:
            self.quality = self.get_quality(self.from_format, self.to_format)

        if self.err.find('Error') > -1:
            self._abort_from_err()
        else:
            self._log_success()
