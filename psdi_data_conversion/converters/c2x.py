"""@file psdi_data_conversion/converters/c2x.py

Created 2025-01-23 by Bryan Gillis.

C2X FileConverter
"""

import os
import subprocess

from psdi_data_conversion.converters.base import FileConverter


class C2xFileConverter(FileConverter):
    """File Converter specialized to use Atomsk for conversions
    """

    def _convert(self):
        atomsk = subprocess.run(['sh', 'psdi_data_conversion/scripts/c2x.sh',
                                 self.in_filename, self.out_filename], capture_output=True, text=True)

        self.out = atomsk.stdout
        self.err = atomsk.stderr

        if self.err.find('Error') > -1:
            self._abort_from_err()

        self.in_size, self.out_size = self._check_file_size()

        if self.file_to_convert != 'file':   # Website only (i.e., not command line option)
            if self.delete_input:
                os.remove(self.in_filename)
            self.from_format = self.form['from_full']
            self.to_format = self.form['to_full']
            self.quality = self.form['success']
        else:
            self.quality = self.get_quality(self.from_format, self.to_format)

        self._log_success()
