"""@file psdi_data_conversion/converters/atomsk.py

Created 2025-01-23 by Bryan Gillis.

Atomsk FileConverter
"""

import subprocess

from psdi_data_conversion.converters.base import FileConverter


class AtoFileConverter(FileConverter):
    """File Converter specialized to use Atomsk for conversions
    """

    def _convert(self):
        atomsk = subprocess.run(['sh', 'psdi_data_conversion/scripts/atomsk.sh',
                                 self.in_filename, self.out_filename], capture_output=True, text=True)

        self.out = atomsk.stdout
        self.err = atomsk.stderr
