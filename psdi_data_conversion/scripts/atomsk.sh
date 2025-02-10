#!/bin/bash

#   atomsk.sh
#   Version 1.1 17th December 2024

#   This shell script allows an Atomsk conversion to be carried out from the Python Flask app.

# Required positional arguments:
# arg2 is the fully-qualified input file name (e.g., /path/to/nacl.cif)
# arg3 is the fully-qualified output file name

psdi_data_conversion/bin/atomsk <<EOD
read $2
write $3
clear
quit
EOD

# Cleanup intermediate files that might have been created
rm -f quit