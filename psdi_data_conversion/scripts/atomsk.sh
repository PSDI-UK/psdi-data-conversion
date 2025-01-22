#!/bin/bash

#   atomsk.sh
#   Version 1.1 17th December 2024

#   This shell script allows an Atomsk conversion to be carried out from the Python Flask app.

# Required positional arguments:
# arg1 is the fully-qualified input file name (e.g., /path/to/nacl.cif)
# arg2 is the fully-qualified output file name

psdi_data_conversion/bin/atomsk <<EOD
read $1
write $2
clear
quit
EOD
