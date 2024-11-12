#!/bin/bash

#   atomsk.sh
#   Version 1.0 29th October 2024

#   This shell script allows an Atomsk conversion to be carried out from the Python Flask app.

# Required positional arguments:
# arg1 is the full inlet file name (e.g., nacl.cif)
# arg2 is the file name, not including its extension (e.g., nacl)
# arg3 is the outlet file's extension without the dot (e.g., pdb)

bin/atomsk <<EOD
read static/uploads/$1
write static/downloads/$2.$3
clear
quit
EOD
