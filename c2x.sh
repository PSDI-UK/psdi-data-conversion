#!/bin/bash

#   c2x.sh
#   Version 1.0 10th January 2025

#   This shell script allows a c2x conversion to be carried out from the Python Flask app.

# Required positional arguments:
# arg1 is the fully-qualified input file name (e.g., /path/to/nacl.cif)
# arg2 is the fully-qualified output file name

cd bin
c2x ../$1 ../$2
cd ..

