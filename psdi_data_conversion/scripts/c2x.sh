#!/bin/bash

#   c2x.sh
#   Version 1.0 10th January 2025

#   This shell script allows a c2x conversion to be carried out from the Python Flask app.

# Required positional arguments:
# arg1 is the output format preceded by two dashes (e.g., --cml)
# arg2 is the fully-qualified input file name (e.g., /path/to/nacl.cif)
# arg3 is the fully-qualified output file name

DEFAULT_DIST=linux

# The ennvar DIST can be used to set the distribution, indicating the subdirectory to search in for the binary
if [ -z $DIST ]; then
  DIST=$DEFAULT_DIST
fi

# If the envvar BIN_PATH is set, we use that path directly. Otherwise we use the binary distributed with this package
if [ -z $BIN_PATH ]; then
  BIN_PATH=psdi_data_conversion/bin/$DIST/atomsk
fi

$BIN_PATH $1 $2 $3 $4 $5

