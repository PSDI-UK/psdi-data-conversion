#!/bin/bash

# The envvar MAX_FILESIZE can be used to set the maximum allowed filesize in MB - 0 indicates no maximum
if [ -z $MAX_FILESIZE ]; then
  export MAX_FILESIZE=0
fi

# The envvar LOGGING can be used to set how logs are stored. Allowed values are:
# Full - multi-file logging, only recommended when running as a public web app
# Simple - logs saved to one file
# None - output only to stdout/stderr
if [ -z $LOGGING ]; then
  export LOGGING=Simple
fi

# Uncomment the following line to enable debug mode
# export FLASK_ENV=development

# Execute a local run of the application from the proper path

PACKAGE_PATH=`python -c "import psdi_data_conversion; print(psdi_data_conversion.__path__[0])"`
cd $PACKAGE_PATH/..
python -m flask --app psdi_data_conversion/app.py run