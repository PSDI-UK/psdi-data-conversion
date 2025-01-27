#!/bin/bash

# The envvar MAX_FILESIZE can be used to set the maximum allowed filesize in MB - 0 indicates no maximum
if [ -z $MAX_FILESIZE ]; then
  export MAX_FILESIZE=0
fi

export LOGGING=Full

# Uncomment the following line to enable debug mode
# export FLASK_ENV=development

# Execute a local run of the application from the proper path

PACKAGE_PATH=`python -c "import psdi_data_conversion; print(psdi_data_conversion.__path__[0])"`
cd $PACKAGE_PATH/..
python -m flask --app psdi_data_conversion/app.py run