#!/bin/bash

# The envvar MAX_FILESIZE can be used to set the maximum allowed filesize in MB - 0 indicates no maximum
if [ -z $MAX_FILESIZE ]; then
  export MAX_FILESIZE=0
fi

# Logging control - "full" sets server-style logging, which is necessary to produce the output logs with the expected
# names. This should not be changed, or else errors will occur
export LOG_MODE=full

# The level to log at. Leave blank for defaults, which logs at INFO level for user output and ERROR level for the server
# log and stdout. If set to a different value here (e.g. DEBUG), all these channels will be set to that level
export LOG_LEVEL=

# The envvar SERVICE_MODE can be set to "true" to make this behave as if it's running as the public web service -
# uncomment the following line to enable that
# export SERVICE_MODE=true

# The envvar PRODUCTION_MODE is set to "true" to hide any dev-only aspects from the GUI, such as the label of the
# latest SHA. Comment out the following line to disable this behaviour and expose those elements
export PRODUCTION_MODE=true

# Uncomment the following line to enable debug mode
# export FLASK_ENV=development

# Execute a local run of the application from the proper path

PACKAGE_PATH=`python -c "import psdi_data_conversion; print(psdi_data_conversion.__path__[0])"`
cd $PACKAGE_PATH/..
python -m flask --app psdi_data_conversion/app.py run