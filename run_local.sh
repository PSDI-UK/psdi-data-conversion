#!/bin/bash

# Execute a local run of the application from the proper path

PACKAGE_PATH=`python -c "import psdi_data_conversion; print(psdi_data_conversion.__path__[0])"`
cd $PACKAGE_PATH/..
python -m flask --app psdi_data_conversion/app.py run