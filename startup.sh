#!/bin/sh

apt update
apt-get install -y libxrender1 libxext6

export SERVICE_MODE=true
export LOGGING=full
export MAX_FILESIZE=1 # Megabyte

gunicorn psdi_data_conversion.app:app
