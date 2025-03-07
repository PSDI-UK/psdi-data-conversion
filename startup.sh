#!/bin/sh

apt update
apt-get install -y libxrender1 libxext6

export SERVICE_MODE=true
export LOG_MODE=full
export LOG_LEVEL=
export MAX_FILESIZE=1 # Megabyte

gunicorn psdi_data_conversion.app:app
