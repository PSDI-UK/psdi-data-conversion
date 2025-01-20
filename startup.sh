#!/bin/sh

apt update
apt-get install -y libxrender1 libxext6

export AUTH=true

gunicorn psdi_data_conversion.app:app
