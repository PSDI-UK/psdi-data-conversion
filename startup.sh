#!/bin/sh

apt update
apt-get install -y libxrender1 libxext6

gunicorn psdi_data_conversion.app:app
