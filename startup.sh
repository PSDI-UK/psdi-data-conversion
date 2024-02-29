#!/bin/sh

apt update
apt-get install -y libxrender1 libxext6

gunicorn app:app
