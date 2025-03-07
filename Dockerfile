#
# Dockerfile for image containerising PSDI's data conversion service
#
#
# Building the image
# ------------------
#
# 1. Download the service repo containing the source code
#    from https://github.com/PSDI-UK/psdi-data-conversion/tree/main
# 2. Copy this file into the main directory of the repo
# 3. From within the main directory of the repo, use docker to build the image
#    'psdi-data-conversion' via the command
#    ``docker build -t psdi-data-conversion .``
#
# To check that the image has been build run ``docker images``, which should
# list and image called 'psdi-data-conversion'.
#
#
# Running the service
# -------------------
#
# The command ``docker run -p 8000:8000 psdi-data-conversion`` will run the
# service on port 8000 of localhost, with logs output to stdout. To access the
# service visit http://localhost:8000 in your browser.
#

FROM python:3.12-slim-bookworm

RUN apt update
RUN apt-get -y install libxrender1 libxext6

# Install Python packages (including openbabel-wheel)
RUN pip install --upgrade pip
RUN pip install gunicorn
RUN pip install flask

WORKDIR /app
COPY requirements.txt /app
COPY psdi_data_conversion /app/psdi_data_conversion

RUN pip install -r requirements.txt

ENV PYTHONPATH="."
ENV SERVICE_MODE=true
ENV MAX_FILESIZE=1
ENV LOG_MODE=full

# Set LOG_LEVEL to a desired level (e.g. "debug") to force all logging to be at that level
ENV LOG_LEVEL= 

EXPOSE 8000

RUN mkdir /app/psdi_data_conversion/static/uploads
RUN mkdir /app/psdi_data_conversion/static/downloads

CMD ["gunicorn", "-b", "0.0.0.0:8000", "psdi_data_conversion.app:app"]
