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
# The command ``docker run -p 5000:5000 psdi-data-conversion`` will run the
# service on port 5000 of localhost, with logs output to stdout. To access the
# service visit http://localhost:5000 in your browser.
#

FROM python:3.12-slim-bookworm
WORKDIR /app
COPY . .

# Install Python packages (including openbabel-wheel)
RUN pip install --upgrade pip
RUN pip install -r requirements.txt
RUN pip install flask
RUN pip install gunicorn

RUN apt update
RUN apt-get -y install libxrender1 libxext6

ENV FLASK_APP=app.py
ENV FLASK_RUN_HOST=0.0.0.0
ENV PYTHONPATH="."

EXPOSE 5000

CMD ["/bin/bash", "-c", "gunicorn app:app & python -m flask run"]

