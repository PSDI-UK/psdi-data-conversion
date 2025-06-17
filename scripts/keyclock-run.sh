#!/bin/bash

export KEYCLOAK_URL="https://auth.psdi-dev.vbox"
export KEYCLOAK_REALM="labtrove"
export KEYCLOAK_CLIENT_ID="data-conversion-service"
export KEYCLOAK_SECRET="WfHgzbu43Cppde48PbMWBLfL7kTovqxi"
export KEYCLOAK_REDIRECT_URL="https://dcs.psdi-dev.vbox/oidc_callback"

export SESSION_TIMEOUT_SECONDS="30"

export REQUESTS_CA_BUNDLE="/etc/pki/tls/certs/psdi/psdi-dev.crt"

export GUNICORN="/home/dgc/.local/bin/gunicorn"

$GUNICORN -b 0.0.0.0:5000 psdi_data_conversion.app:app --timeout 90