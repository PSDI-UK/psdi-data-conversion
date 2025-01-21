#!/bin/bash

ROOTDIR=$(dirname -- $(readlink -f $BASH_SOURCE))

# Source the configuration file to get settings in its envvars, if it exists
CONF_FILE=$ROOTDIR/fetch-common-style.conf
if [ -f $CONF_FILE ]; then
    echo "Sourcing configuration from $CONF_FILE"
    source $CONF_FILE
else
    echo "No configuration file found at $CONF_FILE; configuration will be controlled by environmental variables"
fi

# Files to make copies of will be found in the HTML directory
ORIG_HEADER_FILENAME=$TARGET_HTML_DIR/psdi-common-header.html
ALT_HEADER_FILENAME=$TARGET_HTML_DIR/psdi-common-header-index.html
cp $ORIG_HEADER_FILENAME $ALT_HEADER_FILENAME

ORIG_FOOTER_FILENAME=$TARGET_HTML_DIR/psdi-common-footer.html
ALT_FOOTER_FILENAME=$TARGET_HTML_DIR/psdi-common-footer-index.html
cp $ORIG_FOOTER_FILENAME $ALT_FOOTER_FILENAME

# Update paths in both header and footer
for FILENAME in $ALT_HEADER_FILENAME $ALT_FOOTER_FILENAME; do
    
done