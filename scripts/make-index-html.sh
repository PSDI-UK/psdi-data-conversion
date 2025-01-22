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

if [ -z $TARGET_BASE_DIR ]; then
  TARGET_BASE_DIR=$DEFAULT_TARGET_BASE_DIR
fi
if [ -z $TARGET_HTML_DIR ]; then
  TARGET_HTML_DIR=$TARGET_BASE_DIR/html
fi

INDEX_HTML_DIR=$TARGET_HTML_DIR/index-versions

# The address that clicking on the PSDI brand on the left side of the header should link to. Default: "./", which will
# most likely be the home page of the project
export BRAND_LINK="."

# The location of a file containing HTML links to be added to the header, which will appear on the right side of it.
# Default: (no links will be added)
export HEADER_LINKS_SOURCE="$INDEX_HTML_DIR/header-links.html"

# The address of a directory containing the images referenced by the header and footer. By default, this will link
# to the live versions at https://psdi-uk.github.io/psdi-common-style/img/. This can be set to e.g. "./img" to instead
# link to local copies of these images
export IMG_LOC="$TARGET_BASE_DIR/img"

# Make new copies 
TARGET_DIR=$INDEX_HTML_DIR $SCRIPTS/copy_html.sh