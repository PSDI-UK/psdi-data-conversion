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

# Set alternate variables for the data that will differ for the index version of files
INDEX_HTML_DIR=$TARGET_HTML_DIR/index-versions
export BRAND_LINK="."
export HEADER_LINKS_SOURCE="$INDEX_HTML_DIR/header-links.html"
export IMG_LOC="$TARGET_BASE_DIR/img"

# Find the location of the copy_html.sh script
PACKAGE_FILENAME=$ROOTDIR/psdi-assets.tar.gz
ASSET_SUBDIR=`tar tf $PACKAGE_FILENAME | head -n 1`
SCRIPTS=$ROOTDIR/$ASSET_SUBDIR/scripts

# Make new copies 
TARGET_DIR=$INDEX_HTML_DIR $SCRIPTS/copy_html.sh