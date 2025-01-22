# CHEMISTRY FILE FORMAT CONVERSION DATABASE

This is the repository for the Pathfinder 2 Chemistry File Format Conversion source code.

## Directory structure of the Python Flask app website:

```
CONTRIBUTING.md
LICENSE
psdi_data_conversion
    __init__.py
    app.py
    atomsk.sh
    bin
        atomsk
    converter.py
    log_utility.py
    static
        content
            index-versions
                header-links.htm
                psdi-common-footer.html
                psdi-common-header.html
            accessibility.htm
            convert.htm
            convertato.htm
            documentation.htm
            feedback.htm
            header-links.htm
            psdi-common-footer.html
            psdi-common-header.html
            report.htm
        downloads (created by app.py if not extant)
        img
            (misc. image assets)
        javascript
            convert.js
            convert.py
            convertato.js
            format.js
            load_accessibility.js
            psdi-common.js
            report.js
        styles
            format.css
            psdi-common.css
        uploads (created by app.py if not extant)
    templates
        index.htm
pyproject.toml
requirements.txt
  (azure-functions, requests, openbabel-wheel:
   installed by yaml workflow file on GitHub)
startup.sh
scripts
    fetch-common-style.conf
    fetch-common-style.sh
    make-index-html.sh
tests
    logging_test.py
```

## Running the Python Flask app hosted on the Microsoft Azure site

Enter https://psdidev2.azurewebsites.net in a browser.

## Testing

Install the package requirements locally (ideally within a virtual environment) and test with pytest:

```bash
source .venv/bin/activate # Create a venv first if necessary with `python -m venv .venv`
pip install .[gui,test]
pytest
```

## Running the Python/Flask app locally

Python, Open Babel, and Flask must be installed.

Install the package and its requirements, including the optional requirements used to run the GUI locally via:

```bash
pip install .[gui]
```

To enable debug mode, if required, enter:

```bash
export FLASK_ENV=development
```

If you've cloned this repository, you can use the `run_local.sh` bash script to run the application. Otherwise (e.g. if you've installed from a wheel or PyPI), copy and paste the following into a script:

```bash
#!/bin/bash

# The envvar MAX_FILESIZE can be used to set the maximum allowed filesize in MB - 0 indicates no maximum
if [ -z $MAX_FILESIZE ]; then
  export MAX_FILESIZE=0
fi

# The envvar LOGGING can be used to set how logs are stored. Allowed values are:
# Full - multi-file logging, only recommended when running as a public web app
# Simple - logs saved to one file
# None - output only to stdout/stderr
if [ -z $LOGGING ]; then
  export LOGGING=Simple
fi

# Uncomment the following line to enable debug mode
# export FLASK_ENV=development

# Execute a local run of the application from the proper path

PACKAGE_PATH=`python -c "import psdi_data_conversion; print(psdi_data_conversion.__path__[0])"`
cd $PACKAGE_PATH/..
python -m flask --app psdi_data_conversion/app.py run
```

If desired, you can modify the environmental variables set in this script to modify the operation - see the comments on each for details. Running this script will start the server. You can then access the website by going to <http://127.0.0.1:5000> in a browser (this will also be printed in the terminal, and you can CTRL+click it there to open it in your default browser).

The database can only be accessed from the University of Southampton or when using Global Connect. The current IP address must be added to the database's firewall rules on the Azure site.

In case of problems when using Chrome, try opening Chrome from the command line:
open -a "Google Chrome.app" --args --allow-file-access-from-files

## Using the website

Guidance on usage is given on each page of the website.

## Dependencies

In addition to the dependencies listed in the `pyproject.toml` file, this project uses the assets made public by PSDI's common style project at https://github.com/PSDI-UK/psdi-common-style. The latest versions of these assets are copied to this project periodically (using the scripts in the `scripts` directory). In case a future release of these assets causes a breaking change in this project, the file `fetch-common-style.conf` can be modified to set a previous fixed version to download and use until this project is updated to work with the latest version of the assets.
