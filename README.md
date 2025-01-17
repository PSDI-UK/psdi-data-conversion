# CHEMISTRY FILE FORMAT CONVERSION DATABASE

This is the repository for the Pathfinder 2 Chemistry File Format Conversion source code.

## Directory structure of the Python Flask app website:

```
LICENSE
psdi_data_conversion
    __init__.py
    app.py
    atomsk.sh
    bin
        atomsk
        c2x
    converter.py
    log_utility.py
    main.py
    static
        content
            accessibility.htm
            convert.htm
            convertato.htm
            documentation.htm
            feedback.htm
            header-links.htm
            index-header-links.htm
            report.htm
        downloads (created by app.py if not extant)
        javascript
            convert.js
            convert.py
            convertato.js
            format.js
            load_accessibility.js
            report.js
        styles
            format.css
        uploads (created by app.py if not extant)
    templates
        index.htm
pyproject.toml
requirements.txt
  (azure-functions, requests, openbabel-wheel:
   installed by yaml workflow file on GitHub)
startup.sh
tests
    cli_test.py
    conftest.py
    converter_test.py
    logging_test.py
```

## Running the Python command-line interface

### Installation

This package is not yet available on PyPI, and so must be installed locally. This can be done most easily with:

```bash
pip install .
```

executed from this project's directory. You can also replace the '.' in this command with the path to this project's directory to install it from elsewhere.

### Execution

Once installed, the command-line script `psdi-data-convert` will be made available, which can be called to either perform a data conversion or to get information about possible conversions or converters (the latter TODO). You can see the full options for it by calling:

```bash
psdi-data-convert -h
```

This script has two modes of execution: Data conversion, and requesting information on possible conversions.

#### Data Conversion

Data conversion is the default mode of the script. At its most basic, the syntax for it will look like:

```bash
psdi-data-convert filename.ext1 -t ext2
```

This will convert the file 'filename.ext1' to format 'ext2' using the default converter (Open Babel). A list of files can also be provided, and they will each be converted in turn.

The full possible syntax for the script is:

```
psdi-data-convert <input file 1> [<input file 2>, <input file 3>, ...] -t/--to <output format> [-f/--from <input file format>] [-i/--in <input file location>] [-a/--at <location for output files>] [-w/--with <converter>] [-d] [--from-flags '<flags to be provided to the converter for reading input>'] [--to-flags '<flags to be provided to the converter for writing output>'] [--coord-gen <coordinate generation options] [-q/--quiet] [-l/--log-file <log file name] [--log-level <level>]
```

Call `psdi-data-convert -h` for details on each of these options.

#### Requesting information on possible conversions

The script can also be used to get information on possible conversions by providing the `-l/--list` argument:

```bash
psdi-data-convert -l
```

Without any further arguments, the script will list converters available for use.

Further functionality planned for this script, but yet to be implemented:

- If the name of a converter is provided as an argument, it should provide information on the converter, such as what flags it will accept
- If the names of two formats are provided as arguments, it should provide information on the possible converters that can be used for this conversion and the expected quality of the conversion

## Running the Python Flask app hosted on the Microsoft Azure site

Enter https://psdidev2.azurewebsites.net in a browser.

## Testing

Install the package requirements locally (ideally within a virtual environment) and test with pytest by executing the following commands from this project's directory:

```bash
source .venv/bin/activate # Create a venv first if necessary with `python -m venv .venv`
pip install .[test]
pytest
```

## Running the Python/Flask app locally

Python and Open Babel must be installed.

Install the package and its requirements by executing the following command from this project's directory:

```bash
pip install .
```

To enable debug mode, if required, enter:

```bash
export FLASK_ENV=development
```

If you've cloned this repository, you can then execute the `run_local.sh` bash script to run the application. Otherwise (e.g. if you've installed from a wheel or PyPI), copy and paste the following into a script and then execute it at the command-line:

```bash
#!/bin/bash

# Uncomment the following line to enable debug mode
# export FLASK_ENV=development

# Execute a local run of the application from the proper path

PACKAGE_PATH=`python -c "import psdi_data_conversion; print(psdi_data_conversion.__path__[0])"`
cd $PACKAGE_PATH/..
python -m flask --app psdi_data_conversion/app.py run
```

This will start the server. You can then access the website by going to <http://127.0.0.1:5000> in a browser (this will also be printed in the terminal, and you can CTRL+click it there to open it in your default browser).

The database can only be accessed from the University of Southampton or when using Global Connect. The current IP address must be added to the database's firewall rules on the Azure site.

In case of problems when using Chrome, try opening Chrome from the command line:
open -a "Google Chrome.app" --args --allow-file-access-from-files

## Using the website

Guidance on usage is given on each page of the website.

## Dependencies

In addition to the dependencies listed in the `pyproject.toml` file, this project depends on the assets made public by PSDI's common style project at https://github.com/PSDI-UK/psdi-common-style.. Any changes to these assets will be reflected in this project's web pages, and this project should ideally be tested with any changes before they're made live. An issue with retrieving these assets will appear as the website appearing unstyled and missing its header and footer.

In case these assets become no longer available for some reason, the commit `f1908b3627addfe5072c1e2ad4a648203bd8dee7` can be used as a reference to restore local versions of them.
