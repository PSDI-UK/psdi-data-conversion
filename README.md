# PSDI Data Conversion

This is the repository for the PSDI PF2 Chemistry File Format Conversion project. The goal of this project is to provide utilities to assist in converting files between the many different file formats used in chemistry, providing information on what converters are available for a given conversion and the expected quality of it, and providing multiple interfaces to perform these conversions. These interfaces are:

- Online web service, available at https://psdidev2.azurewebsites.net
- Version of the web app you can download and run locally (e.g. if you need to convert files which exceed the online app's file size limit)
- Command-line interface, to run conversions from a terminal

## Table of Contents

- [Project Structure](#project-structure)
- [Requirements](#requirements)
  - [Python](#python)
  - [Other Dependencies](#other-dependencies)
- [Command-Line Interface](#command-line-interface)
  - [Installation](#installation)
  - [Execution](#execution)
    - [Data Conversion](#data-conversion)
    - [Requesting Information on Possible Conversions](#requesting-information-on-possible-conversions)
- [Using the Online Conversion Service](#using-the-online-conversion-service)
- [Running the Python/Flask app locally](#running-the-pythonflask-app-locally)
  - [Installation and Setup](#installation-and-setup)
  - [Running the App](#running-the-app)
- [Testing](#testing)

## Project Structure

- `psdi_data_conversion`
  - `bin`
    - (Precompiled binaries for running file format converters)
  - `static`
    - `content`
      - (HTML assets for the web app)
    - `downloads` (created by app.py if not extant)
    - `img`
      - (image assets for the web app)
    - `javascript`
      - (JavaScript code for the web app)
    - `styles`
      - (CSS stylesheets for the web app)
    - `uploads` (created by app.py if not extant)
  - `templates`
    - (HTML assets rendered by Flask for the web app)
  - `__init.py__`
  - (Python modules and scripts)
- `scripts`
  - (Scripts used for project maintenance)
- `tests`
  - (Unit tests for the project)
- `CONTRIBUTING.md` (Guidelines for contributors to the project)
- `LICENSE` (Apache Licence version 2.0)
- `pyproject.toml` (Python project metadata and settings)
- `README.md` (This file)
- `requirements.txt` (Requirements for the Azure web app deployment of this project)
- `run_local.sh` (Helper script to run the web app locally)
- `startup.sh` (Startup script for the Azure web app)

## Requirements

### Python

Any local installation of this project requires Python 3.10 or greater. The best way to do this is dependant on your system, and you are likely to find the best tailored instructions by searching the web for e.g. "install Python 3.10 <your-os-or-distribution>". Some standard options are:

For Windows and MacOS: Download and run the installer for the latest version from the official site: https://www.python.org/downloads/

For Linux systems, Python is most readily installed with your distributions package manager. For Ubuntu/Debian-based systems, this is `apt`, and the following series of commands can be used to install the latest version of Python compatible with your system:

```bash
sudo apt update # Make sure the package manager has access to the latest versions of all packages
sudo apt upgrade # Update all installed packages
sudo apt install python3 # Install the latest possible version of Python
```

Check the version of Python installed with one of the following:

```bash
python --version
python3 --version
```

Usually `python` will be set up as an alias to python3, but if you already have an older version installed on your system, this might not be the case.

Also check that this installed Python's package manager, `pip`, on your system:

```bash
pip --version
```

If it didn't, you can manually install it with:

```bash
sudo apt install python3-pip
```

If this doesn't work, or the version installed is too low, an alternative is to install Python via the Anaconda package manager. For this, see the guide here: https://www.askpython.com/python/examples/install-python-with-conda

### Other Dependencies

This project depends on other projects available via pip, which will be installed automatically as required:

Required for all installations:

- `py`
- `openbabel-wheel`

Required to run the web app locally for a GUI experience:

- `Flask`
- `requests`

Required to run unit tests:

- `pytest`

In addition to the dependencies listed above, this project uses the assets made public by PSDI's common style project at https://github.com/PSDI-UK/psdi-common-style. The latest versions of these assets are copied to this project periodically (using the scripts in the `scripts` directory). In case a future release of these assets causes a breaking change in this project, the file `fetch-common-style.conf` can be modified to set a previous fixed version to download and use until this project is updated to work with the latest version of the assets.

## Command-Line Interface

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

#### Requesting Information on Possible Conversions

The script can also be used to get information on possible conversions by providing the `-l/--list` argument:

```bash
psdi-data-convert -l
```

Without any further arguments, the script will list converters available for use.

Further functionality planned for this script, but yet to be implemented:

- If the name of a converter is provided as an argument, it should provide information on the converter, such as what flags it will accept
- If the names of two formats are provided as arguments, it should provide information on the possible converters that can be used for this conversion and the expected quality of the conversion

## Using the Online Conversion Service

Enter https://psdidev2.azurewebsites.net in a browser. Guidance on usage is given on each page of the website.

## Running the Python/Flask app locally

### Installation and Setup

Install the package and its requirements, including the optional requirements used to run the GUI locally, by executing the following command from this project's directory:

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

If desired, you can modify the environmental variables set in this script to modify the operation - see the comments on each for details.

### Running the App

Run the `run_local.sh` script to start the server. You can then access the website by going to <http://127.0.0.1:5000> in a browser (this will also be printed in the terminal, and you can CTRL+click it there to open it in your default browser). Guidance for using the app is given on each page of it.

In case of problems when using Chrome, try opening Chrome from the command line:
open -a "Google Chrome.app" --args --allow-file-access-from-files

## Testing

Install the package requirements locally (ideally within a virtual environment) and test with pytest by executing the following commands from this project's directory:

```bash
source .venv/bin/activate # Create a venv first if necessary with `python -m venv .venv`
pip install .[gui,test]
pytest
```

## Contributors

- Ray Whorley
- Don Cruickshank
- Samantha Pearman-Kanza (s.pearman-kanza@soton.ac.uk)
- Bryan Gillis (7204836+brgillis@users.noreply.github.com)
