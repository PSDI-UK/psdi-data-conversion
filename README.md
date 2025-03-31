# PSDI Data Conversion

Version: Pre-release 2024-02-27

This is the repository for the PSDI PF2 Chemistry File Format Conversion project. The goal of this project is to provide utilities to assist in converting files between the many different file formats used in chemistry, providing information on what converters are available for a given conversion and the expected quality of it, and providing multiple interfaces to perform these conversions. These interfaces are:

- Online web service, available at https://data-conversion.psdi.ac.uk/
- Version of the web app you can download and run locally (e.g. if you need to convert files which exceed the online app's file size limit)
- Command-line application, to run conversions from a terminal
- Python library

## Table of Contents

- [Project Structure](#project-structure)
- [Requirements](#requirements)
  - [Python](#python)
  - [Other Dependencies](#other-dependencies)
- [Command-Line Application](#command-line-application)
  - [Installation](#installation)
  - [Execution](#execution)
    - [Data Conversion](#data-conversion)
    - [Requesting Information on Possible Conversions](#requesting-information-on-possible-conversions)
- [Python Library](#python-library)
  - [Installation](#installation-1)
  - [Use](#use)
    - [`run_converter`](#run_converter)
    - [`get_converter`](#get_converter)
    - [`constants`](#constants)
    - [`database`](#database)
  - [Further Information](#further-information)
- [Using the Online Conversion Service](#using-the-online-conversion-service)
- [Running the Python/Flask app locally](#running-the-pythonflask-app-locally)
  - [Installation and Setup](#installation-and-setup)
  - [Running the App](#running-the-app)
- [Testing](#testing)
- [Licencing](#licencing)
- [Contributors](#contributors)
- [Funding](#funding)

## Project Structure

- `.github`
  - `workflows`
    - (Automated workflows for various tasks related to project maintenance)
- `deploy`
  - (Files used as part of the deployment to STFC infrastructure)
- `psdi_data_conversion` (Primary source directory)
  - `bin`
    - (Precompiled binaries for running file format converters)
  - `static` (Static code and assets for the web app)
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
  - (Python packages, modules, and scripts)
- `scripts`
  - (Scripts used for project maintenance)
- `test_data`
  - (Files used for testing the project)
- `tests`
  - (Unit tests for the project)
- `CHANGELOG.md` (Updates since initial public release)
- `CONTRIBUTING.md` (Guidelines and information for contributors to the project)
- `DOCKERFILE` (Dockerfile for image containerising PSDI's data conversion service)
- `LICENSE` (Apache Licence version 2.0)
- `pyproject.toml` (Python project metadata and settings)
- `README.md` (This file)
- `requirements.txt` (Requirements for the web app deployment of this project)
- `run_local.sh` (Helper script to run the web app locally)

## Requirements

### Python

Any local installation of this project requires Python 3.12 or greater. The best way to do this is dependant on your system, and you are likely to find the best tailored instructions by searching the web for e.g. "install Python 3.12 <your-os-or-distribution>". Some standard options are:

For Windows and MacOS: Download and run the installer for the latest version from the official site: https://www.python.org/downloads/

For Linux systems, Python is most readily installed with your distribution's package manager. For Ubuntu/Debian-based systems, this is `apt`, and the following series of commands can be used to install the latest version of Python compatible with your system:

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

Usually `python` will be set up as an alias to python3, but if you already have an older version installed on your system, this might not be the case. You may be able to set this behaviour up by installing the `python-is-python3` package:

```bash
sudo apt install python-is-python3
```

Also check that this process installed Python's package manager, `pip`, on your system:

```bash
pip --version
```

If it didn't, you can manually install it with:

```bash
sudo apt install python3-pip
```

If this doesn't work, or the version installed is too low, an alternative is to install Python via the Anaconda package manager. For this, see the guide here: https://www.askpython.com/python/examples/install-python-with-conda. If you already have an earlier version of Python installed with Anaconda, you can install and activate a newer version with a command such as:

```bash
conda create --name converter python=3.12 anaconda # Where 'converter' is a possible conda environment name
conda activate converter
```

You can also install a newer version of Python if you wish by substituting "3.12" in the above with e.g. "3.13".

### Other Dependencies

This project depends on other projects available via pip, which will be installed automatically as required:

Required for all installations (`pip install .`):

- `py`
- `openbabel-wheel`

Required to run the web app locally for a GUI experience (`pip install .[gui]`):

- `Flask`
- `requests`

Required to run unit tests (`pip install .[test]`):

- `pytest`
- `coverage`

Required to run unit tests on the web app (`pip install .[gui-test]`):

- (all web app and test requirements listed above)
- `selenium`
- `webdriver_manager`

In addition to the dependencies listed above, this project uses the assets made public by PSDI's common style project at https://github.com/PSDI-UK/psdi-common-style. The latest versions of these assets are copied to this project periodically (using the scripts in the `scripts` directory). In case a future release of these assets causes a breaking change in this project, the file `fetch-common-style.conf` can be modified to set a previous fixed version to download and use until this project is updated to work with the latest version of the assets.

## Command-Line Application

### Installation

The CLA and Python library are installed together. This package is not yet available on PyPI, and so must be installed locally. This can be done most easily with:

```bash
pip install .
```

executed from this project's directory. You can also replace the '.' in this command with the path to this project's directory to install it from elsewhere.

Depending on your system, it may not be possible to install packages in this manner without creating a virtual environment to do so in. You can do this by first installing the `venv` module for Python3 with e.g.:

```bash
sudo apt install python3-venv # Or equivalent for your distribution
```

You can then create and activate a virtual environment with:

```bash
python -m venv .venv # ".venv" here can be replaced with any name you desire for the virtual environment
source .venv/bin/activate
```

You should then be able to install this project. When you wish to deactivate the virtual environment, you can do so with the `deactivate` command.

### Execution

Once installed, the command-line script `psdi-data-convert` will be made available, which can be called to either perform a data conversion or to get information about possible conversions and converters. You can see the full options for it by calling:

```bash
psdi-data-convert -h
```

This script has two modes of execution: Data conversion, and requesting information on possible conversions and converters.

#### Data Conversion

Data conversion is the default mode of the script. At its most basic, the syntax for it will look like:

```bash
psdi-data-convert filename.ext1 -t ext2
```

This will convert the file 'filename.ext1' to format 'ext2' using the default converter (Open Babel). A list of files can also be provided, and they will each be converted in turn.

The full possible syntax for the script is:

```
psdi-data-convert <input file 1> [<input file 2> <input file 3> ...] -t/--to <output format> [-f/--from <input file
format>] [-i/--in <input file location>] [-o/--out <location for output files>] [-w/--with <converter>] [--delete-input]
[--from-flags '<flags to be provided to the converter for reading input>'] [--to-flags '<flags to be provided to the
converter for writing output>'] [--from-options '<options to be provided to the converter for reading input>']
[--to-options '<options to be provided to the converter for writing output>'] [--coord-gen <coordinate generation
options] [-s/--strict] [--nc/--no-check] [-q/--quiet] [-g/--log-file <log file name] [--log-level <level>] [--log-mode
<mode>]
```

Call `psdi-data-convert -h` for details on each of these options.

#### Requesting Information on Possible Conversions

The script can also be used to get information on possible conversions by providing the `-l/--list` argument:

```bash
psdi-data-convert -l
```

Without any further arguments, the script will list converters available for use and file formats supported by at least one converter. More detailed information about a specific converter or conversion can be obtained through providing more information.

To get more information about a converter, call:

```
psdi-data-convert -l <converter name>
```

This will print general information on this converter, including what flags and options it accepts for all conversions, plus a table of what file formats it can handle for input and output.

To get information about which converters can handle a given conversion, call:

```
psdi-data-convert -l -f <input format> -t <output format>
```

This will provide a list of converters which can handle this conversion, and notes on the degree of success for each.

To get information on input/output flags and options a converter supports for given input/output file formats, call:

```
psdi-data-convert -l <converter name> [-f <input format>] [-t <output format>]
```

If an input format is provided, information on input flags and options accepted by the converter for this format will be provided, and similar for if an output format is provided.

## Python Library

### Installation

The CLA and Python library are installed together. See the [above instructions for installing the CLA](#installation), which will also install the Python library.

### Use

Once installed, this project's library can be imported through the following within Python:

```python
import psdi_data_conversion
```

The most useful modules and functions within this package to know about are:

- `psdi_data_conversion`
  - `converter`
    - `run_converter`
  - `constants`
  - `database`

#### `run_converter`

This is the standard method to run a file conversion. This method may be imported via:

```python
from psdi_data_conversion.converter import run_converter
```

For a simple conversion, this can be used via:

```python
run_converter(filename, to_format, name=name, data=data)
```

Where `filename` is the name of the file to convert (either fully-qualified or relative to the current directory), `to_format` is the desired format to convert to (e.g. `"pdb"`), `name` is the name of the converter to use (default "Open Babel"), and `data` is a dict of any extra information required by the specific converter being used, such as flags for how to read/write input/output files (default empty dict).

See the method's documentation via `help(run_converter)` after importing it for further details on usage.

#### `constants`

This package defines most constants used in the package. It may be imported via:

```python
from psdi_data_conversion import constants
```

Of the constants not defined in this package, the most notable are the names of available converters. Each converter has its own name defined in its module within the `psdi_data_conversion.converters` package (e.g. `psdi_data_conversion.converters.atomsk.CONVERTER_ATO`), and these are compiled within the `psdi_data_conversion.converter` module into:

- `D_SUPPORTED_CONVERTERS` - A dict which relates the names of all converters supported by this package to their classes
- `D_REGISTERED_CONVERTERS` - As above, but limited to those converters which can be run on the current machine (e.g. a converter may require a precompiled binary which is only available for certain platforms, and hence it will be in the "supported" dict but not the "registered" dict)
- `L_SUPPORTED_CONVERTERS`/`L_REGISTERED_CONVERTERS` - Lists of the names of supported/registered converters

#### `database`

The `database` module provides classes and methods to interface with the database of converters, file formats, and known possible conversions. This database is distributed with the project at `psdi_data_conversion/static/data/data.json`, but isn't user-friendly to read. The methods provided in this module provide a more user-friendly way to make common queries from the database:

- `get_converter_info` - This method takes the name of a converter and returns an object containing the general information about it stored in the database (note that this doesn't include file formats it can handle - use the `get_possible_formats` method for that)
- `get_format_info` - This method takes the name of a file format (its extension) and returns an object containing the general information about it stored in the database
- `get_degree_of_success` - This method takes the name of a converter, the name of an input file format (its extension), and the name of an output file format, and provides the degree of success for this conversion (`None` if not possible, otherwise a string describing it)
- `get_possible_converters` - This method takes the names of an input and output file format, and returns a list of converters which can perform the desired conversion and their degree of success
- `get_possible_formats` - This method takes the name of a converter and returns a list of input formats it can accept and a list of output formats it can produce. While it's usually a safe bet that a converter can handle any combination between these lists, it's best to make sure that it can with the `get_degree_of_success` method
- `get_in_format_args` and `get_out_format_args` - These methods take the name of a converter and the name of an input/output file format, and return a list of info on flags accepted by the converter when using this format for input/output
- `get_conversion_quality` - Provides information on the quality of a conversion from one format to another with a given converter. If conversion isn't possible, returns `None`. Otherwise returns a short string describing the quality of the conversion, a string providing information on possible issues with the conversion, and a dict providing details on property support between the input and output formats

### Further Information

The code documentation for the Python library is published online at https://psdi-uk.github.io/psdi-data-conversion/. Information on modules, classes, and methods in the package can also be obtained through standard Python methods such as `help()` and `dir()`.

## Using the Online Conversion Service

Enter https://data-conversion.psdi.ac.uk/ in a browser. Guidance on usage is given on each page of the website.

## Running the Python/Flask app locally

### Installation and Setup

Install the package and its requirements, including the optional requirements used to run the GUI locally, by executing the following command from this project's directory:

```bash
pip install .'[gui]'
```

If your system does not allow installs in this manner, it may be necessary to set up a virtual environment. See the instructions in the [command-line application installation](#installation) section above for how to do that, and then try to install again once you've set one up and activated it.

If you've cloned this repository, you can use the `run_local.sh` bash script to run the application. Otherwise (e.g. if you've installed from a wheel or PyPI), copy and paste the following into a script:

```bash
#!/bin/bash

# The envvar MAX_FILESIZE can be used to set the maximum allowed filesize in MB - 0 indicates no maximum
if [ -z $MAX_FILESIZE ]; then
  export MAX_FILESIZE=0
fi

# Logging control - "full" sets server-style logging, which is necessary to produce the output logs with the expected
# names. This should not be changed, or else errors will occur
export LOG_MODE=full

# The level to log at. Leave blank for defaults, which logs at INFO level for user output and ERROR level for the server
# log and stdout. If set to a different value here (e.g. DEBUG), all these channels will be set to that level
export LOG_LEVEL=

# The envvar SERVICE_MODE can be set to "true" to make this behave as if it's running as the public web service -
# uncomment the following line to enable that
# export SERVICE_MODE=true

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

## Extending Functionality

The Python library and CLA are written to make it easy to extend the functionality of this package to use other file format converters. This can be done by downloading or cloning the project's source from it's GitHub Repository (https://github.com/PSDI-UK/psdi-data-conversion), editing the code to add your converter following the guidance in the "[Adding File Format Converters](https://github.com/PSDI-UK/psdi-data-conversion/blob/main/CONTRIBUTING.md#adding-file-format-converters)" section of CONTRIBUTING.md to integrate it with the Python code, and installing the modified package on your system via:

```bash
pip install --editable .'[test]'
```

(This command uses the `--editable` option and optional `test` dependencies to ease the process of testing and debugging your changes.)

Note that when adding a converter in this manner, information on its possible conversions will not be added to the database, and so these will not show up when you run the CLA with the `-l/--list` option. You will also need to add the `--nc/--no-check` option when running conversions to skip the database check that the conversion is allowed.

## Testing

To test the CLA and Python library, install the optional testing requirements locally (ideally within a virtual environment) and test with pytest by executing the following commands from this project's directory:

```bash
pip install .'[test]'
pytest
```

To test the local version of the web app, install the GUI testing requirements locally (which also include the standard GUI requirements and standard testing requirements), start the server, and test by executing the GUI test script:

```bash
pip install .'[gui-test]'
./run_local.sh & # Start the server for the web app in the background
cd tests/selenium
./run.sh
kill %1 # Stop the web server - it may have a different job ID. If you don't know the job ID, you can alternatively call "fg" to bring the job to the foreground, then type CTRL+c to stop it
```

## Licencing

This project is provided under the Apache License version 2.0, the terms of which can be found in the file `LICENSE`.

This project redistributes compiled binaries for the Atomsk and c2x converters. These are both licenced under the
GNU General Public License version 3 and are redistributed per its terms. Any further redistribution of these binaries,
including redistribution of this project as a whole, including them, must also follow the terms of this licence.

This requires conspicuously displaying notice of this licence, providing the text of of the licence (provided here in
the files `psdi_data_conversion/bin/LICENSE_C2X` and `psdi_data_conversion/bin/LICENSE_ATOMSK`), and appropriately
conveying the source code for each of these. Their respective source code may be found at:

- Atomsk: https://github.com/pierrehirel/atomsk/
- c2x: https://www.c2x.org.uk/downloads/

## Contributors

- Ray Whorley
- Don Cruickshank
- Samantha Pearman-Kanza (s.pearman-kanza@soton.ac.uk)
- Bryan Gillis (7204836+brgillis@users.noreply.github.com)
- Tom Underwood

## Funding

PSDI acknowledges the funding support by the EPSRC grants EP/X032701/1, EP/X032663/1 and EP/W032252/1
