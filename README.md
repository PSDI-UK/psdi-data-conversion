# PSDI Data Conversion

[![License Badge](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
![Coverage Badge](https://img.shields.io/endpoint?url=https://gist.githubusercontent.com/brgillis/dbd938192dc4de9b7779978e515c0e79/raw/covbadge.json)

Release date: 2024-04-29

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
- [Troubleshooting](#troubleshooting)
  - [OSError: [Errno 24] Too many open files](#oserror-errno-24-too-many-open-files)
  - [Errors running c2x or Atomsk converters](#errors-running-c2x-or-atomsk-converters)
  - [A supported conversion fails](#a-supported-conversion-fails)
    - [Input file is malformatted or corrupt](#input-file-is-malformatted-or-corrupt)
    - [Input file's format is misidentified](#input-files-format-is-misidentified)
    - [Other known issues](#other-known-issues)
- [Licensing](#licensing)
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
  - `gui`
    - (Unit tests for the GUI, aka the local version of the web app)
  - `python`
    - (Unit tests for the Python library and command-line application)
- `CHANGELOG.md` (Updates since initial public release)
- `CONTRIBUTING.md` (Guidelines and information for contributors to the project)
- `DOCKERFILE` (Dockerfile for image containerising PSDI's data conversion service)
- `LICENSE` (Apache License version 2.0)
- `pyproject.toml` (Python project metadata and settings)
- `README.md` (This file)
- `requirements.txt` (Requirements for the web app deployment of this project)

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

The CLA and Python library are installed together. This project is available on PyPI, and so can be installed via pip with:

```bash
pip install psdi-data-conversion
```

If you wish to install from source, this can be done most easily by cloning the project and then executing:

```bash
pip install .
```

from this project's directory. You can also replace the '.' in this command with the path to this project's directory to install it from elsewhere.

**Note:** This project uses git to determine the version number. If you clone the repository, you won't have to do anything special here, but if you get the source e.g. by extracting a release archive, you'll have to do one additional step before running the command above. If you have git installed, simply run `git init` in the project directory and it will be able to install. Otherwise, edit the project's `pyproject.toml` file to uncomment the line that sets a fixed version, and comment out the lines that set it up to determine the version from git - these are pointed out in the comments there.

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

Note that some requested conversions may involve ambiguous formats which share the same extension. In this case, the application will print a warning and list possible matching formats, with a disambiguating name that can be used to specify which one. For instance, the `c2x` converter can convert into two variants of the `pdb` format, and if you ask it to convert to `pdb` without specifying which one, you'll see:

```
ERROR: Extension 'pdb' is ambiguous and must be defined by ID. Possible formats and their IDs are:
9: pdb-0 (Protein Data Bank)
259: pdb-1 (Protein Data Bank with atoms numbered)
```

This provides the IDs ("9" and "259") and disambiguating names ("pdb-0" and "pdb-1") for the matching formats. Either can be used in the call to the converter, e.g.:

```bash
psdi-data-conversion nacl.cif -t 9 -w c2x
# Or equivalently:
psdi-data-conversion nacl.cif -t pdb-0 -w c2x
```

The "<format>-0" pattern can be used with any format, even if it's unambiguous, and will be interpreted as the first instance of the format in the database with valid conversions. Note that as the database expands in future versions and more valid conversions are added, these disambiguated names may change, so it is recommended to use the format's ID in scripts and with the library to ensure consistency between versions of this package.

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

Note that as with running the application through the command-line, some extra care may be needed in the case that the input or output format is ambiguous - see the [Data Conversion](#data-conversion) section above for more details on this. As with running through the command-line, a format's ID or disambiguated name must be used in the case of ambiguity.

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

This project is available on PyPI, and so can be installed via pip, including the necessary dependencies for the GUI, with:

```bash
pip install psdi-data-conversion'[gui]'
```

If you wish to install the project locally from source, this can be done most easily by cloning the project and then executing:

```bash
pip install .'[gui]'
```

**Note:** This project uses git to determine the version number. If you clone the repository, you won't have to do anything special here, but if you get the source e.g. by extracting a release archive, you'll have to do one additional step before running the command above. If you have git installed, simply run `git init` in the project directory and it will be able to install. Otherwise, edit the project's `pyproject.toml` file to uncomment the line that sets a fixed version, and comment out the lines that set it up to determine the version from git - these are pointed out in the comments there.

If your system does not allow installation in this manner, it may be necessary to set up a virtual environment. See the instructions in the [command-line application installation](#installation) section above for how to do that, and then try to install again once you've set one up and activated it.

### Running the App

Once installed, the command-line script `psdi-data-convert-gui` will be made available, which can be called to start the server. You can then access the website by going to <http://127.0.0.1:5000> in a browser (this will also be printed in the terminal, and you can CTRL+click it there to open it in your default browser). Guidance for using the app is given on each page of it. When you're finished with the app, key CTRL+C in the terminal where you called the script to shut down the server, or, if the process was backgrounded, kill the appropriate process.

In case of problems when using Chrome, try opening Chrome from the command line:
open -a "Google Chrome.app" --args --allow-file-access-from-files

The local version has some customisable options for running it, which can can be seen by running `psdi-data-convert-gui --help`. Most of these are only useful for development, but one notable setting is `--max-file-size-ob`, which sets the maximum allowed filesize for conversions with the Open Babel converter in megabytes. This is set to 1 MB by default, since Open Babel has a known bug which causes it to hang indefinitely for some conversions over this size (such as from large `mmcif` files). This can be set to a higher value (or to 0 to disable the limit) if the user wishes to disable this safeguard.

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
pytest tests/python
```

To test the local version of the web app, install the GUI testing requirements locally (which also include the standard GUI requirements and standard testing requirements), start the server, and test by executing the GUI test script:

```bash
pip install .'[gui-test]'
pytest tests/gui
```

Both of these sets of tests can also be run together if desired through:

```bash
pip install .'[gui-test]'
pytest
```

## Troubleshooting

This section presents solutions for commonly-encountered issues.

### OSError: [Errno 24] Too many open files

You may see the error:

```
OSError: [Errno 24] Too many open files
```

while running the command-line application, using the Python library, or running tests This error is caused by a program hitting the limit of the number of open filehandles allowed by the OS. This limit is typically set to 1024 on Linux systems and 256 on MacOS systems, and thus this issue occurs much more often on the latter. You can see what your current limit is by running the command:

```bash
ulimit -a | grep "open files"
```

This limit can be temporarily increased for the current terminal session by running the command:

```bash
ulimit -n 1024 # Or another, higher number
```

First, try increasing the limit and then redo the operation which caused this error to see if this resolves it. If this does, you can make this change permanent in a few ways, the easiest of which is to add this command to your `.bashrc` file so it will be set for every new terminal session.

If you see this error when the filehandle limit is already set to a high value such as 1024, this may indicate the presence of a bug in the project which causes a leaking of filehandles, so please open an issue about it, pasting the error you get and the details of your system, particularly including your current filehandle limit.

### Errors running c2x or Atomsk converters

We provide support for the c2x and Atomsk converters by packing binaries which support common Linux and MacOS platforms with this project, but we cannot guarantee universal support for these binaries. In particular, they may rely on dynamically-linked libraries which aren't installed on your system.

Look through the error message you received for messages such as "Library not loaded" or "no such file", and see if they point to the name of a library which you can try installing. For instance, if you see that it's searching for `libquadmath.0.dylib` but not finding it, you can try installing this library. In this case, this library can be installed through apt with:

```bash
sudo apt install libquadmath0
```

or through brew via:

```bash
brew install gcc
```

Alternatively, you can run your own versions of the `c2x` and `atomsk` binaries with this project. Compile them yourself however you wish - see the projects at https://github.com/codenrx/c2x and https://github.com/pierrehirel/atomsk and follow their instructions to build a binary on your system. Once you've done so, add the binary to your `$PATH`, and this project will pick that up and use it in preference to the prepackaged binary.

On the other hand, it's possible that an error of this sort will occur if you have a non-working binary of one of these converters in your `$PATH`. If this might be the case, you can try removing it and see if the prepackaged binary works for you, or else recompile it to try to fix errors.

### A supported conversion fails

Here we'll go over some of the most common reasons that a supported conversion might fail, and what can be done to fix the issue.

#### Input file is malformatted or corrupt

Usually if there is a problem with the input file, the error message you receive should indicate some difficulty reading it. If the error message indicates this might be the issue, try the following:

Check the validity of the input file, ideally using another tool which can read in a file of its type, and confirm that it can be read successfully. This doesn't guarantee that the file is valid, as some utilities are tolerant to some formatting errors, but if you get an error here, then you know the issue is with the file. If the file can be read by another utility, see if the conversion you're attempting is supported by another converter - it might be that the file has a negligible issue that another converter is able to ignore.

If you've confirmed that the input file is malformatted or corrupt, see if it's possible to regenerate it or fix it manually. There may be a bug in the program which generated it - if this is under your control, check the format's documentation to help fix it. Otherwise, you can see if you can use the format's documentation as a guide to manually fix the file.

#### Input file's format is misidentified

If you've followed the steps in the previous section and confirmed that the input file is valid, but you're still having issues with it, one possibility is that this application is misidentifying the file's format. This can happen if you've given the file an extension which isn't expected of its format, or in rare cases where an extension is shared by multiple formats.

To remedy this, try explicitly specifying the format, rather than letting the application guess it from the extension. You can see all supported formats by running `psdi-data-convert -l`, and then get details on one with `psdi-data-convert -l -f <format>` to confirm that it's the correct format. You can then call the conversion script with the argument `-f <format>`, or within Python make a call to the library with `run_converter(..., from_format=<format>)` to specify the format.

`<format>` here can be the standard extension of the format (in the case of unambiguous extensions), its ID, or its disambiguated name. To give an example which explains what each of these are, let's say you have an MDL MOL file you wish to convert to XYZ, so you get information about it and possible converters with `psdi-data-convert -l -f mol -t xyz`:

```base
$ psdi-data-convert -l -f mol
WARNING: Format 'mol' is ambiguous and could refer to multiple formats. It may be necessary to explicitly specify which
you want to use when calling this script, e.g. with '-f mol-0' - see the disambiguated names in the list below:

18: mol-0 (MDL MOL)
216: mol-1 (MOLDY)

20: xyz (XYZ cartesian coordinates)

The following registered converters can convert from mol-0 to xyz:

    Open Babel
    c2x

For details on input/output flags and options allowed by a converter for this conversion, call:
psdi-data-convert -l <converter name> -f mol-0 -t xyz

The following registered converters can convert from mol-1 to xyz:

    Atomsk

For details on input/output flags and options allowed by a converter for this conversion, call:
psdi-data-convert -l <converter name> -f mol-1 -t xyz
```

This output indicates that the application is aware of two formats which share the `mol` extension: MDL MOL and MOLDY. It lists the ID, disambiguated name, and description of each: ID `18` and disambiguated name `mol-0` for MDL MOL, and ID `216` and disambiguated name `mol-1` for MOLDY. The XYZ format, on the other hand, is unambiguous, and only lists the standard extension for it as its disambiguated name (although `xyz-0` will be accepted without error as well).

The program then lists converters which can handle the requested conversion, revealing a potential pitfall: The Open Babel and c2x converters can convert from MDL MOL to XYZ, which the Atomsk converter can convert from MOLDY to XYZ. If you don't specify which format you're converting from, the script might assume you meant to use the other one, if that's the only one compatible with the converter you've requested (or with the default converter, Open Babel, if you didn't explicitly request one). So to be careful here, it's best to specify this input format unambiguously.

Since in this example you have an MDL MOL file, you would use `-f 18` or `-f mol-0` to explicitly specify it in the command-line, or similarly provide one of these to the `from_format` argument of `run_converter` within Python. The application will then properly handle it, including alerting you if you request a conversion that isn't supported by your requested converter (e.g. if you request a conversion of this MDL MOL file to XYZ with Atomsk).

Important note: The disambiguated name is generated dynamically and isn't stored in the database, and in rare cases may change for some formats in future versions of this application which expand support to more formats and conversions. For uses which require forward-compatibility with future versions of this application, the ID should be used instead.

#### Other known issues

Through testing, we've identified some other conversion issues, which we list here:

- Open Babel will indefinitely hang when attempting to convert large files (more than ~1 MB) of certain types (such as `mmcif`). This is an issue with the converter itself and not our application, which we hope will be fixed in a future version. If this occurs, the job will have to be forcibly terminated. CTRL+C will fail to terminate it, but it can be stopped with CTRL+Z, then terminated with `kill %N`, where N is the number listed beside the job when it is stopped (usually 1). The conversion should then be attempted with another supported converter.

## Licensing

This project is provided under the Apache License version 2.0, the terms of which can be found in the file `LICENSE`.

This project redistributes compiled binaries for the Atomsk and c2x converters. These are both licensed under the
GNU General Public License version 3 and are redistributed per its terms. Any further redistribution of these binaries,
including redistribution of this project as a whole, including them, must also follow the terms of this license.

This requires conspicuously displaying notice of this license, providing the text of of the license (provided here in
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
