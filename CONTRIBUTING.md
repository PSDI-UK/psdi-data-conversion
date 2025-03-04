# Contributing

## Git Workflow

This project uses a version of [GitLab Flow](https://about.gitlab.com/topics/version-control/what-is-gitlab-flow/) for its workflow. This uses two primary long-lived branches, `main` and `release`, plus short-lived feature and release candidate branches:

- `main` - This is the default and main working branch of the project. Day-to-day changes start from this branch and are merged back into it when ready. `main` is protected so that it can only be merged to via Pull Request.
- `release` - This is the branch that is used for deployments of the project. It is periodically updated by a release candidate branch being created from `main`, tested, and merged into this when approved. `release` is protected so that it can only be merged to via Pull Request.
- Feature branches - These are used for day-to-day work. Feature branches are branched off of `main` and then worked on (usually by a single developer). They may be discarded if it's decided they are unneeded, or else if they're approved, they're merged back into `main` via a Pull Request. To avoid issues which can arise from feature branches being open too long and diverging from `main`, it is a good idea to aim for each to only represent a small, discrete change. A good target is approximately one day's worth of work, and only exceptionally should a feature branch be open for more than a week. The recommended naming convention for feature branches is `<initials>-<description-of-change>`, e.g. `brg-new-feature`.
- Release candidate branches - These branches are used to bridge `main` and `release`. When it's time to prepare for a release, a release candidate branch is branched off of `main` to isolate it from any further changes. It then undergoes full testing, including manual testing of the web app. If any issues are identified which require fixing, the branch can be updated until all tests pass. Once this is the case, the branch is merged into `release` as well as back into `main` (so that any fixes made in it are also reflected there). The recommended naming convention for release candidate branches is `rc-<target-version-number>`, e.g. `rc-0.1.4`.

## Release Checklist and Procedure

The following tasks should be completed before merging a release candidate branch to `release`:

- Open a draft pull request from this release branch to `release`, which indicates the target version number
- Ensure that all automated tests and checks pass - these should be run automatically on the PR opened above
- Manually test the local web interface
  - If there have been any changes to the Python backend, run a test that a file can be converted successfully and produces a proper log
  - If there have been any changes to the web frontend, check the appearance of the site to ensure that it looks as desired. Test the Accessibility page to ensure that changes there work properly, are saved when requested and apply to other pages
- Check that `CHANGELOG.md` is up-to-date with all changes in this version. Any subsections for categories with no changes in this version can be removed to keep the file concise
- Check that the project version is updated to the desired new version in all places it appears:
  - `CHANGELOG.md` (The top section should reflect the new version)
  - `pyproject.toml` (if/when the version number is added to it - currently it gets the version automatically from the last tag)

If any of these tasks fail and require changes, make the needed changes and then recheck that all other tasks still pass. E.g. if testing the local web interface reveals a bug in the Python backend that needs to be fixed, ensure that all automated tests still pass after doing so

Then, follow the following steps to make the release:

1. Merge the pull request to `release`. The release candidate branch can be safely deleted
2. Merge `release` into `main` via PR (obviously don't delete `release` - if it even gives you the option to, something has gone wrong in the project rulesets, so report this)
3. Create a new release for the project on GitHub from the `release` branch, and tag the latest commit to `release` in the process. The release and tag should both be named `v<version-number>`, e.g. `v0.1.4`. The description can be used to highlight any important changes.

## Changelog

Explanation of possible sections in the CHANGELOG:

- **Breaking Changes:** Any changes which would result in previous integration with this project breaking and no longer being compatible with the new version (e.g. a command-line option is removed)
- **Deprecated Features:** Any features which are planned to be removed or changed in a breaking way in a future release (e.g. a command-line option will be removed soon, and perhaps an alternative way to achieve the same goal already exists)
- **New and Changed Functionality:** Anything new the project can do for the user
- **Bugfixes:** Any issues fixed
- **Testing Changes:** New unit tests or significant changes to existing tests, not including any incidental updates to tests which are necessary for other features already noted. E.g. if a test is implemented to ensure that a bugfix noted in this section works properly, that doesn't need to be noted
- **Documentation Changes:** Any notable updates to `README.md`, `CONTRIBUTING.md`, docstrings, and comments. This doesn't need to include incidental changes - for instance, if a new command-line option is added, it isn't necessary to make a separate note that documentation for it is added. However, if documentation was previously missing for an option and is now added, that should be noted
- **Formatting and Refactoring Changes:** Any non-functional changes to the code (e.g. "Changed all python variables to `snake_case`", "Refactored set of separate functions to use the same common code")
- **Stylistic Changes:** Non-functional changes to the aesthetic appearance of the web app or the formatting of text displayed to the user by the CLI (including logs)
- **Miscellaneous Changes:** Anything that doesn't fit into one of the above categories, such as project meta changes (e.g. "Implemented new GitHub action to lint JavaScript code")

The below can be used as a template for new sections to be added to `CHANGELOG.md` with each new release. When a release is made, any sections without any entries can be removed to help keep the file concise:

### Breaking Changes

-

### Deprecated Features

-

### New and Changed Functionality

-

### Bugfixes

-

### Testing Changes

-

### Documentation Changes

-

### Formatting and Refactoring Changes

-

### Stylistic Changes

-

### Miscellaneous Changes

-

## Editing Advice

### Adding File Format Converters

If you wish to make a converter accessible by the CLI, you only need to follow the steps in the Python Integration section below. If you also wish to make it available in the web app, you will additionally have to follow this instructions in the Web App Integration section.

#### Python Integration

In the Python layer of the code, each file format converter is defined in its own module in the `psdi_data_conversion.converters` package, as a subclass of the `FileConverter` class defined in `psdi_data_conversion.converters.base`. A new converter can be integrated with the Python layer of the code by:

1. Create a new module in the `psdi_data_conversion.converters` for the converter
2. Define a new class for this converter, as a subclass of `FileConverter`
3. Set the class variable `name` for the converter to be the name you want to use for this converter. This will be what needs to be specified in the command-line to request this converter
4. Implement the `_convert(self)` method of the class with a call to the converter. This method must create the converted file at the location specified by the variable `self.out_filename` (which is provided fully qualified) and set the variables `self.out` and `self.err` with normal output and error output respectively (at minimum they must be set to empty strings)
5. After defining the converter's class, set the module-level variable `converter` to the class

This will look something like:

```python
from psdi_data_conversion.converters.base import FileConverter

CONVERTER_MY = 'My Converter'

class MyFileConverter(FileConverter):
    """File Converter specialized to use my method for conversions
    """

    name = CONVERTER_MY

    def _convert(self):

        # Run whatever steps are necessary to perform the conversion
        create_my_converted_file_at(self.out_filename)

        self.out = "Standard output goes here"
        self.err = "Errors go here"

# Assign this converter to the `converter` variable - this lets the psdi_data_conversion.converter module detect and
# register it, making it available for use by the CLI and web app
converter = MyFileConverter
```

That's all you need to do! The `psdi_data_conversion.converter` module parses all modules in the `converters` package to find converters, so if you've done everything correctly, it will find the new converter and register it for you. You can test that it is properly registered by using the CLI to run:

```bash
psdi-data-convert -l
```

Your new converter should appear, or else you will probably see an error message which will detail an exception raised when trying to register it.

For file converters which can be run with a call to a script, this can be streamlined even further by taking advantage of the `ScriptFileConverter` subclass. With this, the converter's subclass can be defined even more succinctly:

```python
from psdi_data_conversion.converters.base import ScriptFileConverter

CONVERTER_MY_SCRIPT = 'My Script Converter'

class MyScriptFileConverter(ScriptFileConverter):
    """File Converter specialized to use my script for conversions
    """

    name = CONVERTER_MY_SCRIPT
    script = "my_script.sh"

converter = MyScriptFileConverter
```

When a converter is defined this way, the `_convert(self)` method will be defined to execute a subprocess call to run the script defined in the class's `script` class variable, searching for it in the `psdi_data_conversion/scripts` directory. Typically this will be a wrapper to call a binary to perform the conversion, which should be placed in the `psdi_data_conversion/bin` directory. It will pass to it the fully-qualified input filename (`self.in_filename`) as the first argument, the fully-qualified output filename (`self.out_filename`) as the second argument, and then any flags defined in `self.data["to_flags"]` and `self.data["from_flags"]`.

Finally, it's good practice to add a unit test of the converter. You can do this by following the example of tests in `tests/converter_test.py`. If necessary, add a (small) file it can convert to the `test_data` folder, and implement a test that it can convert it to another format by adding a new method to the `TestConverter` class in this file. At its simplest, this method should look something like:

```python
    def test_c2x(self):
        """Run a test of the C2X converter on a straightforward `.pdb` to `.cif` conversion
        """

        self.get_input_info(filename="hemoglobin.pdb",
                            to="cif")

        # "from" is a reserved word so we can't set it as a kwarg in the function call above
        self.mock_form["from"] = "pdb"

        self.run_converter(name=CONVERTER_C2X)

        # Check that the input file has been deleted and the output file exists where we expect it to
        self.check_file_status(input_exist=False, output_exist=True)
```

Ensure that the method you add starts with `test_` so that it will be detected by `pytest`. The basic check here that the output file exists can be extended to check that the details of it are as expected.

It may also be useful to add a test that the converter fails when you expect it to. This can be done e.g. with a test method that looks like:

```python
    def test_xyz_to_inchi_err(self):
        """Run a test of the converter on an `.xyz` to `.inchi` conversion we expect to fail
        """

        self.get_input_info(filename="quartz_err.xyz",
                            to="inchi")

        # "from" is a reserved word so we can't set it as a kwarg in the function call above
        self.mock_form["from"] = "xyz"

        # Pass the `expect_code` argument to the call to run the converter. This causes it to check that when it runs,
        # the conversion process aborts with the provided error code
        self.run_converter(name=CONVERTER_OB,
                           expect_code=const.STATUS_CODE_GENERAL)

        # Check that the input and output files have properly been deleted
        self.check_file_status(input_exist=False, output_exist=False)
```

If the test is more complicated that this, you can implement a modified version of `self.run_converter` within the test method to perform the desired test. You can also check that output logs include the desired information by either opening the log filenames or using PyTest's `capsys` feature, which captures output to logs, stdout, and stderr.

You can then run the any tests you added, plus the existing test suite through running the following commands from the project's root directory:

```bash
source .venv/bin/activate # Create a venv first if necessary with `python -m venv .venv`
pip install --editable .'[test]'
pytest
```

This installs the project in a virtual environment in 'editable' mode (which means the source files will be used from where they are rather than being copied, so any changes to them will be directly reflected in tests and uses of the CLI) and then calls `pytest` to run the unit tests in the project. `pytest` will automatically pick up any extra tests you add and run them as well.

#### Web App Integration

TODO: Write guide for adding converter to the web app.

List of necessary steps:

- Update 'converters' table in database.
- Update 'formats' table in database.
- Update 'converts_to' table in database.
- Find/compile suitable Linux binary and upload it
- Write script to call binary
- New HTML file for conversion page.
- New associated JS file.

### Debugging

For debugging python issues, it's recommended to install the package in editable mode via pip. This sets it up so that the python source files are used in-place rather than copied to a separate install directory, meaning that changes to them will be reflected in runs without need for a new installation. This can be done through the following command (which also installs all optional packages):

```bash
pip install --editable .'[gui-test]'
```
