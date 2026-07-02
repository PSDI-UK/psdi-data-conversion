# Changelog for PSDI Data Conversion

## v0.4.0

### **Breaking Changes**

- The IDs of all formats, converters, and arguments have been changed from indices to UUIDs. This is necessary to allow new converter plugins to be added without fearing ID collisions. This has the following direct and indirect impacts which may require changes in code using the Python library or CLI:
  - All IDs of formats, converts, and arguments have been changed. The `doc` folder contains three `.json` files which list what these changes were, providing dicts of the old IDs to the new IDs, so these can be referenced to convert any IDs used to the new UUIDs
  - Some attributes of the `DataConversionDatabase` class returned by the method `psdi_data_conversion.database.get_database` have been deprecated, since the change from indices to UUIDs makes them non-functional or misleading. These are:
    - `d_converter_info` -> Renamed to `d_converter_info_from_name` (since now there's also a dict from ID). Previous name is still functional for now, but will give a deprecation warning
    - `l_converter_info` -> Fully deprecated. Functionality now replaced by `d_converter_info_from_id` (to look up by UUID) and `l_unsorted_converter_info` (to get an unsorted list, similar to using `list(self.d_converter_info_from_id.values())`)
    - `d_format_info` -> Renamed to `d_format_info_from_name` (since now there's also a dict from ID). Previous name is still functional for now, but will give a deprecation warning
    - `l_format_info` -> Fully deprecated. Functionality now replaced by `d_format_info_from_id` (to look up by UUID) and `l_unsorted_format_info` (to get an unsorted list, similar to using `list(self.d_format_info_from_id.values())`)
  - Some attributes of the `ConverterInfo` class returned by various methods in `psdi_data_conversion.database` to get information on a converter have been deprecated, since the change from indices to UUIDs makes them non-functional. These are:
    - `l_in_flag_info` -> Fully deprecated. Functionality now replaced by `d_in_flag_info` (to look up by UUID) and `l_unsorted_in_flag_info` (to get an unsorted list, similar to using `list(self.d_in_flag_info.values())`)
    - `l_out_flag_info` -> Ditto, replaced by `d_out_flag_info` and `l_unsorted_out_flag_info`
    - `l_in_option_info` -> Ditto, replaced by `d_in_option_info` and `l_unsorted_in_option_info`
    - `l_out_option_info` -> Ditto, replaced by `d_out_option_info` and `l_unsorted_out_option_info`
  - Since the graphs of the `ConversionsTable` class (`graph`, `supported_graph`, and `registered_graph`) don't support UUIDs as vertex IDs, they internally use indices for these vertices, which are generated at runtime. To support converting between UUIDs and vertex indices, the following dicts have been added to this class:
    - `d_indices_from_uuids`
    - `d_uuids_from_indices`
- The `ConversionQualityInfo` object returned by the method `psdi_data_conversion.database.get_conversion_quality` will now have attributes `in_format` and `out_format` always be `FormatInfo` objects, rather than the types of these matching the types of the formats input to this method

### New and Changed Functionality

- The database method `get_converter_info` can now be called without a `name` argument, and will return a list of info on all converters
- The converter info provided by queries to the database methods now contains member variables `supported` and `registered` indicating the status of the converter:
  - Both `False`: The converter is known to exist, but we currently provide no support for it
  - `supported==True`, `registered==False`: This package has a plugin ready to support this converter, but it cannot currently be used due to e.g. a required binary being missing which must be provided by the user
  - Both `True`: The converter is supported and ready to use
- `ConverterInfo`, `FormatInfo`, `FlagInfo`, and `OptionInfo` now have a `uuid` property which provides the ID converted to a UUID class

### Bugfixes

- Fixed a CLI bug where if the user requested info about a conversion where one of the formats is ambiguous, they would be shown a stack trace instead of a helpful message
- Fixed a CLI bug where if requesting info on both a converter and format, the user would always be told that conversion to/from the format with this converter was not possible

### Documentation Changes

- Fixed incorrect output type in documentation for database methods `get_in_format_args` and `get_out_format_args`
- Improved the error output for the database method `get_converter_info` if the converter name is not recognised, providing a list of known names

## v0.3.23

### Bugfixes

- Fixed light-mode/dark-mode toggle button functionality, which was causing the site to break when toggling into dark mode

### Testing Changes

- Updated test data to use newest version of Open Babel

## v0.3.19

### Documentation Changes:

- Updating text on package documentation intro page to no longer state it's a WIP and explain the purpose of the package doc
- Misc. fixes to code samples in documentation

## v0.3.17

### Miscellaneous Changes:

- Lowered Python requirement to 3.11, since this allows all functionality a user doing a local install should need. Still using 3.12 for deployment, since it has better security for unpacking archive files

## v0.3.16

### Documentation Changes

- Update README to display quick links at top

## v0.3.13

### Miscellaneous Changes

- Display code repo in project details now that it will be public

## v0.3.11

### Bugfixes

- In deployment, set all logs to be output unbuffered so they should be properly captured by Kibana

### Miscellaneous Changes

- Updated to use new deployment workflow with private deployment repo connected to STFC runners

## v0.3.7

### Bugfixes

- Fixed bug where conversions weren't being logged in service mode
- Fixed improper logging of input/output options when Open Babel conversions are called through the GUI
- Fixed checks on validity of Open Babel input/output flags/options

### Miscellaneous Changes

- Fixed bug with previous certificates which was resulting them being flagged as invalid by stricter browsers

## v0.3.6

### Miscellaneous Changes

- Updated certificates for web deployment (no change to pypi deployment)

## v0.3.5

### Bugfixes

- Fixed functionality of format-selection page on mobile browsers

## v0.3.4

### Stylistic Changes

- Updated page titles to new standardised PSDI format

## v0.3.3

### Stylistic Changes

- Each page of the GUI now has a unique title

## v0.3.2

### Bugfixes

- Removed incorrect underline under "archives" on Conversion page when it doesn't have a tooltip

## v0.3.1

### New and Changed Functionality

- The light-mode/dark-mode toggle button and Accessibility settings now play more nicely together, with the latest used taking precedence
- Added login functionality via Keycloak using a link in the GUI header (only in service mode). When running in service mode, logged out users will now be unable to convert archives of files and will have a reduced file size limit. This is also accompanied by appropriate text updates in the guidance, including tooltips for the use of archives and maximum file size

### Bugfixes

- Fixed requirements for generating documentation

### Formatting and Refactoring Changes

- All pages are now served via templates rendered with Flask, using inheritance to reduce code duplication
- Service mode and Production mode are now handled by variables used by Flask, not CSS stylings

### Bugfixes

- On the Accessibility page of the GUI, the "Default" options now properly show the default stylings rather than the current
- Fixed potential error with building the Docker image locally if the uploads/downloads directories already exist
- Fixed bug in testing where elements weren't properly scrolled into view

## v0.2.4

### Bugfixes

- Fixed false warnings which said that every format option was unsupported

## v0.2.3

### New and Changed Functionality

- When listing formats supported by a given converter in the command-line application, the description of each format will also be shown in the table
- A warning will now be printed to stderr and logged if an unrecognised format flag or option is provided for conversion with Open Babel

### Bugfixes

- Fixed coordinate generation quality not being properly logged

### Documentation Changes

- Fixed help for the "--from-flags", "--from-options" etc. command-line options to properly describe how values should be provided for them
- Add note to README about how to submit feedback and missing formats/conversion
- Updated README discussion of format IDs and disambiguated names, and provided more information about how to get IDs when formats are listed or when an ambiguous conversion is requested

## v0.2.2

### Bugfixes

- Fixed bug where c2x and Atomsk converters would fail if the current working directory wasn't the base directory of the project

### Testing Changes

- Disabled automated MacOS testing, which started failing due to an update on GitHub's end, while we decide how to fix it

## v0.2.1

### Bugfixes

- Fixed bug where when a conversion pathway is requested which turns out to be impossible, an exception is thrown instead of `None` being returned
- The logging level in the production deployment will now properly be INFO, while it will be DEBUG in the dev deployment
- Fixed the label for formats supporting 3D coordinates, which was unintentionally a duplicate of the 2D label
- Fixed crash when requesting info on a conversion which is impossible even with chained conversions

### Documentation Changes

- Added file `doc/conversion_chaining.md`, which explains the thought process behind the algorithm we (intend to) use for finding the best chained conversion

## v0.2.0

### New and Changed Functionality

- Changed the keyword arguments `upload_dir` and `download_dir` to `input_dir` and `output_dir` respectively
- Formats can now be specified case-insensitively
- When requesting details on a format through the command-line interface, details will be provided on which molecular properties it supports (e.g. whether or not it supports connections information)
- Added function `database.get_conversion_pathway` which can be used to get possible conversion routes between formats a direct conversion isn't possible with any converter
- When requesting details on two formats through the command-line interface and a direct conversion between them is not possible, a possible chain conversion will now be recommended

### Bugfixes

- Fixed bug where the `input_dir` keyword argument for `run_converter` was being ignored
- Fixed bug where the local-mode-only text was incorrectly appearing on the report page in service mode

### Testing Changes

- Excluded GUI modules from the calculating unit test coverage which can't be measured by the tool
- Added automated test that the production deployment is working on a schedule and after deploying to it

### Documentation Changes

- The Documentation page of the GUI now shows the mode that's being run, the most recent tag, and the SHA of the most recent commit (if this isn't the latest tagged commit)
- Updated release procedure and checklist in `CONTRIBUTING.md` to reflect current procedure

### Formatting and Refactoring Changes

- Changed Documentation and Accessibility pages of the GUI to work as Flask templates
- Cleaned up Flask files to not be all in one module
- Changed the database functionality to store possible conversions as a graph instead of a table
- Dockerfile now builds from `pyproject.toml`, with the now-unused `requirements.txt` removed

### Stylistic Changes

- Reformatted pages of the GUI/web app to use a two panel display, with instructions for components in boxes alongside them

## v0.1.7

### New and Changed Functionality

- Version, SHA, and service/prod modes now always shown in new About section on the Documentation page

### Documentation Changes

- Added information about deployment to CONTRIBUTING.md

### Bugfixes

- Environmental variable indicating dev or production mode should now be properly set for the deployed service

## v0.1.6

### New and Changed Functionality

- SHA banner at the bottom of home page now preferentially shows the version, only showing the SHA if the current version doesn't match the last tag

### Bugfixes

- Fixed bug which was blocking deployment to production

## v0.1.0

Initial public release. Features included:

- Online server functionality
- Locally-hosted server
- Command-line interface
- Python library
