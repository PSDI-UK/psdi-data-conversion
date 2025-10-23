# Changelog for PSDI Data Conversion

## v0.3.13

### Miscellaneous Changes

- Uisplay code repo in project details now that it will be public

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
