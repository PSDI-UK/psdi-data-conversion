# Changelog for PSDI Data Conversion

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
