# Changelog for PSDI Data Conversion

## Since v0.1.6

### New and Changed Functionality

- Changed the keyword arguments `upload_dir` and `download_dir` to `input_dir` and `output_dir` respectively
- Formats can now be specified case-insensitively

### Bugfixes

- Fixed bug where the `input_dir` keyword argument for `run_converter` was being ignored

### Testing Changes

- Excluded GUI modules from the calculating unit test coverage which can't be measured by the tool

### Documentation Changes

- The Documentation page of the GUI now shows the mode that's being run, the most recent tag, and the SHA of the most recent commit (if this isn't the latest tagged commit)

### Formatting and Refactoring Changes

- Changed Documentation and Accessibility pages of the GUI to work as Flask templates
- Cleaned up Flask files to not be all in one module

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
