# Converter Scripts

This folder contains Bash scripts which are called by converter plugins and themselves wrap calls to converter binaries. These scripts aren't limited to one script per converter; the calling of scripts is handled by the plugin class defined in `psdi_data_conversion/converters/<name>/converter.py`, and this class may rely on multiple scripts if it's useful to do so.

The existing scripts may be referenced as examples for writing a new script. In particular, attention should be paid to the use of the `DIST` and `BIN_PATH` environmental variables to determine which executable should be used. E.g., from the `c2x.sh` script:

```bash
DEFAULT_DIST=linux

# The ennvar DIST can be used to set the distribution, indicating the subdirectory to search in for the binary
if [ -z $DIST ]; then
  DIST=$DEFAULT_DIST
fi

# If the envvar BIN_PATH is set, we use that path directly. Otherwise we use the binary distributed with this package
if [ -z $BIN_PATH ]; then
  BIN_PATH=psdi_data_conversion/bin/$DIST/c2x
fi

$BIN_PATH $1 $2 $3 $4 $5
```

The arguments provided to the script are determined by the converter plugin class. If the ready-made `ScriptFileConverter` base class is used, it will by default provide arguments to the script, in order:

- Name of output format, prefixed with "--", e.g. "--pdb"
- Input filename (fully-qualified)
- Output filename (fully-qualified)
- Input format flags (concatenated)
- Output format flags (concatenated)
- Input format options (space-separated)
- Output format options (space-separated)

If this doesn't meet your needs, the class's methods can be overridden to customise it as desired before these are passed to the script.
