#!/usr/bin/env python3

"""psdi_data_conversion/main_create_plugin.py
=============

Created 2026-07-01 by Bryan Gillis.

Entry-point file for the script to create a new converter plugin.
"""

import importlib
import os
import re
import shutil
import sys
from argparse import ArgumentParser

from psdi_data_conversion.converters import base as converters_base
from psdi_data_conversion.utils import TextColors

PLUGIN_EXAMPLEDIR = "example"
PLUGIN_TEMPLATEDIR = "template"
PLUGIN_SCRIPT_TEMPLATEDIR = "script_template"
PLUGIN_PYFILE = "converter.py"
PLUGIN_DATAFILE = "data.json"

NON_SNAKE_CASE_CHAR_RE = re.compile(r"[^\w_]")


def import_from_path(module_name, file_path):
    spec = importlib.util.spec_from_file_location(module_name, file_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def get_argument_parser():
    """Get an argument parser for this script.

    Returns
    -------
    parser : ArgumentParser
        An argument parser set up with the allowed command-line arguments for this script.
    """

    parser = ArgumentParser()

    parser.add_argument("plugin_name", type=str, nargs="+",
                        help="The name of the plugin to be created, e.g. 'Open Babel'")

    parser.add_argument("--label", type=str, default=None,
                        help="The label for the package (i.e. Python-compatible package name). By default, will "
                        "convert `name` to snake_case (e.g. 'Open Babel' -> 'open_babel')")

    parser.add_argument("--script", action="store_true",
                        help="If set, will create the plugin using the 'ScriptFileConverter' base class, which uses "
                        "a script to run the conversion. The script will by default be named `{label}.sh`")

    return parser


def parse_args():
    """Parses arguments for this executable.

    Returns
    -------
    args : Namespace
        The parsed arguments.
    """

    parser = get_argument_parser()

    args = parser.parse_args()

    return args


def run_from_args(args):
    """Workhorse function to perform primary execution of this script, using the provided parsed arguments.

    Parameters
    ----------
    args : Namespace
        The parsed arguments for this script.
    """

    # Get the name and label from input arguments
    name: str = " ".join(args.plugin_name)
    name = name.strip("'\"")

    label: str | None = args.label
    if label:
        # Check that the label appears to be properly in snake_case
        if label != label.lower().replace(" ", "_") or NON_SNAKE_CASE_CHAR_RE.search(label):
            print(f"{TextColors.FAIL}ERROR:{TextColors.ENDC} Label '{label}' is invalid. The label should be in "
                  "snake_case (all lower-case with underscores in place of spaces), containing only letters and "
                  "underscores", file=sys.stderr)
            exit(1)
    else:
        # Create the label by converting the name to snake_case and stripping invalid characters
        label = NON_SNAKE_CASE_CHAR_RE.sub("", name.lower().replace(" ", "_"))
        if not label:
            print(f"{TextColors.FAIL}ERROR:{TextColors.ENDC} A valid label could not be generated from converter name "
                  f"'{name}'. Please specify a label directly with '--label ...'. The label should be in "
                  "snake_case (all lower-case with underscores in place of spaces), containing only letters and "
                  "underscores", file=sys.stderr)
            exit(1)

    # Determine the PascalCase name of the converter (to be used for classes)
    l_name_words = name.split(" ")
    pascal_name = "".join(map(lambda x: x.capitalize(), l_name_words))

    # Get the directory containing converter plugins
    conv_path = os.path.split(os.path.realpath(converters_base.__file__))[0]

    # Check for any name clashes with existing converters
    for dir in os.listdir(conv_path):
        qual_dir = os.path.join(conv_path, dir)
        if not os.path.isdir(qual_dir) or not os.path.isfile(os.path.join(qual_dir, "__init__.py")):
            continue
        if label == dir:
            print(f"{TextColors.FAIL}ERROR:{TextColors.ENDC} Label '{label}' clashes with the label of an existing "
                  "converter plugin. Please choose a different label (or different name if this was determined from "
                  "the name)", file=sys.stderr)
            exit(1)
        conv_module = import_from_path(label, os.path.join(qual_dir, PLUGIN_PYFILE))
        if name == conv_module.converter.meta.name:
            print(f"{TextColors.FAIL}ERROR:{TextColors.ENDC} Name '{name}' clashes with the name of an existing "
                  "converter plugin. Please choose a different name", file=sys.stderr)
            exit(1)

    # Determine which template directory to use and other related info
    if args.script:
        template_path = os.path.join(conv_path, PLUGIN_SCRIPT_TEMPLATEDIR)
        template_str: str = "ScriptTemplate"
    else:
        template_path = os.path.join(conv_path, PLUGIN_TEMPLATEDIR)
        template_str: str = "Template"

    # Create a new folder for the plugin by copying/modifying from the template appropriately
    plugin_path = os.path.join(conv_path, label)
    os.makedirs(plugin_path)
    shutil.copy(os.path.join(template_path, "__init__.py"), plugin_path)
    for filename in (PLUGIN_PYFILE, PLUGIN_DATAFILE):
        text = open(os.path.join(template_path, filename)).read()

        # Replace template info as appropriate
        if filename == PLUGIN_DATAFILE:
            text = text.replace(template_str, name.replace(r'"', r'\"'))
        else:
            text = text.replace(template_str, pascal_name)

        open(os.path.join(plugin_path, filename), "w").write(text)

    print(f"Success! The plugin has been created at {plugin_path}. Next steps:\n"
          f"- Edit the '{PLUGIN_PYFILE}' and '{PLUGIN_DATAFILE}' files in this directory to contain all necessary "
          "information about this converter and how to run it\n"
          f"- Run the script `psdi-data-convert-install-plugin {label}` to install it (RODO: script under "
          "development)\n"
          "- If this script highlights that formats provided by this plugin may already be in the database,"
          "follow the provided instructions to resolve this and then run it again")


def main():
    """Standard entry-point function for this script.
    """

    args = parse_args()

    run_from_args(args)


if __name__ == "__main__":

    main()
