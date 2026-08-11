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

from psdi_data_conversion.testing.constants import TEST_PATH_KEY
from psdi_data_conversion.testing.utils import get_test_data_loc
from psdi_data_conversion.utils import TextColors as TC
from psdi_data_conversion.utils import confirm_editable_mode, get_project_path, print_wrap

PLUGIN_EXAMPLEDIR = "example"
PLUGIN_TEMPLATEDIR = "template"
PLUGIN_SCRIPT_TEMPLATEDIR = "script_template"
PLUGIN_PYFILE = "converter.py"
PLUGIN_DATAFILE = "data.json"

NON_SNAKE_CASE_CHAR_RE = re.compile(r"[^a-zA-Z0-9_]")

# Environmental variable to alter functionality of this script to retrieve a premade test converter plugin from the
# `test_data` folder. In this mode, the required `plugin_name` argument will instead be used as the name of the
# subfolder for the test plugin within `test_data`, e.g. `TEST_DATA=true psdi-data-convert-create-plugin test_converter`
# will retrieve the test plugin at `test_data/test_converter`."
TEST_DATA_KEY = "TEST_DATA"


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

    # Check if the package is installed in editable mode
    confirm_editable_mode()

    # Get the name and label from input arguments
    name: str = " ".join(args.plugin_name)
    name = name.strip("'\"").strip()

    label: str | None = args.label
    if label:
        # Check that the label appears to be properly in snake_case
        if label != label.lower().replace(" ", "_") or NON_SNAKE_CASE_CHAR_RE.search(label):
            print_wrap(f"{TC.FAIL}ERROR:{TC.ENDC} Label '{label}' is invalid. The label should be in snake_case (all "
                       "lower-case with underscores in place of spaces), containing only letters, digits, "
                       "and underscores", err=True)
            exit(1)
    else:
        # Create the label by converting the name to snake_case and stripping invalid characters
        label = NON_SNAKE_CASE_CHAR_RE.sub("", name.lower().replace(" ", "_"))
        if not label:
            print_wrap(f"{TC.FAIL}ERROR:{TC.ENDC} A valid label could not be generated from converter name '{name}'. "
                       f"Please specify a label directly with '{TC.WARNING}--label LABEL{TC.ENDC}'. The label should "
                       "be in snake_case (all lower-case with underscores in place of spaces), containing only "
                       "letters, digits, and underscores", err=True)
            exit(1)

    # Determine the PascalCase name of the converter (to be used for classes)
    l_name_words = name.split(" ")
    pascal_name = "".join(map(lambda x: x.capitalize(), l_name_words))

    # Get the project path to use based on if we're using a test path or not
    if os.environ[TEST_PATH_KEY]:
        project_path: str = os.path.realpath(os.environ[TEST_PATH_KEY])
        if not os.path.isdir(project_path):
            print_wrap(f"{TC.FAIL}ERROR:{TC.ENDC} When running this script with '{TC.WARNING}--test-path TEST_PATH" +
                       f"{TC.ENDC}', the provided path ({TC.OKCYAN}{project_path}{TC.ENDC}) must already exist.",
                       err=True)
            exit(1)
    else:
        project_path = get_project_path()

    # Get the directory containing converter plugins
    conv_path = os.path.join(project_path, "psdi_data_conversion", "converters")

    # Check for any name clashes with existing converters
    for dir in os.listdir(conv_path):
        qual_dir = os.path.join(conv_path, dir)
        if not os.path.isdir(qual_dir) or not os.path.isfile(os.path.join(qual_dir, "__init__.py")):
            continue
        if label == dir:
            if os.environ.get(TEST_DATA_KEY):
                # In test mode, don't worry if it already exists, just delete it so it doesn't clash
                shutil.rmtree(qual_dir)
                break
            else:
                print_wrap(f"{TC.FAIL}ERROR:{TC.ENDC} Label '{label}' clashes with the label of an "
                           "existing converter plugin. Please choose a different label (or different name if this was "
                           "determined from the name)", err=True)
                exit(1)
        conv_module = import_from_path(label, os.path.join(qual_dir, PLUGIN_PYFILE))
        if name == conv_module.converter.meta.name:
            print_wrap(f"{TC.FAIL}ERROR:{TC.ENDC} Name '{name}' clashes with the name of an existing "
                       "converter plugin. Please choose a different name", err=True)
            exit(1)

    # Determine which template directory to use and other related info
    if os.environ.get(TEST_DATA_KEY):
        template_path = os.path.join(get_test_data_loc(os.environ[TEST_DATA_KEY]), name)
        template_str: str = "Template"
    elif args.script:
        template_path = os.path.join(conv_path, PLUGIN_SCRIPT_TEMPLATEDIR)
        template_str: str = "ScriptTemplate"
    else:
        template_path = os.path.join(conv_path, PLUGIN_TEMPLATEDIR)
        template_str: str = "Template"

    # Create a new folder for the plugin by copying/modifying from the template appropriately
    plugin_path = os.path.join(conv_path, label)
    os.makedirs(plugin_path, exist_ok=True)
    shutil.copy(os.path.join(template_path, "__init__.py"), plugin_path)
    for filename in (PLUGIN_PYFILE, PLUGIN_DATAFILE):
        text = open(os.path.join(template_path, filename)).read()

        if not os.environ.get(TEST_DATA_KEY):
            # Replace template info as appropriate
            if filename == PLUGIN_DATAFILE:
                text = text.replace(template_str, name.replace(r'"', r'\"'))
            else:
                # Where the template string appears alone as part of a word, replace it with the name, otherwise
                # replace it with the Pascal-case name
                text = re.sub(f"\b{template_str}\b", name, text)
                text = text.replace(template_str, pascal_name)

        open(os.path.join(plugin_path, filename), "w").write(text)

    print(f"{TC.OKGREEN}Success!{TC.ENDC} The plugin has been created at "
          f"{TC.OKCYAN}{plugin_path}{TC.ENDC}\nNext steps:\n")
    print_wrap(f"- Edit the '{TC.OKCYAN}{PLUGIN_PYFILE}{TC.ENDC}' and "
               f"'{TC.OKCYAN}{PLUGIN_DATAFILE}{TC.ENDC}' files in this directory to contain all "
               "necessary information about this converter and how to run it\n",
               initial_indent="", subsequent_indent=" "*2)
    print_wrap(f"- Run the script '{TC.WARNING}psdi-data-convert-install-plugins{TC.ENDC}' to install "
               "it\n",
               initial_indent="", subsequent_indent=" "*2)
    print_wrap("- If this script highlights that formats provided by this plugin may already be in the database,"
               "follow the provided instructions to resolve this and then run it again",
               initial_indent="", subsequent_indent=" "*2)


def main():
    """Standard entry-point function for this script.
    """

    args = parse_args()

    run_from_args(args)


if __name__ == "__main__":

    main()
