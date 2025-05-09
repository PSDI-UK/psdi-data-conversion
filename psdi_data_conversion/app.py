"""app.py

Version 1.0, 8th November 2024

This script acts as a server for the PSDI Data Conversion Service website.
"""

from argparse import ArgumentParser

from psdi_data_conversion import constants as const
from psdi_data_conversion.gui.env import update_env
from psdi_data_conversion.gui.setup import get_app, limit_upload_size, start_app
from psdi_data_conversion.main import print_wrap

app = get_app()


def main():
    """Standard entry-point function for this script.
    """

    parser = ArgumentParser()

    parser.add_argument("--use-env-vars", action="store_true",
                        help="If set, all other arguments and defaults for this script are ignored, and environmental "
                        "variables and their defaults will instead control execution. These defaults will result in "
                        "the app running in production server mode.")

    parser.add_argument("--max-file-size", type=float, default=const.DEFAULT_MAX_FILE_SIZE/const.MEGABYTE,
                        help="The maximum allowed filesize in MB - 0 (default) indicates no maximum")

    parser.add_argument("--max-file-size-ob", type=float, default=const.DEFAULT_MAX_FILE_SIZE_OB/const.MEGABYTE,
                        help="The maximum allowed filesize in MB for the Open Babel converter, taking precendence over "
                        "the general maximum file size when Open Babel is used - 0 indicates no maximum. Default 1 MB.")

    parser.add_argument("--service-mode", action="store_true",
                        help="If set, will run as if deploying a service rather than the local GUI")

    parser.add_argument("--dev-mode", action="store_true",
                        help="If set, will expose development elements, such as the SHA of the latest commit")

    parser.add_argument("--debug", action="store_true",
                        help="If set, will run the Flask server in debug mode, which will cause it to automatically "
                        "reload if code changes and show an interactive debugger in the case of errors")

    parser.add_argument("--log-mode", type=str, default=const.LOG_FULL,
                        help="How logs should be stored. Allowed values are: \n"
                        "- 'full' - Multi-file logging, not recommended for the CLI, but allowed for a compatible "
                        "interface with the public web app"
                        "- 'simple' - Logs saved to one file"
                        "- 'stdout' - Output logs and errors only to stdout"
                        "- 'none' - Output only errors to stdout")

    parser.add_argument("--log-level", type=str, default=None,
                        help="The desired level to log at. Allowed values are: 'DEBUG', 'INFO', 'WARNING', 'ERROR, "
                             "'CRITICAL'. Default: 'INFO' for logging to file, 'WARNING' for logging to stdout")

    # Set global variables for settings based on parsed arguments, unless it's set to use env vars
    args = parser.parse_args()

    if not args.use_env_vars:
        # Overwrite the values from environmental variables with the values from the command-line arguments
        update_env(args)

    # Set the upload limit based on provided arguments now
    limit_upload_size()

    print_wrap("Starting the PSDI Data Conversion GUI. This GUI is run as a webpage, which you can open by "
               "right-clicking the link below to open it in your default browser, or by copy-and-pasting it into your "
               "browser of choice.")

    start_app()


if __name__ == "__main__":

    main()
