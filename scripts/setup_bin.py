"""@file scripts/setup_helper.py

Created 2025-03-04 by Bryan Gillis.

Copy appropriate binaries into the project's files before packaging
"""
from typing import Any
from hatchling.builders.hooks.plugin.interface import BuildHookInterface  # type: ignore

import sys
import os
import shutil


class SetupBin(BuildHookInterface):
    PLUGIN_NAME = 'setup binaries'

    def initialize(self, version: str, build_data: dict[str, Any]) -> None:

        platform = sys.platform
        if platform.startswith('linux'):
            platform_folder_name = "linux"
        elif platform.startswith('win'):
            platform_folder_name = "windows"
        elif platform.startswith('darwin'):
            platform_folder_name = "mac"
        else:
            platform_folder_name = None

        source_root_dir = "bin"
        dest_root_dir = "psdi_data_conversion/bin"

        if platform_folder_name is None:
            shutil.copytree(source_root_dir, dest_root_dir, dirs_exist_ok=True)
        else:
            shutil.copytree(os.path.join(source_root_dir, platform_folder_name),
                            os.path.join(dest_root_dir, platform_folder_name), dirs_exist_ok=True)
