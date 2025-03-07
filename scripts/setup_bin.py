"""@file scripts/setup_helper.py

Created 2025-03-04 by Bryan Gillis.

Copy appropriate binaries into the project's files before packaging
"""
from tempfile import TemporaryDirectory
from typing import Any
from hatchling.builders.hooks.plugin.interface import BuildHookInterface  # type: ignore

import sys
import os
import shutil

base_dir = os.path.abspath(os.path.split(__file__)[0] + "/..")
source_bin_dir = os.path.join(base_dir, "psdi_data_conversion/bin")

L_PLATFORMS = ["linux", "windows", "mac"]


class HideUnusedBins(BuildHookInterface):
    """A Hatch plugin which hides binaries for platforms other than the current platform from the build, to reduce the
    size of the generated sdists and wheels. It does this by moving the binaries away to a temporary directory during
    the build so they won't be discovered, then moving them back at the end.
    """

    PLUGIN_NAME = 'Hide unused binaries'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.tmpdir = TemporaryDirectory("bin-storage")
        self.platform: str | None = None

    def __del__(self):
        self.restore_bins()

    def initialize(self, version: str, build_data: dict[str, Any]) -> None:

        platform = sys.platform
        if platform.startswith('linux'):
            self.platform = "linux"
        elif platform.startswith('win'):
            self.platform = "windows"
        elif platform.startswith('darwin'):
            self.platform = "mac"
        else:
            self.platform = None

        if self.platform is not None:
            for platform in L_PLATFORMS:
                source_platform_bin_dir = os.path.join(source_bin_dir, platform)
                if platform != self.platform and os.path.exists(source_platform_bin_dir):
                    shutil.move(source_platform_bin_dir,
                                os.path.join(self.tmpdir.name, platform))

    def finalize(self, version: str, build_data: dict[str, Any], artifact_path: str) -> None:
        self.restore_bins()

    def restore_bins(self):
        if self.platform is None or self.tmpdir is None:
            return
        for platform in L_PLATFORMS:
            source_platform_bin_dir = os.path.join(self.tmpdir.name, platform)
            if platform != self.platform and os.path.exists(source_platform_bin_dir):
                shutil.move(source_platform_bin_dir,
                            os.path.join(source_bin_dir, platform))
