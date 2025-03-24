"""
# conversion_test_specs.py

This module contains conversion test specifications, which define conversions to be run and how the results should be
checked. These test specs can be used to test the same conversion in each of the Python library, command-line
application, and GUI.
"""

from psdi_data_conversion import constants as const
from psdi_data_conversion.converters.atomsk import CONVERTER_ATO
from psdi_data_conversion.converters.c2x import CONVERTER_C2X
from psdi_data_conversion.converters.openbabel import CONVERTER_OB
from psdi_data_conversion.testing.conversion_callbacks import (CheckLogContents, CheckLogContentsSuccess,
                                                               CheckOutputStatus, MultiCallback)
from psdi_data_conversion.testing.utils import ConversionTestSpec

default_ob_data = {"coordinates": "neither",
                   "coordOption": "medium",
                   "upload_file": "true"}

basic_tests = ConversionTestSpec(filename=["1NE6.mmcif", "standard_test.cdxml",
                                           "hemoglobin.pdb", "nacl.cif",
                                           "hemoglobin.pdb", "nacl.cif"],
                                 to_format=["pdb", "inchi",
                                            "cif", "xyz",
                                            "cif", "xyz"],
                                 name=[CONVERTER_OB, CONVERTER_OB,
                                       CONVERTER_ATO, CONVERTER_ATO,
                                       CONVERTER_C2X, CONVERTER_C2X],
                                 conversion_kwargs=[{"data": default_ob_data}, {"data": default_ob_data},
                                                    {}, {},
                                                    {}, {},],
                                 post_conversion_callback=MultiCallback([CheckOutputStatus(),
                                                                         CheckLogContentsSuccess()]))
"""A basic set of test conversions which we expect to succeed without issue, running two conversions with each of the
Open Babel, Atomsk, and c2x converters"""

quality_note_callback = CheckLogContentsSuccess(
    l_strings_to_find=["WARNING",
                       const.QUAL_NOTE_OUT_MISSING.format(const.QUAL_2D_LABEL),
                       const.QUAL_NOTE_OUT_MISSING.format(const.QUAL_3D_LABEL),
                       const.QUAL_NOTE_IN_MISSING.format(const.QUAL_CONN_LABEL)])
quality_note_test = ConversionTestSpec(filename="quartz.xyz",
                                       to_format="inchi",
                                       conversion_kwargs={"data": default_ob_data},
                                       post_conversion_callback=quality_note_callback)
"""A test conversion which we expect to produce a warning for conversion quality issues, where the connections property
isn't present in the input and has to be extrapolated, and the 2D and 3D coordinates properties aren't present in the
output and will be lost"""

max_size_test = ConversionTestSpec(filename=["1NE6.mmcif",
                                             "caffeine-smi.tar.gz"],
                                   to_format="pdb",
                                   conversion_kwargs=[{"data": default_ob_data, "max_file_size": 0.0001},
                                                      {"data": default_ob_data, "max_file_size": 0.0005}],
                                   expect_success=False,
                                   post_conversion_callback=CheckLogContents(
                                       l_strings_to_find=["Output file exceeds maximum size"]))
"""A set of test conversion that the maximum size constraint is properly applied. In the first test, the input file
will be greater than the maximum size, and the test should fail as soon as it checks it. In the second test, the input
archive is smaller than the maximum size, but the unpacked files in it are greater, so it should fail midway through."""
