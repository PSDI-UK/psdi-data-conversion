"""
# conversion_test_specs.py

This module contains conversion test specifications, which define conversions to be run and how the results should be
checked. These test specs can be used to test the same conversion in each of the Python library, command-line
application, and GUI.
"""

from psdi_data_conversion.converters.atomsk import CONVERTER_ATO
from psdi_data_conversion.converters.c2x import CONVERTER_C2X
from psdi_data_conversion.converters.openbabel import CONVERTER_OB
from psdi_data_conversion.testing.conversion_callbacks import CheckOutputStatus
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
                                 post_conversion_callback=CheckOutputStatus())
