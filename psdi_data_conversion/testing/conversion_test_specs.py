"""
# conversion_test_specs.py

This module contains conversion test specifications, which define conversions to be run and how the results should be
checked. These test specs can be used to test the same conversion in each of the Python library, command-line
application, and GUI.
"""

from psdi_data_conversion import constants as const
from psdi_data_conversion.converters.atomsk import CONVERTER_ATO
from psdi_data_conversion.converters.base import (FileConverterAbortException, FileConverterHelpException,
                                                  FileConverterInputException, FileConverterSizeException)
from psdi_data_conversion.converters.c2x import CONVERTER_C2X
from psdi_data_conversion.converters.openbabel import CONVERTER_OB, COORD_GEN_KEY, COORD_GEN_QUAL_KEY
from psdi_data_conversion.testing.conversion_callbacks import (CheckArchiveContents, CheckException, CheckLogContents,
                                                               CheckLogContentsSuccess, CheckFileStatus,
                                                               CheckStderrContents, CheckStdoutContents,
                                                               MatchOutputFile, MultiCallback as MCB)
from psdi_data_conversion.testing.utils import ConversionTestSpec as Spec


l_all_test_specs: list[Spec] = []
"""All test specs defined in this file"""

l_all_test_specs.append(Spec(name="Basic",
                             filename=["1NE6.mmcif", "standard_test.cdxml",
                                                     "hemoglobin.pdb", "nacl.cif",
                                                     "hemoglobin.pdb", "nacl.cif",
                                                     "hemoglobin.pdb", "nacl.cif"],
                             to_format=["pdb", "inchi", "cif", "xyz",
                                        "cif", "xyz",
                                        "cif", "xyz"],
                             converter_name=[CONVERTER_OB, CONVERTER_OB,
                                             CONVERTER_OB, CONVERTER_OB,
                                             CONVERTER_ATO, CONVERTER_ATO,
                                             CONVERTER_C2X, CONVERTER_C2X],
                             callback=MCB(CheckFileStatus(),
                                          CheckLogContentsSuccess()),
                             ))
"""A basic set of test conversions which we expect to succeed without issue, running two conversions with each of the
Open Babel, Atomsk, and c2x converters"""

archive_callback = MCB(CheckFileStatus(),
                       CheckArchiveContents(l_filename_bases=["caffeine-no-flags",
                                                              "caffeine-ia",
                                                              "caffeine-ia-ox",
                                                              "caffeine-ia-okx",
                                                              "caffeine-ia-okx-oof4",
                                                              "caffeine-ia-okx-oof4l5",],
                                            to_format="inchi"))

l_all_test_specs.append(Spec(name="Archive",
                             filename=["caffeine-smi.zip",
                                       "caffeine-smi.tar",
                                       "caffeine-smi.tar.gz"],
                             to_format="inchi",
                             callback=archive_callback,
                             ))
"""A test of converting a archives of files"""

l_all_test_specs.append(Spec(name="Archive (wrong format)",
                             filename="caffeine-smi.zip",
                             to_format="inchi",
                             conversion_kwargs=[{"from_format": "pdb"},
                                                {"from_format": "pdb", "strict": True}],
                             expect_success=[True, False],
                             callback=[CheckStderrContents(const.ERR_WRONG_EXTENSIONS),
                                       CheckException(ex_type=FileConverterInputException,
                                                      ex_message=const.ERR_WRONG_EXTENSIONS)],
                             ))
"""A test that if the user provides the wrong input format for files in an archive, and error will be output to stderr
"""

l_all_test_specs.append(Spec(name="Log mode",
                             conversion_kwargs=[{"log_mode": const.LOG_NONE},
                                                {"log_mode": const.LOG_STDOUT},
                                                {"log_mode": const.LOG_SIMPLE},
                                                {"log_mode": const.LOG_FULL},
                                                {"log_mode": const.LOG_FULL_FORCE},],
                             callback=[CheckFileStatus(expect_log_exists=False,
                                                       expect_global_log_exists=False),
                                       CheckFileStatus(expect_log_exists=False,
                                                       expect_global_log_exists=False),
                                       CheckFileStatus(expect_log_exists=True,
                                                       expect_global_log_exists=False),
                                       CheckFileStatus(expect_log_exists=True,
                                                       expect_global_log_exists=True),
                                       CheckFileStatus(expect_log_exists=True,
                                                       expect_global_log_exists=True)],
                             compatible_with_gui=False,
                             ))
"""Tests that the different log modes have the desired effects on logs

Not compatible with GUI tests, since the GUI requires the log mode to always be "Full"
"""

l_all_test_specs.append(Spec(name="Stdout",
                             conversion_kwargs={"log_mode": const.LOG_STDOUT},
                             callback=CheckStdoutContents(l_strings_to_exclude=["ERROR", "exception",
                                                                                "Exception"],
                                                          l_regex_to_find=[r"File name:\s*nacl",
                                                                           const.DATETIME_RE_RAW]
                                                          ),
                             compatible_with_gui=False,
                             ))
"""Test that the log is output to stdout when requested

Not compatible with GUI tests, since the GUI requires the log mode to always be "Full"
"""

l_all_test_specs.append(Spec(name="Quiet",
                             conversion_kwargs={"log_mode": const.LOG_NONE},
                             callback=CheckStdoutContents(l_regex_to_exclude=r"."),
                             compatible_with_gui=False,
                             ))
"""Test that nothing is output to stdout when quiet mode is enabled

Not compatible with GUI tests, since the GUI doesn't support quiet mode
"""

l_all_test_specs.append(Spec(name="Open Babel Warning",
                             filename="1NE6.mmcif",
                             to_format="pdb",
                             callback=CheckLogContentsSuccess(["Open Babel Warning",
                                                               "Failed to kekulize aromatic bonds",])
                             ))
"""A test that confirms expected warnings form Open Babel are output and captured in the log"""

invalid_converter_callback = MCB(CheckFileStatus(expect_output_exists=False,
                                                 expect_log_exists=False),
                                 CheckException(ex_type=FileConverterInputException,
                                                ex_message="Converter {} not recognized"))
l_all_test_specs.append(Spec(name="Invalid Converter",
                             converter_name="INVALID",
                             expect_success=False,
                             callback=invalid_converter_callback,
                             ))
"""A test that a proper error is returned if an invalid converter is requested"""

quartz_quality_note_callback = CheckLogContentsSuccess(["WARNING",
                                                        const.QUAL_NOTE_OUT_MISSING.format(const.QUAL_2D_LABEL),
                                                        const.QUAL_NOTE_OUT_MISSING.format(const.QUAL_3D_LABEL),
                                                        const.QUAL_NOTE_IN_MISSING.format(const.QUAL_CONN_LABEL)])
ethanol_quality_note_callback = CheckLogContentsSuccess(["WARNING",
                                                         "Potential data loss or extrapolation",
                                                         const.QUAL_NOTE_IN_MISSING.format(const.QUAL_CONN_LABEL)])
hemoglobin_quality_note_callback = CheckLogContentsSuccess(["WARNING",
                                                            "Potential data loss or extrapolation",
                                                            const.QUAL_NOTE_OUT_MISSING.format(const.QUAL_CONN_LABEL)])
l_all_test_specs.append(Spec(name="Quality note",
                             filename=["quartz.xyz", "ethanol.xyz", "hemoglobin.pdb"],
                             to_format=["inchi", "cml", "xyz"],
                             callback=[quartz_quality_note_callback,
                                       ethanol_quality_note_callback,
                                       hemoglobin_quality_note_callback],
                             ))
"""A test conversion which we expect to produce a warning for conversion quality issues, where the connections property
isn't present in the input and has to be extrapolated, and the 2D and 3D coordinates properties aren't present in the
output and will be lost"""

l_all_test_specs.append(Spec(name="Cleanup input",
                             conversion_kwargs={"delete_input": True},
                             callback=CheckFileStatus(expect_input_exists=False),
                             ))
"""A test that the input file to a conversion is deleted when cleanup is requested"""

cant_read_xyz_callback = MCB(CheckFileStatus(expect_output_exists=False,
                                             expect_log_exists=None),
                             CheckException(ex_type=FileConverterAbortException,
                                            ex_message="Problems reading an XYZ file"))
invalid_conversion_callback = MCB(CheckFileStatus(expect_output_exists=False,
                                                  expect_log_exists=None),
                                  CheckException(ex_type=FileConverterHelpException,
                                                 ex_message="is not supported"))
wrong_type_callback = MCB(CheckFileStatus(expect_output_exists=False,
                                          expect_log_exists=None),
                          CheckException(ex_type=FileConverterAbortException,
                                         ex_message="not a valid {} file"))
l_all_test_specs.append(Spec(name="Failed Conversion",
                             filename=["quartz_err.xyz", "hemoglobin.pdb", "1NE6.mmcif"],
                             to_format=["inchi", "pdb", "cif"],
                             conversion_kwargs=[{}, {}, {"from_format": "pdb"}],
                             expect_success=False,
                             callback=[cant_read_xyz_callback, invalid_conversion_callback,
                                       wrong_type_callback],
                             ))
"""A test that a conversion which fails due to an invalid input file will properly fail"""

max_size_callback = MCB(CheckFileStatus(expect_output_exists=False),
                        CheckLogContents("file exceeds maximum size"),
                        CheckException(ex_type=FileConverterSizeException,
                                       ex_message="exceeds maximum size",
                                       ex_status_code=const.STATUS_CODE_SIZE))
l_all_test_specs.append(Spec(name="Max size exceeded",
                             filename=["1NE6.mmcif", "caffeine-smi.tar.gz"],
                             to_format="pdb",
                             conversion_kwargs=[{"max_file_size": 0.0001}, {"max_file_size": 0.0005}],
                             expect_success=False,
                             callback=max_size_callback,
                             compatible_with_cla=False,
                             ))
"""A set of test conversion that the maximum size constraint is properly applied. In the first test, the input file
will be greater than the maximum size, and the test should fail as soon as it checks it. In the second test, the input
archive is smaller than the maximum size, but the unpacked files in it are greater, so it should fail midway through.

Not compatible with CLA tests, since the CLA doesn't allow the imposition of a maximum size.
"""


l_all_test_specs.append(Spec(name="Format args",
                             filename="caffeine.inchi",
                             to_format="smi",
                             conversion_kwargs=[{},
                                                {"data": {"from_flags": "a"}},
                                                {"data": {"from_flags": "a", "to_flags": "x"}},
                                                {"data": {"from_flags": "a", "to_flags": "kx"}},
                                                {"data": {"from_flags": "a", "to_flags": "kx",
                                                                        "to_options": "f4"}},
                                                {"data": {"from_flags": "a", "to_flags": "kx",
                                                                        "to_options": "f4 l5"}}
                                                ],
                             callback=[MatchOutputFile("caffeine.smi"),
                                       MatchOutputFile("caffeine_a_in.smi"),
                                       MatchOutputFile("caffeine_a_in_x_out.smi"),
                                       MatchOutputFile("caffeine_a_in_kx_out.smi"),
                                       MatchOutputFile("caffeine_a_in_kx_f4_out.smi"),
                                       MatchOutputFile("caffeine_a_in_kx_f4_l5_out.smi")
                                       ]
                             ))
"""A set of tests which checks that format args (for how to read from and write to specific file formats) are processed
correctly, by matching tests using them to expected output files"""


l_all_test_specs.append(Spec(name="Coord gen",
                             filename="caffeine.inchi",
                             to_format="xyz",
                             conversion_kwargs=[{},
                                                {"data": {COORD_GEN_KEY: "Gen2D",
                                                          COORD_GEN_QUAL_KEY: "fastest"}},
                                                {"data": {COORD_GEN_KEY: "Gen3D",
                                                          COORD_GEN_QUAL_KEY: "best"}}
                                                ],
                             callback=[MatchOutputFile("caffeine.xyz"),
                                       MatchOutputFile("caffeine-2D-fastest.xyz"),
                                       MatchOutputFile("caffeine-3D-best.xyz"),
                                       ]
                             ))
"""A set of tests which checks that coordinate generation options are processed correctly, by matching tests using them
to expected output files"""

l_library_test_specs = [x for x in l_all_test_specs if x.compatible_with_library]
"""All test specs which are compatible with being run on the Python library"""

l_cla_test_specs = [x for x in l_all_test_specs if x.compatible_with_cla]
"""All test specs which are compatible with being run on the command-line application"""

l_gui_test_specs = [x for x in l_all_test_specs if x.compatible_with_gui]
"""All test specs which are compatible with being run on the GUI"""
