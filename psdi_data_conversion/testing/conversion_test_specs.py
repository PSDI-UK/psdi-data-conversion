"""
# conversion_test_specs.py

This module contains conversion test specifications, which define conversions to be run and how the results should be
checked. These test specs can be used to test the same conversion in each of the Python library, command-line
application, and GUI.
"""

from psdi_data_conversion import constants as const
from psdi_data_conversion.converters.atomsk import CONVERTER_ATO
from psdi_data_conversion.converters.base import (FileConverterAbortException, FileConverterInputException,
                                                  FileConverterSizeException)
from psdi_data_conversion.converters.c2x import CONVERTER_C2X
from psdi_data_conversion.converters.openbabel import CONVERTER_OB, COORD_GEN_KEY, COORD_GEN_QUAL_KEY
from psdi_data_conversion.database import FileConverterDatabaseException
from psdi_data_conversion.testing.conversion_callbacks import (CheckArchiveContents, CheckException, CheckFileStatus,
                                                               CheckLogContents, CheckLogContentsSuccess,
                                                               CheckStderrContents, CheckStdoutContents,
                                                               MatchOutputFile)
from psdi_data_conversion.testing.conversion_callbacks import MultiCallback as MCB
from psdi_data_conversion.testing.utils import ConversionTestSpec as Spec

l_all_test_specs: list[Spec] = []
"""All test specs defined in this file"""

l_all_test_specs.append(Spec(name="Standard Single Test",
                             filename="standard_test.cdxml",
                             to_format="inchi",
                             callback=MCB(CheckFileStatus(),
                                          CheckLogContentsSuccess(),
                                          MatchOutputFile("standard_test.inchi")),
                             ))
"""A quick single test, functioning mostly as a smoke test for things going right in the simplest case"""

simple_success_callback = MCB(CheckFileStatus(), CheckLogContentsSuccess())
l_all_test_specs.append(Spec(name="Standard Multiple Tests",
                             filename=["1NE6.mmcif",
                                       "hemoglobin.pdb", "aceticacid.mol", "nacl.cif",
                                       "hemoglobin.pdb", "hemoglobin.pdb", "nacl.cif",
                                       "hemoglobin.pdb", "hemoglobin.pdb", "nacl.cif",
                                       "ethanol.xyz"],
                             to_format=["pdb-0",
                                        "cif", "mol2", "xyz",
                                        "cif", "xyz", "xyz",
                                        "cif", "xyz-0", "xyz-0",
                                        "cml"],
                             from_format=[None,
                                          None, None, None,
                                          None, None, None,
                                          "pdb-0", "pdb-0", None,
                                          None],
                             converter_name=[CONVERTER_OB,
                                             CONVERTER_OB, CONVERTER_OB, CONVERTER_OB,
                                             CONVERTER_ATO, CONVERTER_ATO, CONVERTER_ATO,
                                             CONVERTER_C2X, CONVERTER_C2X, CONVERTER_C2X,
                                             CONVERTER_OB],
                             callback=simple_success_callback,
                             ))
"""A basic set of test conversions which we expect to succeed without issue, running conversions with each of the
Open Babel, Atomsk, and c2x converters"""

l_all_test_specs.append(Spec(name="c2x Formats Tests",
                             to_format=["res", "abi", "POSCAR", "cml"],
                             converter_name=CONVERTER_C2X,
                             callback=simple_success_callback,
                             compatible_with_gui=False,
                             ))
"""Test converting with c2x to a few different formats which require special input. This test isn't run in the GUI
solely to save on resources, since there are unlikely to be an GUI-specific issues raised by this test that aren't
caught in others."""

l_all_test_specs.append(Spec(name="Converter Name Sensitivity Tests",
                             converter_name=["open babel", "oPeNbaBEL", "C2X", "atomsk"],
                             to_format="xyz-0",
                             callback=simple_success_callback,
                             compatible_with_gui=False,
                             ))
"""Tests that converters can be specified case- and space-insensitively in the library and CLI"""

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
                             from_format="smi",
                             to_format="inchi",
                             callback=archive_callback,
                             ))
"""A test of converting a archives of files"""

l_all_test_specs.append(Spec(name="Archive (wrong format) - Library and CLA",
                             filename="caffeine-smi.zip",
                             to_format="inchi",
                             from_format=["pdb-0", "pdb-0"],
                             conversion_kwargs=[{}, {"strict": True}],
                             expect_success=[True, False],
                             callback=[CheckStderrContents(const.ERR_WRONG_EXTENSIONS),
                                       CheckException(ex_type=FileConverterInputException,
                                                      ex_message=const.ERR_WRONG_EXTENSIONS)],
                             compatible_with_gui=False,
                             ))
"""A test that if the user provides the wrong input format for files in an archive, and error will be output to stderr
"""

l_all_test_specs.append(Spec(name="Archive (wrong format) - GUI",
                             filename="caffeine-smi.zip",
                             to_format="inchi",
                             from_format=["pdb-0", "pdb-0"],
                             conversion_kwargs=[{}, {"strict": True}],
                             expect_success=[False, False],
                             callback=CheckException(ex_type=FileConverterInputException,
                                                     ex_message=const.ERR_WRONG_EXTENSIONS),
                             compatible_with_library=False,
                             compatible_with_cla=False,
                             ))
"""A test that if the user provides the wrong input format for files in an archive - variant for the GUI test, which is
more strict
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
                             to_format="pdb-0",
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
                             compatible_with_gui=False,
                             ))
"""A test that a proper error is returned if an invalid converter is requested

Not compatible with GUI tests, since the GUI only offers valid converters to choose from
"""

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
                             compatible_with_gui=False,
                             ))
"""A test that the input file to a conversion is deleted when cleanup is requested.

Not compatible with the GUI, since the GUI can't forcibly delete files uploaded from the user's computer
"""

l_all_test_specs.append(Spec(name="Failed conversion - bad input file",
                             filename=["quartz_err.xyz", "cyclopropane_err.mol"],
                             to_format=["inchi", "xyz-0"],
                             from_format=[None, "mol-0"],
                             expect_success=False,
                             converter_name=[CONVERTER_OB, CONVERTER_C2X],
                             callback=[MCB(CheckFileStatus(expect_output_exists=False,
                                                           expect_log_exists=None),
                                           CheckException(ex_type=FileConverterAbortException,
                                                          ex_message="Problems reading an XYZ file")),
                                       MCB(CheckFileStatus(expect_output_exists=False,
                                                           expect_log_exists=None),
                                           CheckException(ex_type=FileConverterAbortException,
                                                          ex_message="Aborting: mol file contains no atoms"))],
                             ))
"""A test that a conversion that fails due to an unreadable input file will properly fail"""

quartz_error_ob_callback = CheckLogContents(["ERROR",
                                             "Problems reading an XYZ file: Could not read line #11, file error"])
l_all_test_specs.append(Spec(name="Errors in logs - Library and CLA",
                             filename="quartz_err.xyz",
                             to_format="inchi",
                             converter_name=CONVERTER_OB,
                             expect_success=False,
                             callback=quartz_error_ob_callback,
                             compatible_with_gui=False,
                             ))
"""A test that when a conversion fails in the library or CLA, logs are still produced and contain the expected error
message"""

l_all_test_specs.append(Spec(name="Errors in logs - GUI",
                             filename="quartz_err.xyz",
                             to_format="inchi",
                             converter_name=CONVERTER_OB,
                             expect_success=False,
                             callback=CheckException(ex_type=FileConverterAbortException,
                                                     ex_message=("Problems reading an XYZ file: Could not read line "
                                                                 "#11, file error")),
                             compatible_with_library=False,
                             compatible_with_cla=False,
                             ))
"""A test that when a conversion fails in the GUI, the log message is output to the alert box"""

l_all_test_specs.append(Spec(name="Failed conversion - invalid conversion",
                             filename=["Fapatite.ins", "nacl.mol"],
                             from_format=["ins", "mol-0"],
                             to_format=["cml", "xyz"],
                             expect_success=False,
                             converter_name=[CONVERTER_C2X, CONVERTER_ATO],
                             callback=MCB(CheckFileStatus(expect_output_exists=False,
                                                          expect_log_exists=None),
                                          CheckException(ex_type=FileConverterDatabaseException,
                                                         ex_message="is not supported")),
                             compatible_with_gui=False,
                             ))
"""A test that a conversion that fails due an unsupported conversion will properly fail.

Not compatible with the GUI, since the GUI only offers valid conversions.
"""

l_all_test_specs.append(Spec(name="Blocked conversion - wrong input type",
                             filename="1NE6.mmcif",
                             to_format="cif",
                             from_format="pdb-0",
                             conversion_kwargs={"strict": True},
                             expect_success=False,
                             callback=MCB(CheckFileStatus(expect_output_exists=False,
                                                          expect_log_exists=None),
                                          CheckException(ex_type=FileConverterInputException,
                                                         ex_message=("The file extension is not {} or a zip or tar "
                                                                     "archive extension"))),
                             compatible_with_library=False,
                             compatible_with_cla=False,
                             ))
"""A test that a conversion which is blocked in the GUI"""

l_all_test_specs.append(Spec(name="Failed conversion - wrong input type",
                             filename="1NE6.mmcif",
                             to_format="cif",
                             from_format="pdb-0",
                             conversion_kwargs={"strict": False},
                             expect_success=False,
                             callback=MCB(CheckFileStatus(expect_output_exists=False,
                                                          expect_log_exists=None),
                                          CheckException(ex_type=FileConverterAbortException,
                                                         ex_message=("not a valid {} file"))),
                             ))
"""A test that a conversion which fails due to the wrong input file type will properly fail"""

l_all_test_specs.append(Spec(name="Large files - Library and CLA",
                             filename=["ch3cl-esp.cub", "benzyne.molden", "periodic_dmol3.outmol",
                                       "fullRhinovirus.pdb"],
                             to_format=["cdjson", "dmol", "mol", "cif"],
                             from_format=[None, None, None, "pdb-0"],
                             conversion_kwargs=[{}, {}, {}, {"strict": False}],
                             converter_name=[CONVERTER_OB, CONVERTER_OB, CONVERTER_OB, CONVERTER_C2X],
                             callback=CheckFileStatus(),
                             compatible_with_gui=False,
                             ))
"""Test that the library and CLA can process large files properly"""

l_all_test_specs.append(Spec(name="Large files - GUI",
                             filename=["ch3cl-esp.cub", "benzyne.molden",
                                       "periodic_dmol3.outmol", "fullRhinovirus.pdb"],
                             to_format=["cdjson", "dmol", "mol", "cif"],
                             from_format=[None, None, None, "pdb-0"],
                             converter_name=[CONVERTER_OB, CONVERTER_OB, CONVERTER_OB, CONVERTER_C2X],
                             expect_success=[False, False, False, True],
                             callback=[CheckException(ex_type=FileConverterInputException),
                                       CheckException(ex_type=FileConverterInputException),
                                       CheckException(ex_type=FileConverterInputException),
                                       CheckFileStatus()],
                             compatible_with_library=False,
                             compatible_with_cla=False,
                             ))
"""Test that the GUI will refuse to process large files with OB, but will with other converters"""

max_size_callback = MCB(CheckFileStatus(expect_output_exists=False),
                        CheckLogContents("file exceeds maximum size"),
                        CheckException(ex_type=FileConverterSizeException,
                                       ex_message="exceeds maximum size",
                                       ex_status_code=const.STATUS_CODE_SIZE))
l_all_test_specs.append(Spec(name="Max size exceeded",
                             filename=["1NE6.mmcif", "caffeine-smi.tar.gz"],
                             to_format="pdb-0",
                             conversion_kwargs=[{"max_file_size": 0.0001}, {"max_file_size": 0.0005}],
                             expect_success=False,
                             callback=max_size_callback,
                             compatible_with_cla=False,
                             compatible_with_gui=False,
                             ))
"""A set of test conversion that the maximum size constraint is properly applied. In the first test, the input file
will be greater than the maximum size, and the test should fail as soon as it checks it. In the second test, the input
archive is smaller than the maximum size, but the unpacked files in it are greater, so it should fail midway through.

Not compatible with CLA tests, since the CLA doesn't allow the imposition of a maximum size.

Not compatible with GUI tests in current setup of test implementation, which doesn't let us set env vars to control
things like maximum size on a per-test basis. May be possible to set up in the future though
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

l_library_test_specs = [x for x in l_all_test_specs if x.compatible_with_library and not x.skip_all]
"""All test specs which are compatible with being run on the Python library"""

l_cla_test_specs = [x for x in l_all_test_specs if x.compatible_with_cla and not x.skip_all]
"""All test specs which are compatible with being run on the command-line application"""

l_gui_test_specs = [x for x in l_all_test_specs if x.compatible_with_gui and not x.skip_all]
"""All test specs which are compatible with being run on the GUI"""
