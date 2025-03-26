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
                                                               MatchOutputFile, MultiCallback)
from psdi_data_conversion.testing.utils import ConversionTestSpec

basic_tests = ConversionTestSpec(filename=["1NE6.mmcif", "standard_test.cdxml",
                                           "hemoglobin.pdb", "nacl.cif",
                                           "hemoglobin.pdb", "nacl.cif",
                                           "hemoglobin.pdb", "nacl.cif"],
                                 to_format=["pdb", "inchi", "cif", "xyz",
                                            "cif", "xyz",
                                            "cif", "xyz"],
                                 name=[CONVERTER_OB, CONVERTER_OB,
                                       CONVERTER_OB, CONVERTER_OB,
                                       CONVERTER_ATO, CONVERTER_ATO,
                                       CONVERTER_C2X, CONVERTER_C2X],
                                 callback=MultiCallback(CheckFileStatus(),
                                                        CheckLogContentsSuccess()))
"""A basic set of test conversions which we expect to succeed without issue, running two conversions with each of the
Open Babel, Atomsk, and c2x converters"""

archive_callback = MultiCallback(CheckFileStatus(),
                                 CheckArchiveContents(l_filename_bases=["caffeine-no-flags",
                                                                        "caffeine-ia",
                                                                        "caffeine-ia-ox",
                                                                        "caffeine-ia-okx",
                                                                        "caffeine-ia-okx-oof4",
                                                                        "caffeine-ia-okx-oof4l5",],
                                                      to_format="inchi"))
archive_tests = ConversionTestSpec(filename=["caffeine-smi.zip",
                                             "caffeine-smi.tar",
                                             "caffeine-smi.tar.gz"],
                                   to_format="inchi",
                                   callback=archive_callback)
"""A test of converting a archives of files"""

archive_wrong_format_test = ConversionTestSpec(filename="caffeine-smi.zip",
                                               to_format="inchi",
                                               conversion_kwargs=[{"from_format": "pdb"},
                                                                  {"from_format": "pdb", "strict": True}],
                                               expect_success=[True, False],
                                               callback=[CheckStderrContents(const.ERR_WRONG_EXTENSIONS),
                                                         CheckException(ex_type=FileConverterInputException,
                                                                        ex_message=const.ERR_WRONG_EXTENSIONS)]
                                               )
"""A test that if the user provides the wrong input format for files in an archive, and error will be output to stderr
"""

log_mode_tests = ConversionTestSpec(conversion_kwargs=[{"log_mode": const.LOG_NONE},
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
                                    )
"""Tests that the different log modes have the desired effects on logs

NOTE: Not compatible with GUI tests, since the GUI requires the log mode to always be "Full"
"""

stdout_test = ConversionTestSpec(conversion_kwargs={"log_mode": const.LOG_STDOUT},
                                 callback=CheckStdoutContents(l_strings_to_exclude=["ERROR", "exception", "Exception"],
                                                              l_regex_to_find=[r"File name:\s*nacl",
                                                                               const.DATETIME_RE_RAW]
                                                              ))
"""Test that the log is output to stdout when requested

NOTE: Not compatible with GUI tests, since the GUI requires the log mode to always be "Full"
"""

quiet_test = ConversionTestSpec(conversion_kwargs={"log_mode": const.LOG_NONE},
                                callback=CheckStdoutContents(l_regex_to_exclude=r"."))
"""Test that nothing is output to stdout when quiet mode is enabled

NOTE: Not compatible with GUI tests, since the GUI doesn't support quiet mode
"""

open_babel_warning_test = ConversionTestSpec(filename="1NE6.mmcif",
                                             to_format="pdb",
                                             callback=CheckLogContentsSuccess(["Open Babel Warning",
                                                                               "Failed to kekulize aromatic bonds",])
                                             )
"""A test that confirms expected warnings form Open Babel are output and captured in the log"""

invalid_converter_callback = MultiCallback(CheckFileStatus(expect_output_exists=False,
                                                           expect_log_exists=False),
                                           CheckException(ex_type=FileConverterInputException,
                                                          ex_message="Converter {} not recognized"))
invalid_converter_test = ConversionTestSpec(name="INVALID",
                                            expect_success=False,
                                            callback=invalid_converter_callback)
"""A test that a proper error is returned if an invalid converter is requested"""

quality_note_callback = CheckLogContentsSuccess(["WARNING",
                                                 const.QUAL_NOTE_OUT_MISSING.format(const.QUAL_2D_LABEL),
                                                 const.QUAL_NOTE_OUT_MISSING.format(const.QUAL_3D_LABEL),
                                                 const.QUAL_NOTE_IN_MISSING.format(const.QUAL_CONN_LABEL)])
quality_note_test = ConversionTestSpec(filename="quartz.xyz",
                                       to_format="inchi",
                                       callback=quality_note_callback)
"""A test conversion which we expect to produce a warning for conversion quality issues, where the connections property
isn't present in the input and has to be extrapolated, and the 2D and 3D coordinates properties aren't present in the
output and will be lost"""

cleanup_input_test = ConversionTestSpec(conversion_kwargs={"delete_input": True},
                                        callback=CheckFileStatus(expect_input_exists=False))
"""A test that the input file to a conversion is deleted when cleanup is requested"""

cant_read_xyz_callback = MultiCallback(CheckFileStatus(expect_output_exists=False,
                                                       expect_log_exists=None),
                                       CheckException(ex_type=FileConverterAbortException,
                                                      ex_message="Problems reading an XYZ file"))
invalid_conversion_callback = MultiCallback(CheckFileStatus(expect_output_exists=False,
                                                            expect_log_exists=None),
                                            CheckException(ex_type=FileConverterHelpException,
                                                           ex_message="is not supported"))
wrong_type_callback = MultiCallback(CheckFileStatus(expect_output_exists=False,
                                                    expect_log_exists=None),
                                    CheckException(ex_type=FileConverterAbortException,
                                                   ex_message="not a valid {} file"))
failed_conversion_test = ConversionTestSpec(filename=["quartz_err.xyz", "hemoglobin.pdb", "1NE6.mmcif"],
                                            to_format=["inchi", "pdb", "cif"],
                                            conversion_kwargs=[{}, {}, {"from_format": "pdb"}],
                                            expect_success=False,
                                            callback=[cant_read_xyz_callback, invalid_conversion_callback,
                                                      wrong_type_callback])
"""A test that a conversion which fails due to an invalid input file will properly fail"""

max_size_callback = MultiCallback(CheckFileStatus(expect_output_exists=False),
                                  CheckLogContents("Output file exceeds maximum size"),
                                  CheckException(ex_type=FileConverterSizeException,
                                                 ex_message="exceeds maximum size",
                                                 ex_status_code=const.STATUS_CODE_SIZE))
max_size_test = ConversionTestSpec(filename=["1NE6.mmcif", "caffeine-smi.tar.gz"],
                                   to_format="pdb",
                                   conversion_kwargs=[{"max_file_size": 0.0001}, {"max_file_size": 0.0005}],
                                   expect_success=False,
                                   callback=max_size_callback)
"""A set of test conversion that the maximum size constraint is properly applied. In the first test, the input file
will be greater than the maximum size, and the test should fail as soon as it checks it. In the second test, the input
archive is smaller than the maximum size, but the unpacked files in it are greater, so it should fail midway through.

NOTE: Not compatible with CLA tests, since the CLA doesn't allow the imposition of a maximum size.
"""


format_args_test = ConversionTestSpec(filename="caffeine.inchi",
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
                                      callback=[MatchOutputFile("caffeine-no-flags.smi"),
                                                MatchOutputFile("caffeine-ia.smi"),
                                                MatchOutputFile("caffeine-ia-ox.smi"),
                                                MatchOutputFile("caffeine-ia-okx.smi"),
                                                MatchOutputFile("caffeine-ia-okx-oof4.smi"),
                                                MatchOutputFile("caffeine-ia-okx-oof4l5.smi")
                                                ]
                                      )
"""A set of tests which checks that format args (for how to read from and write to specific file formats) are processed
correctly, by matching tests using them to expected output files"""


coord_gen_test = ConversionTestSpec(filename="caffeine.inchi",
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
                                    )
"""A set of tests which checks that coordinate generation options are processed correctly, by matching tests using them
to expected output files"""
