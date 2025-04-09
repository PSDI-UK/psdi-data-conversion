"""
# gui.py

Utilities to aid in testing of the GUI
"""


from dataclasses import dataclass
import os
from tempfile import TemporaryDirectory

import pytest
from psdi_data_conversion.constants import LOG_FULL
from psdi_data_conversion.converters.openbabel import (COORD_GEN_KEY, COORD_GEN_QUAL_KEY, DEFAULT_COORD_GEN,
                                                       DEFAULT_COORD_GEN_QUAL)
from psdi_data_conversion.testing.utils import (ConversionTestInfo, ConversionTestSpec, SingleConversionTestSpec,
                                                get_input_test_data_loc)


@dataclass
class GuiConversionTestInfo(ConversionTestInfo):
    """Information about a tested conversion, specifically for when it was tested through the GUI (the local version of
    the web app)
    """


def run_test_conversion_with_gui(test_spec: ConversionTestSpec):
    """Runs a test conversion or series thereof through the GUI. Note that this requires the server to be started before
    this is called.

    Parameters
    ----------
    test_spec : ConversionTestSpec
        The specification for the test or series of tests to be run
    """
    # Make temporary directories for the input and output files to be stored in
    with TemporaryDirectory("_input") as input_dir, TemporaryDirectory("_output") as output_dir:
        # Iterate over the test spec to run each individual test it defines
        for single_test_spec in test_spec:
            _run_single_test_conversion_with_gui(test_spec=single_test_spec,
                                                 input_dir=input_dir,
                                                 output_dir=output_dir)


def _run_single_test_conversion_with_gui(test_spec: SingleConversionTestSpec,
                                         input_dir: str,
                                         output_dir: str):
    """Runs a single test conversion through the GUI.

    Parameters
    ----------
    test_spec : _SingleConversionTestSpec
        The specification for the test to be run
    input_dir : str
        A directory which can be used to store input data before uploading
    output_dir : str
        A directory which can be used to create output data after downloading
    """

    # Symlink the input file to the input directory
    qualified_in_filename = os.path.realpath(os.path.join(input_dir, test_spec.filename))
    try:
        os.symlink(os.path.join(get_input_test_data_loc(), test_spec.filename),
                   qualified_in_filename)
    except FileExistsError:
        pass

    success = run_converter_through_gui(filename=qualified_in_filename,
                                        to_format=test_spec.to_format,
                                        name=test_spec.converter_name,
                                        input_dir=input_dir,
                                        output_dir=output_dir,
                                        log_file=os.path.join(output_dir, test_spec.log_filename),
                                        **test_spec.conversion_kwargs)

    # Compile output info for the test and call the callback function if one is provided
    if test_spec.callback:
        test_info = GuiConversionTestInfo(test_spec=test_spec,
                                          input_dir=input_dir,
                                          output_dir=output_dir,
                                          success=success)
        callback_msg = test_spec.callback(test_info)
        assert not callback_msg, callback_msg


def run_converter_through_gui(filename: str,
                              to_format: str,
                              name: str,
                              input_dir: str,
                              output_dir: str,
                              log_file: str,
                              **conversion_kwargs):
    """_summary_

    Parameters
    ----------
    filename : str
        The (unqualified) name of the input file to be converted
    to_format : str
        The format to convert the input file to
    name : str
        The name of the converter to use
    input_dir : str
        The directory which contains the input file
    output_dir : str
        The directory which contains the output file
    log_file : str
        The desired name of the log file
    """

    # Default options for conversion
    from_format = os.path.splitext(filename)[1]
    strict = True
    from_flags = None
    to_flags = None
    from_options = None
    to_options = None
    coord_gen = DEFAULT_COORD_GEN
    coord_gen_qual = DEFAULT_COORD_GEN_QUAL

    # For each argument in the conversion kwargs, interpret it as the appropriate option for this conversion, overriding
    # defaults set above
    for key, val in conversion_kwargs.items():
        if key == "from_format":
            from_format = val
        elif key == "log_mode":
            raise ValueError(f"The conversion kwarg {key} is not valid with conversions through the GUI")
        elif key == "delete_input":
            raise ValueError(f"The conversion kwarg {key} is not valid with conversions through the GUI")
        elif key == "strict":
            strict = val
        elif key == "max_file_size":
            raise ValueError(f"The conversion kwarg {key} is not valid with conversions through the GUI")
        elif key == "data":
            for subkey, subval in val.items():
                if subkey == "from_flags":
                    from_flags = subval
                elif subkey == "to_flags":
                    to_flags = subval
                elif subkey == "from_options":
                    from_options = subval
                elif subkey == "to_options":
                    to_options = subval
                elif subkey == COORD_GEN_KEY:
                    coord_gen = subval
                    if COORD_GEN_QUAL_KEY in val:
                        coord_gen_qual = val[COORD_GEN_QUAL_KEY]
                elif subkey == COORD_GEN_QUAL_KEY:
                    # Handled alongside COORD_GEN_KEY above
                    pass
                else:
                    pytest.fail(f"The key 'data[\"{subkey}\"]' was passed to `conversion_kwargs` but could not be "
                                "interpreted")
        else:
            pytest.fail(f"The key '{key}' was passed to `conversion_kwargs` but could not be interpreted")

    # Cleanup of arguments
    if from_format.startswith("."):
        from_format = from_format[1:]

    return True
