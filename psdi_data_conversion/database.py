"""@file psdi_data_conversion/database.py

Created 2025-02-03 by Bryan Gillis.

Python module provide utilities for accessing the converter database
"""

from __future__ import annotations

from copy import copy
from dataclasses import dataclass, field
import json
from logging import getLogger
import os
from typing import Any

from psdi_data_conversion import constants as const
from psdi_data_conversion.converter import D_REGISTERED_CONVERTERS
from psdi_data_conversion.converters.base import FileConverterException

logger = getLogger(__name__)


class FileConverterDatabaseException(FileConverterException):
    """Class for any exceptions which arise from issues with the database classes and methods
    """
    pass


@dataclass
class ArgInfo:
    """Class providing information on an argument accepted by a converter (whether it accepts a value or not)
    """

    parent: ConverterInfo
    id: int
    flag: str
    description: str
    info: str

    s_in_formats: set[int] = field(default_factory=set)
    s_out_formats: set[int] = field(default_factory=set)


@dataclass
class FlagInfo(ArgInfo):
    """Class providing information on a flag accepted by a converter (an argument which doesn't accept a value)
    """
    pass


@dataclass
class OptionInfo(ArgInfo):
    """Class providing information on an option accepted by a converter (an argument accepts a value)
    """
    brief: str


class ConverterInfo:
    """Class providing information on a converter stored in the PSDI Data Conversion database
    """

    def __init__(self,
                 name: str,
                 parent: DataConversionDatabase,
                 d_single_converter_info: dict[str, int | str],
                 d_data: dict[str, Any]):
        """Set up the class - this will be initialised within a `DataConversionDatabase`, which we set as the parent

        Parameters
        ----------
        name : str
            The name of the converter
        parent : DataConversionDatabase
            The database which this belongs to
        d_data : dict[str, Any]
            The loaded database dict
        """

        self.name = name
        self.parent = parent

        # Get info about the converter from the database
        self.id: int = d_single_converter_info.get(const.DB_ID_KEY, -1)
        self.description: str = d_single_converter_info.get(const.DB_DESC_KEY, "")
        self.url: str = d_single_converter_info.get(const.DB_URL_KEY, "")

        # Get necessary info about the converter from the class
        try:
            self._key_prefix = D_REGISTERED_CONVERTERS[name].database_key_prefix
        except KeyError:
            # We'll get a KeyError for converters in the database that don't yet have their own class, which we can
            # safely ignore
            self._key_prefix = None

        self._arg_info: dict[str, list[dict[str, int | str]]] = {}

        # Placeholders for members that are generated when needed
        self._l_flag_info: list[FlagInfo | None] | None = None
        self._d_flag_info: dict[str, FlagInfo] | None = None
        self._l_option_info: list[OptionInfo | None] | None = None
        self._d_option_info: dict[str, OptionInfo] | None = None
        self._d_arg_info: dict[str, ArgInfo] | None = None

        # If the converter class has no defined key prefix, don't add any extra info for it
        if self._key_prefix is None:
            return
        for key_base in (const.DB_IN_FLAGS_KEY_BASE,
                         const.DB_OUT_FLAGS_KEY_BASE,
                         const.DB_IN_OPTIONS_KEY_BASE,
                         const.DB_OUT_OPTIONS_KEY_BASE,
                         const.DB_IN_FLAGS_FORMATS_KEY_BASE,
                         const.DB_OUT_FLAGS_FORMATS_KEY_BASE,
                         const.DB_IN_OPTIONS_FORMATS_KEY_BASE,
                         const.DB_OUT_OPTIONS_FORMATS_KEY_BASE):
            self._arg_info[key_base] = d_data.get(self._key_prefix + key_base)

    @property
    def d_flag_info(self) -> dict[str, FlagInfo] | None:
        """Generate the flag info dict (index by flag) when needed. Returns None if the converter has no flag info in
        the database
        """
        if self._d_flag_info is None and self._key_prefix is not None:
            # Load from the database
            self._d_flag_info = {}
            for d_single_flag_info in self._arg_info[const.DB_IN_FLAGS_KEY_BASE]:
                flag: str = d_single_flag_info[const.DB_FLAG_KEY]
                flag_id: int = d_single_flag_info[const.DB_ID_KEY]
                flag_info = FlagInfo(parent=self,
                                     id=flag_id,
                                     flag=flag,
                                     description=d_single_flag_info[const.DB_DESC_KEY],
                                     info=d_single_flag_info[const.DB_INFO_KEY])
                self._d_flag_info[flag] = flag_info

                # Get a list of all in and formats applicable to this flag, and add them to the flag info's sets
                l_in_formats = [x[const.DB_FORMAT_ID_KEY]
                                for x in self._arg_info[const.DB_IN_FLAGS_FORMATS_KEY_BASE]
                                if [self._key_prefix + const.DB_IN_FLAGS_ID_KEY_BASE] == flag_id]
                l_out_formats = [x[const.DB_FORMAT_ID_KEY]
                                 for x in self._arg_info[const.DB_OUT_FLAGS_FORMATS_KEY_BASE]
                                 if [self._key_prefix + const.DB_OUT_FLAGS_ID_KEY_BASE] == flag_id]
                flag_info.s_in_formats.update(l_in_formats)
                flag_info.s_out_formats.update(l_out_formats)

        return self._d_flag_info

    @property
    def l_flag_info(self) -> list[FlagInfo | None]:
        """Generate the flag info list (indexed by ID) when needed. Returns None if the converter has no flag info in
        the database
        """
        if self._l_flag_info is None and self._key_prefix is not None:
            # Pre-size a list based on the maximum ID plus 1 (since IDs are 1-indexed)
            max_id: int = max([x.id for x in self.d_flag_info.values()])
            self._l_flag_info: list[FlagInfo | None] = [None] * (max_id+1)

            # Fill the list with all flags in the dict
            for single_flag_info in self.d_flag_info.values():
                self._l_flag_info[single_flag_info.id] = single_flag_info

        return self._l_flag_info

    @property
    def d_option_info(self) -> dict[str, OptionInfo] | None:
        """Generate the option info dict (index by flag) when needed. Returns None if the converter has no option info
        in the database
        """
        if self._d_option_info is None and self._key_prefix is not None:
            # Load from the database
            self._d_option_info = {}
            for d_single_option_info in self._arg_info[const.DB_IN_OPTIONS_KEY_BASE]:
                flag: str = d_single_option_info[const.DB_FLAG_KEY]
                option_id = d_single_option_info[const.DB_ID_KEY]
                option_info = OptionInfo(parent=self,
                                         id=option_id,
                                         flag=flag,
                                         brief=d_single_option_info[const.DB_BRIEF_KEY],
                                         description=d_single_option_info[const.DB_DESC_KEY],
                                         info=d_single_option_info[const.DB_INFO_KEY])
                self._d_option_info[flag] = option_info

                # Get a list of all in and formats applicable to this flag, and add them to the flag info's sets
                l_in_formats = [x[const.DB_FORMAT_ID_KEY]
                                for x in self._arg_info[const.DB_IN_OPTIONS_FORMATS_KEY_BASE]
                                if [self._key_prefix + const.DB_IN_OPTIONS_ID_KEY_BASE] == option_id]
                l_out_formats = [x[const.DB_FORMAT_ID_KEY]
                                 for x in self._arg_info[const.DB_OUT_OPTIONS_FORMATS_KEY_BASE]
                                 if [self._key_prefix + const.DB_OUT_OPTIONS_ID_KEY_BASE] == option_id]
                option_info.s_in_formats.update(l_in_formats)
                option_info.s_out_formats.update(l_out_formats)

        return self._d_option_info

    @property
    def l_option_info(self) -> list[OptionInfo | None]:
        """Generate the option info list (indexed by ID) when needed. Returns None if the converter has no option info
        in the database
        """
        if self._l_option_info is None and self._key_prefix is not None:
            # Pre-size a list based on the maximum ID plus 1 (since IDs are 1-indexed)
            max_id: int = max([x.id for x in self.d_option_info.values()])
            self._l_option_info: list[OptionInfo | None] = [None] * (max_id+1)

            # Fill the list with all options in the dict
            for single_option_info in self.d_option_info.values():
                self._l_option_info[single_option_info.id] = single_option_info

        return self._l_option_info

    @property
    def d_arg_info(self) -> dict[str, ArgInfo] | None:
        """Generate the arg info dict (index by flag) when needed, which is a combination of the dicts for flags and
        options. Returns None if the converter has no flag or option info in the database
        """
        if self._d_arg_info is None and self._key_prefix is not None:
            # Make the arg info dict by combining the flag and option info dicts
            self._d_arg_info = copy(self.d_flag_info)
            self._d_arg_info.update(self.d_option_info)
        return self._d_arg_info


class FormatInfo:
    """Class providing information on a file format from the PSDI Data Conversion database
    """

    def __init__(self,
                 name: str,
                 parent: DataConversionDatabase,
                 d_single_format_info: dict[str, bool | int | str | None]):
        """Set up the class - this will be initialised within a `DataConversionDatabase`, which we set as the parent

        Parameters
        ----------
        name : str
            The name (extension) of the file format
        parent : DataConversionDatabase
            The database which this belongs to
        d_single_format_info : dict[str, bool  |  int  |  str  |  None]
            The dict of info on the format stored in the database
        """

        # Load attributes from input
        self.name = name
        self.parent = parent

        # Load attributes from the database
        self.id: int = d_single_format_info.get(const.DB_ID_KEY, -1)
        self.note: str = d_single_format_info.get(const.DB_FORMAT_NOTE_KEY, "")
        self.composition = bool(d_single_format_info.get(const.DB_FORMAT_COMP_KEY, False))
        self.connections = bool(d_single_format_info.get(const.DB_FORMAT_CONN_KEY, False))
        self.two_dim = bool(d_single_format_info.get(const.DB_FORMAT_2D_KEY, False))
        self.three_dim = bool(d_single_format_info.get(const.DB_FORMAT_3D_KEY, False))


class ConversionsTable:
    """Class providing information on available file format conversions.

    Information on internal data handling of this class:

    The idea here is that we need to be able to get information on whether a converter can handle a conversion from one
    file format to another. This results in 3D data storage, with dimensions: Converter, Input Format, Output Format.
    The most important operations are (in roughly descending order of importance):

    - For a given Converter, Input Format, and Output Format, get whether or not the conversion is possible, and the
    degree of success if it is possible.
    - For a given Input Format and Output Format, list available Converters and their degrees of success
    - For a given Converter, list available Input Formats and Output Formats
    - For a given Input Format, list available Output Formats and Converters, and the degree of success of each

    At date of implementation, the data comprises 9 Converters and 280 Input/Output Formats, for 705,600 possibilities,
    increasing linearly with the number of converters and quadratically with the number of formats. (Self-to-self format
    conversions don't need to be stored, but this may not be a useful optimisation.)

    Conversion data is available for 23,013 Converter, Input, Output values, or ~3% of the total possible conversions.
    While this could currently work as a sparse array, it will likely be filled to become denser over time, so a dense
    representation makes the most sense.

    The present implementation uses a list-of-lists-of-lists approach, to avoid adding NumPy as a dependency
    until/unless efficiency concerns motivate it in the future.
    """

    def __init__(self,
                 l_converts_to: list[dict[str, bool | int | str | None]],
                 parent: DataConversionDatabase):
        """Set up the class - this will be initialised within a `DataConversionDatabase`, which we set as the parent

        Parameters
        ----------
        l_converts_to : list[dict[str, bool  |  int  |  str  |  None]]
            The list of dicts in the database providing information on possible conversions
        parent : DataConversionDatabase
            The database which this belongs to

        Raises
        ------
        FileConverterDatabaseException
            _description_
        """

        self.parent = parent

        # Store references to needed data
        self._l_converts_to = l_converts_to

        # To save data, we store the degree of success in the table as an integer indexed to possible strings
        # `self._l_dos` maps from integer index value to degree of success string, and `self._d_dos` maps from degree of
        # success string to index value
        self._l_dos: list[str | None] = [None] + list(set([x[const.DB_SUCCESS_KEY] for x in l_converts_to]))
        self._d_dos = {dos: index for index, dos in enumerate(self._l_dos)}

        # Build the conversion table, indexed Converter, Input Format, Output Format - note that each of these is
        # 1-indexed, so we add 1 to each of the lengths here
        num_converters = len(parent.converters)
        num_formats = len(parent.formats)

        self._table = [[[0 for k in range(num_formats+1)] for j in range(num_formats+1)]
                       for i in range(num_converters+1)]

        for possible_conversion in l_converts_to:

            try:
                conv_id: int = possible_conversion[const.DB_CONV_ID_KEY]
                in_id: int = possible_conversion[const.DB_IN_ID_KEY]
                out_id: int = possible_conversion[const.DB_OUT_ID_KEY]
            except KeyError:
                raise FileConverterDatabaseException(
                    f"Malformed 'converts_to' entry in database: {possible_conversion}")

            self._table[conv_id][in_id][out_id] = self._d_dos[possible_conversion[const.DB_SUCCESS_KEY]]

    def get_degree_of_success(self,
                              converter_name: str,
                              in_format: str,
                              out_format: str) -> str | None:
        """Get the degree of success for a desired conversion, represented as a string (or else None if not possible)

        Parameters
        ----------
        converter_name : str
            The name of the converter to use
        in_format : str
            The extension of the input file format
        out_format : str
            The extension of the output file format

        Returns
        -------
        str | None
            If the conversion is possible, returns a string describing the degree of success. If the conversion is not
            possible, returns None
        """

        conv_id: int = self.parent.get_converter_info(converter_name).id
        in_id: int = self.parent.get_format_info(in_format).id
        out_id: int = self.parent.get_format_info(out_format).id

        return self._l_dos[self._table[conv_id][in_id][out_id]]

    def get_possible_converters(self,
                                in_format: str,
                                out_format: str) -> list[tuple[str, str]]:
        """Get a list of converters which can perform a conversion from one format to another and the degree of success
        with each of these converters

        Parameters
        ----------
        in_format : str
            The extension of the input file format
        out_format : str
            The extension of the output file format

        Returns
        -------
        list[tuple[str, str]]
            A list of tuples, where each tuple's first item is the name of a converter which can perform this
            conversion, and the second item is the degree of success for the conversion
        """
        in_id: int = self.parent.get_format_info(in_format).id
        out_id: int = self.parent.get_format_info(out_format).id

        # Slice the table to get a list of the success for this conversion for each converter
        l_converter_success = [x[in_id][out_id] for x in self._table]

        # Filter for possible conversions (dos_index > 0) and get the converter name and degree-of-success string
        # for each possible conversion
        l_possible_converters_and_dos = [(self.parent.get_converter_info(converter_id).name,
                                          self._l_dos[dos_index]) for converter_id, dos_index
                                         in enumerate(l_converter_success) if dos_index > 0]

        return l_possible_converters_and_dos

    def get_possible_formats(self, converter_name: str) -> tuple[list[str], list[str]]:
        """Get a list of input and output formats that a given converter supports

        Parameters
        ----------
        converter_name : str
            The name of the converter

        Returns
        -------
        tuple[list[str], list[str]]
            A tuple of a list of the supported input formats and a list of the supported output formats
        """
        conv_id: int = self.parent.get_converter_info(converter_name).id
        ll_in_out_format_dos = self._table[conv_id]

        # Filter for possible input formats by checking if at least one output format for each has a degree of success
        # index greater than 0, and stored the filtered lists where the input format is possible so we only need to
        # check them for possible output formats
        (l_possible_in_format_ids,
         ll_filtered_in_out_format_dos) = zip(*[(i, l_out_format_dos) for i, l_out_format_dos
                                                in enumerate(ll_in_out_format_dos) if sum(l_out_format_dos) > 0])

        # As with input IDs, filter for output IDs where at least one input format has a degree of success index greater
        # than 0. A bit more complicated for the second index, forcing us to do list comprehension to fetch a list
        # across the table before summing
        l_possible_out_format_ids = [j for j, _ in enumerate(ll_filtered_in_out_format_dos[0]) if
                                     sum([x[j] for x in ll_filtered_in_out_format_dos]) > 0]

        # Get the name for each format ID, and return lists of the names
        return ([self.parent.get_format_info(x).name for x in l_possible_in_format_ids],
                [self.parent.get_format_info(x).name for x in l_possible_out_format_ids])


class DataConversionDatabase:
    """Class providing interface for information contained in the PSDI Data Conversion database
    """

    def __init__(self, d_data: dict[str, Any]):
        """Initialise the DataConversionDatabase object

        Parameters
        ----------
        d_data : dict[str, Any]
            The dict of the database, as loaded in from the JSON file
        """

        # Store the database dict internally for debugging purposes
        self._d_data = d_data

        # Store top-level items not tied to a specific converter
        self.formats: list[dict[str, bool | int | str | None]] = d_data[const.DB_FORMATS_KEY]
        self.converters: list[dict[str, bool | int | str | None]] = d_data[const.DB_CONVERTERS_KEY]
        self.converts_to: list[dict[str, bool | int | str | None]] = d_data[const.DB_CONVERTS_TO_KEY]

        # Placeholders for properties that are generated when needed
        self._d_converter_info: dict[str, ConverterInfo] | None = None
        self._l_converter_info: list[ConverterInfo] | None = None
        self._d_format_info: dict[str, FormatInfo] | None = None
        self._l_format_info: list[FormatInfo] | None = None
        self._conversions_table: ConversionsTable | None = None

    @property
    def d_converter_info(self) -> dict[str, ConverterInfo]:
        """Generate the converter info dict (indexed by name) when needed
        """
        if self._d_converter_info is None:
            self._d_converter_info: dict[str, ConverterInfo] = {}
            for d_single_converter_info in self.converters:
                name: str = d_single_converter_info[const.DB_NAME_KEY]
                if name in self._d_converter_info:
                    logger.warning(f"Converter '{name}' appears more than once in the database. Only the first instance"
                                   " will be used.")
                    continue

                self._d_converter_info[name] = ConverterInfo(name=name,
                                                             parent=self,
                                                             d_single_converter_info=d_single_converter_info,
                                                             d_data=self._d_data)
        return self._d_converter_info

    @property
    def l_converter_info(self) -> list[ConverterInfo | None]:
        """Generate the converter info list (indexed by ID) when needed
        """
        if self._l_converter_info is None:
            # Pre-size a list based on the maximum ID plus 1 (since IDs are 1-indexed)
            max_id: int = max([x[const.DB_ID_KEY] for x in self.converters])
            self._l_converter_info: list[ConverterInfo | None] = [None] * (max_id+1)

            # Fill the list with all converters in the dict
            for single_converter_info in self.d_converter_info.values():
                self._l_converter_info[single_converter_info.id] = single_converter_info

        return self._l_converter_info

    @property
    def d_format_info(self) -> dict[str, FormatInfo]:
        """Generate the format info dict when needed
        """
        if self._d_format_info is None:
            self._d_format_info: dict[str, FormatInfo] = {}

            for d_single_format_info in self.formats:
                name: str = d_single_format_info[const.DB_FORMAT_EXT_KEY]

                format_info = FormatInfo(name=name,
                                         parent=self,
                                         d_single_format_info=d_single_format_info)

                if name in self._d_format_info:
                    logger.debug(f"File extension '{name}' appears more than once in the database. Duplicates will use "
                                 "a key appended with an index")
                    loop_concluded = False
                    for i in range(97):
                        test_name = f"{name}-{i+2}"
                        if test_name in self._d_format_info:
                            continue
                        else:
                            self._d_format_info[test_name] = format_info
                            loop_concluded = True
                            break
                    if not loop_concluded:
                        logger.warning("Loop counter exceeded when searching for valid new name for file extension "
                                       f"'{name}'. New entry will not be added to the database to avoid possibility of "
                                       "an infinite loop")
                else:
                    self._d_format_info[name] = format_info
        return self._d_format_info

    @property
    def l_format_info(self) -> list[FormatInfo | None]:
        """Generate the format info list (indexed by ID) when needed
        """
        if self._l_format_info is None:
            # Pre-size a list based on the maximum ID plus 1 (since IDs are 1-indexed)
            max_id: int = max([x[const.DB_ID_KEY] for x in self.formats])
            self._l_format_info: list[FormatInfo | None] = [None] * (max_id+1)

            # Fill the list with all formats in the dict
            for single_format_info in self.d_format_info.values():
                self._l_format_info[single_format_info.id] = single_format_info

        return self._l_format_info

    @property
    def conversions_table(self) -> ConversionsTable:
        """Generates the conversions table when needed
        """
        if self._conversions_table is None:
            self._conversions_table = ConversionsTable(l_converts_to=self.converts_to,
                                                       parent=self)
        return self._conversions_table

    def get_converter_info(self, converter_name_or_id: str | int) -> ConverterInfo:
        """Get a converter's info from either its name or ID
        """
        if isinstance(converter_name_or_id, str):
            try:
                return self.d_converter_info[converter_name_or_id]
            except KeyError:
                raise FileConverterDatabaseException(f"Converter name '{converter_name_or_id}' not recognised")
        elif isinstance(converter_name_or_id, int):
            return self.l_converter_info[converter_name_or_id]
        else:
            raise FileConverterDatabaseException(f"Invalid key passed to `get_converter_info`: '{converter_name_or_id}'"
                                                 f" of type '{type(converter_name_or_id)}'. Type must be `str` or "
                                                 "`int`")

    def get_format_info(self, format_name_or_id: str | int) -> FormatInfo:
        """Get a format's ID info from either its name or ID
        """
        if isinstance(format_name_or_id, str):
            try:
                return self.d_format_info[format_name_or_id]
            except KeyError:
                raise FileConverterDatabaseException(f"Format name '{format_name_or_id}' not recognised")
        elif isinstance(format_name_or_id, int):
            return self.l_format_info[format_name_or_id]
        else:
            raise FileConverterDatabaseException(f"Invalid key passed to `get_format_info`: '{format_name_or_id}'"
                                                 f" of type '{type(format_name_or_id)}'. Type must be `str` or "
                                                 "`int`")


# The database will be loaded on demand when `get_database()` is called
_database: DataConversionDatabase | None = None


def load_database() -> DataConversionDatabase:
    """Load and return a new instance of the data conversion database from the JSON database file in this package. This
    function should not be called directly unless you specifically need a new instance of the database object and can't
    deepcopy the database returned by `get_database()`, as it's expensive to load it in.

    Returns
    -------
    DataConversionDatabase
    """

    # Find and load the database JSON file

    # For an interactive shell, __file__ won't be defined for this module, so use the constants module instead
    reference_file = os.path.realpath(const.__file__)

    qualified_database_filename = os.path.join(os.path.dirname(reference_file), const.DATABASE_FILENAME)
    d_data: dict = json.load(open(qualified_database_filename, "r"))

    return DataConversionDatabase(d_data)


def get_database() -> DataConversionDatabase:
    """Gets the global database object, loading it in first if necessary. Since it's computationally expensive to load
    the database, it's best treated as an immutable singleton.

    Returns
    -------
    DataConversionDatabase
        The global database object
    """
    global _database
    if _database is None:
        # Create the database object and store it globally
        _database = load_database()
    return _database


def get_converter_info(name: str) -> ConverterInfo:
    """Gets the information on a given converter stored in the database

    Parameters
    ----------
    name : str
        The name of the converter

    Returns
    -------
    ConverterInfo
    """

    return get_database().d_converter_info[name]


def get_format_info(name: str) -> FormatInfo:
    """Gets the information on a given file format stored in the database

    Parameters
    ----------
    name : str
        The name (extension) of the form

    Returns
    -------
    FormatInfo
    """

    return get_database().d_format_info[name]


def get_degree_of_success(converter_name: str,
                          in_format: str,
                          out_format: str) -> str | None:
    """Get the degree of success for a desired conversion, represented as a string (or else None if not possible)

    Parameters
    ----------
    converter_name : str
        The name of the converter to use
    in_format : str
        The extension of the input file format
    out_format : str
        The extension of the output file format

    Returns
    -------
    str | None
        If the conversion is possible, returns a string describing the degree of success. If the conversion is not
        possible, returns None
    """

    return get_database().conversions_table.get_degree_of_success(converter_name=converter_name,
                                                                  in_format=in_format,
                                                                  out_format=out_format)


def get_possible_converters(in_format: str,
                            out_format: str) -> list[tuple[str, str]]:
    """Get a list of converters which can perform a conversion from one format to another and the degree of success
    with each of these converters

    Parameters
    ----------
    in_format : str
        The extension of the input file format
    out_format : str
        The extension of the output file format

    Returns
    -------
    list[tuple[str, str]]
        A list of tuples, where each tuple's first item is the name of a converter which can perform this
        conversion, and the second item is the degree of success for the conversion
    """

    return get_database().conversions_table.get_possible_converters(in_format=in_format,
                                                                    out_format=out_format)


def get_possible_formats(converter_name: str) -> tuple[list[str], list[str]]:
    """Get a list of input and output formats that a given converter supports

    Parameters
    ----------
    converter_name : str
        The name of the converter

    Returns
    -------
    tuple[list[str], list[str]]
        A tuple of a list of the supported input formats and a list of the supported output formats
    """
    return get_database().conversions_table.get_possible_formats(converter_name=converter_name)


data = get_database()
ob = get_converter_info('Open Babel')
