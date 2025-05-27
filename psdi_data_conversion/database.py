"""@file psdi_data_conversion/database.py

Created 2025-02-03 by Bryan Gillis.

Python module provide utilities for accessing the converter database
"""

from __future__ import annotations

import json
import os
import sys
from copy import copy
from dataclasses import dataclass, field
from functools import lru_cache
from itertools import product
from logging import getLogger
from typing import Any, Literal, overload
from warnings import catch_warnings

import igraph as ig

from psdi_data_conversion import constants as const
from psdi_data_conversion.converter import (L_REGISTERED_CONVERTERS, L_SUPPORTED_CONVERTERS,
                                            get_registered_converter_class)
from psdi_data_conversion.converters.base import FileConverter, FileConverterException
from psdi_data_conversion.utils import regularize_name

# Keys for top-level and general items in the database
DB_FORMATS_KEY = "formats"
DB_CONVERTERS_KEY = "converters"
DB_CONVERTS_TO_KEY = "converts_to"
DB_ID_KEY = "id"
DB_NAME_KEY = "name"

# Keys for converter general info in the database
DB_DESC_KEY = "description"
DB_INFO_KEY = "further_info"
DB_URL_KEY = "url"

# Keys for format general info in the database
DB_FORMAT_EXT_KEY = "extension"
DB_FORMAT_C2X_KEY = "format"
DB_FORMAT_NOTE_KEY = "note"
DB_FORMAT_COMP_KEY = "composition"
DB_FORMAT_CONN_KEY = "connections"
DB_FORMAT_2D_KEY = "two_dim"
DB_FORMAT_3D_KEY = "three_dim"

# Keys for converts_to info in the database
DB_CONV_ID_KEY = "converters_id"
DB_IN_ID_KEY = "in_id"
DB_OUT_ID_KEY = "out_id"
DB_SUCCESS_KEY = "degree_of_success"

# Key bases for converter-specific items in the database
DB_IN_FLAGS_KEY_BASE = "flags_in"
DB_OUT_FLAGS_KEY_BASE = "flags_out"
DB_IN_OPTIONS_KEY_BASE = "argflags_in"
DB_OUT_OPTIONS_KEY_BASE = "argflags_out"
DB_IN_FLAGS_FORMATS_KEY_BASE = "format_to_flags_in"
DB_OUT_FLAGS_FORMATS_KEY_BASE = "format_to_flags_out"
DB_IN_OPTIONS_FORMATS_KEY_BASE = "format_to_argflags_in"
DB_OUT_OPTIONS_FORMATS_KEY_BASE = "format_to_argflags_out"

# Keys for argument info in the database
DB_FLAG_KEY = "flag"
DB_BRIEF_KEY = "brief"
DB_FORMAT_ID_KEY = "formats_id"
DB_IN_FLAGS_ID_KEY_BASE = "flags_in_id"
DB_OUT_FLAGS_ID_KEY_BASE = "flags_out_id"
DB_IN_OPTIONS_ID_KEY_BASE = "argflags_in_id"
DB_OUT_OPTIONS_ID_KEY_BASE = "argflags_out_id"

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
    # We need to provide a default argument here, since it will come after the sets with default arguments in ArgInfo
    brief: str = ""


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
            The regularized name of the converter
        parent : DataConversionDatabase
            The database which this belongs to
        d_single_converter_info : dict[str, int | str]
            The dict within the database file which describes this converter
        d_data : dict[str, Any]
            The loaded database dict
        """

        self.name = regularize_name(name)
        """The regularized name of the converter"""

        self.converter_class: type[FileConverter]
        """The class used to perform conversions with this converter"""

        self.pretty_name: str
        """The name of the converter, properly spaced and capitalized"""

        try:
            self.converter_class = get_registered_converter_class(self.name)
            self.pretty_name = self.converter_class.name
        except KeyError:
            self.converter_class = None
            self.pretty_name = name

        self.parent = parent
        """The parent database"""

        # Get info about the converter from the database
        self.id: int = d_single_converter_info.get(DB_ID_KEY, -1)
        """The converter's ID"""

        self.description: str = d_single_converter_info.get(DB_DESC_KEY, "")
        """A description of the converter"""

        self.url: str = d_single_converter_info.get(DB_URL_KEY, "")
        """The official URL for the converter"""

        # Get necessary info about the converter from the class
        try:
            self._key_prefix = get_registered_converter_class(name).database_key_prefix
        except KeyError:
            # We'll get a KeyError for converters in the database that don't yet have their own class, which we can
            # safely ignore
            self._key_prefix = None

        self._arg_info: dict[str, list[dict[str, int | str]]] = {}

        # Placeholders for members that are generated when needed
        self._l_in_flag_info: list[FlagInfo] | None = None
        self._l_out_flag_info: list[FlagInfo] | None = None
        self._l_in_option_info: list[OptionInfo] | None = None
        self._l_out_option_info: list[OptionInfo] | None = None

        self._d_in_format_flags: dict[str | int, set[str]] | None = None
        self._d_out_format_flags: dict[str | int, set[str]] | None = None
        self._d_in_format_options: dict[str | int, set[str]] | None = None
        self._d_out_format_options: dict[str | int, set[str]] | None = None

        # If the converter class has no defined key prefix, don't add any extra info for it
        if self._key_prefix is None:
            return
        for key_base in (DB_IN_FLAGS_KEY_BASE,
                         DB_OUT_FLAGS_KEY_BASE,
                         DB_IN_OPTIONS_KEY_BASE,
                         DB_OUT_OPTIONS_KEY_BASE,
                         DB_IN_FLAGS_FORMATS_KEY_BASE,
                         DB_OUT_FLAGS_FORMATS_KEY_BASE,
                         DB_IN_OPTIONS_FORMATS_KEY_BASE,
                         DB_OUT_OPTIONS_FORMATS_KEY_BASE):
            self._arg_info[key_base] = d_data.get(self._key_prefix + key_base)

    def _create_l_arg_info(self, subclass: type[ArgInfo]) -> tuple[list[ArgInfo], list[ArgInfo]]:
        """Creates either the flag or option info list
        """

        # Set values based on whether we're working with flags or options
        if issubclass(subclass, FlagInfo):
            in_key_base = DB_IN_FLAGS_KEY_BASE
            out_key_base = DB_OUT_FLAGS_KEY_BASE
            in_formats_key_base = DB_IN_FLAGS_FORMATS_KEY_BASE
            in_args_id_key_base = DB_IN_FLAGS_ID_KEY_BASE
            out_formats_key_base = DB_OUT_FLAGS_FORMATS_KEY_BASE
            out_args_id_key_base = DB_OUT_FLAGS_ID_KEY_BASE
        elif issubclass(subclass, OptionInfo):
            in_key_base = DB_IN_OPTIONS_KEY_BASE
            out_key_base = DB_OUT_OPTIONS_KEY_BASE
            in_formats_key_base = DB_IN_OPTIONS_FORMATS_KEY_BASE
            in_args_id_key_base = DB_IN_OPTIONS_ID_KEY_BASE
            out_formats_key_base = DB_OUT_OPTIONS_FORMATS_KEY_BASE
            out_args_id_key_base = DB_OUT_OPTIONS_ID_KEY_BASE
        else:
            raise FileConverterDatabaseException(f"Unrecognised subclass passed to `_create_l_arg_info`: {subclass}")

        for key_base, in_or_out in ((in_key_base, "in"),
                                    (out_key_base, "out")):

            max_id = max([x[DB_ID_KEY] for x in self._arg_info[key_base]])
            l_arg_info: list[ArgInfo] = [None]*(max_id+1)

            for d_single_arg_info in self._arg_info[key_base]:
                name: str = d_single_arg_info[DB_FLAG_KEY]
                arg_id: int = d_single_arg_info[DB_ID_KEY]
                brief = d_single_arg_info.get(DB_BRIEF_KEY)
                optional_arg_info_kwargs = {}
                if brief is not None:
                    optional_arg_info_kwargs["brief"] = brief
                arg_info = subclass(parent=self,
                                    id=arg_id,
                                    flag=name,
                                    description=d_single_arg_info[DB_DESC_KEY],
                                    info=d_single_arg_info[DB_INFO_KEY],
                                    **optional_arg_info_kwargs)
                l_arg_info[arg_id] = arg_info

                # Get a list of all in and formats applicable to this flag, and add them to the flag info's sets
                if in_or_out == "in":
                    l_in_formats = [x[DB_FORMAT_ID_KEY]
                                    for x in self._arg_info[in_formats_key_base]
                                    if x[self._key_prefix + in_args_id_key_base] == arg_id]
                    arg_info.s_in_formats.update(l_in_formats)
                else:
                    l_out_formats = [x[DB_FORMAT_ID_KEY]
                                     for x in self._arg_info[out_formats_key_base]
                                     if x[self._key_prefix + out_args_id_key_base] == arg_id]
                    arg_info.s_out_formats.update(l_out_formats)

            if in_or_out == "in":
                l_in_arg_info = l_arg_info
            else:
                l_out_arg_info = l_arg_info

        return l_in_arg_info, l_out_arg_info

    @property
    def l_in_flag_info(self) -> list[FlagInfo | None]:
        """Generate the input flag info list (indexed by ID) when needed. Returns None if the converter has no flag info
        in the database
        """
        if self._l_in_flag_info is None and self._key_prefix is not None:
            self._l_in_flag_info, self._l_out_flag_info = self._create_l_arg_info(FlagInfo)
        return self._l_in_flag_info

    @property
    def l_out_flag_info(self) -> list[FlagInfo | None]:
        """Generate the output flag info list (indexed by ID) when needed. Returns None if the converter has no flag
        info in the database
        """
        if self._l_out_flag_info is None and self._key_prefix is not None:
            self._l_in_flag_info, self._l_out_flag_info = self._create_l_arg_info(FlagInfo)
        return self._l_out_flag_info

    @property
    def l_in_option_info(self) -> list[OptionInfo | None]:
        """Generate the input option info list (indexed by ID) when needed. Returns None if the converter has no option
        info in the database
        """
        if self._l_in_option_info is None and self._key_prefix is not None:
            self._l_in_option_info, self._l_out_option_info = self._create_l_arg_info(OptionInfo)
        return self._l_in_option_info

    @property
    def l_out_option_info(self) -> list[OptionInfo | None]:
        """Generate the output option info list (indexed by ID) when needed. Returns None if the converter has no option
        info in the database
        """
        if self._l_out_option_info is None and self._key_prefix is not None:
            self._l_in_option_info, self._l_out_option_info = self._create_l_arg_info(OptionInfo)
        return self._l_out_option_info

    def _create_d_format_args(self,
                              subclass: type[ArgInfo],
                              in_or_out: str) -> dict[str | int, set[int]]:
        """Creates either the flag or option format args dict
        """

        if in_or_out not in ("in", "out"):
            raise FileConverterDatabaseException(
                f"Unrecognised `in_or_out` value passed to `_create_d_format_args`: {in_or_out}")

        # Set values based on whether we're working with flags or options, and input or output
        if issubclass(subclass, FlagInfo):
            l_arg_info = self.l_in_flag_info if in_or_out == "in" else self.l_out_flag_info
        elif issubclass(subclass, OptionInfo):
            l_arg_info = self.l_in_option_info if in_or_out == "in" else self.l_out_option_info
        else:
            raise FileConverterDatabaseException(
                f"Unrecognised subclass passed to `_create_d_format_args`: {subclass}")

        d_format_args: dict[str | int, set[ArgInfo]] = {}
        l_parent_format_info = self.parent.l_format_info

        # If the converter doesn't provide argument info, set l_arg_info to an empty list so it can be iterated in
        # the next step, rather than None
        if not l_arg_info:
            l_arg_info = []

        for arg_info in l_arg_info:

            if arg_info is None:
                continue

            if in_or_out == "in":
                s_formats = arg_info.s_in_formats
            else:
                s_formats = arg_info.s_out_formats
            l_format_info = [l_parent_format_info[format_id] for format_id in s_formats]
            for format_info in l_format_info:
                format_name = format_info.name
                format_id = format_info.id

                # Add an empty set for this format to the dict if it isn't yet there, otherwise add to the set
                if format_name not in d_format_args:
                    d_format_args[format_name] = set()
                    # Keying by ID will point to the same set as keying by name
                    d_format_args[format_id] = d_format_args[format_name]

                d_format_args[format_name].add(arg_info.id)

        return d_format_args

    @property
    def d_in_format_flags(self) -> dict[str | int, set[int]]:
        """Generate the dict of flags for an input format (keyed by format name/extension or format ID) when needed.
        The format will not be in the dict if no flags are accepted
        """
        if self._d_in_format_flags is None:
            self._d_in_format_flags = self._create_d_format_args(FlagInfo, "in")
        return self._d_in_format_flags

    @property
    def d_out_format_flags(self) -> dict[str | int, set[int]]:
        """Generate the dict of flags for an output format (keyed by format name/extension or format ID) when needed.
        The format will not be in the dict if no options are accepted
        """
        if self._d_out_format_flags is None:
            self._d_out_format_flags = self._create_d_format_args(FlagInfo, "out")
        return self._d_out_format_flags

    @property
    def d_in_format_options(self) -> dict[str | int, set[int]]:
        """Generate the dict of options for an input format (keyed by format name/extension or format ID) when needed.
        The format will not be in the dict if no options are accepted
        """
        if self._d_in_format_options is None:
            self._d_in_format_options = self._create_d_format_args(OptionInfo, "in")
        return self._d_in_format_options

    @property
    def d_out_format_options(self) -> dict[str | int, set[int]]:
        """Generate the dict of options for an output format (keyed by format name/extension or format ID) when needed.
        The format will not be in the dict if no options are accepted
        """
        if self._d_out_format_options is None:
            self._d_out_format_options = self._create_d_format_args(OptionInfo, "out")
        return self._d_out_format_options

    def get_in_format_args(self, name: str) -> tuple[list[FlagInfo], list[OptionInfo]]:
        """Get the input flags and options supported for a given format (provided as its extension)

        Parameters
        ----------
        name : str
            The file format name (extension)

        Returns
        -------
        tuple[set[FlagInfo], set[OptionInfo]]
            A set of info for the allowed flags, and a set of info for the allowed options
        """
        l_flag_ids = list(self.d_in_format_flags.get(name, set()))
        l_flag_ids.sort()
        l_flag_info = [self.l_in_flag_info[x] for x in l_flag_ids]

        l_option_ids = list(self.d_in_format_options.get(name, set()))
        l_option_ids.sort()
        l_option_info = [self.l_in_option_info[x] for x in l_option_ids]

        return l_flag_info, l_option_info

    def get_out_format_args(self, name: str) -> tuple[list[FlagInfo], list[OptionInfo]]:
        """Get the output flags and options supported for a given format (provided as its extension)

        Parameters
        ----------
        name : str
            The file format name (extension)

        Returns
        -------
        tuple[set[FlagInfo], set[OptionInfo]]
            A set of info for the allowed flags, and a set of info for the allowed options
        """
        l_flag_ids = list(self.d_out_format_flags.get(name, set()))
        l_flag_ids.sort()
        l_flag_info = [self.l_out_flag_info[x] for x in l_flag_ids]

        l_option_ids = list(self.d_out_format_options.get(name, set()))
        l_option_ids.sort()
        l_option_info = [self.l_out_option_info[x] for x in l_option_ids]

        return l_flag_info, l_option_info


class FormatInfo:
    """Class providing information on a file format from the PSDI Data Conversion database
    """

    D_PROPERTY_ATTRS = {const.QUAL_COMP_KEY: const.QUAL_COMP_LABEL,
                        const.QUAL_CONN_KEY: const.QUAL_CONN_LABEL,
                        const.QUAL_2D_KEY: const.QUAL_2D_LABEL,
                        const.QUAL_3D_KEY: const.QUAL_3D_LABEL}
    """A dict of attrs of this class which describe properties that a format may or may not have"""

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
        """The name of this format"""

        self.parent = parent
        """The database which this format belongs to"""

        # Load attributes from the database
        self.id: int = d_single_format_info.get(DB_ID_KEY, -1)
        """The ID of this format"""

        self.c2x_format: str = d_single_format_info.get(DB_FORMAT_C2X_KEY)
        """The name of this format as the c2x converter expects it"""

        self.note: str = d_single_format_info.get(DB_FORMAT_NOTE_KEY, "")
        """The description of this format"""

        self.composition = d_single_format_info.get(DB_FORMAT_COMP_KEY)
        """Whether or not this format stores composition information"""

        self.connections = d_single_format_info.get(DB_FORMAT_CONN_KEY)
        """Whether or not this format stores connections information"""

        self.two_dim = d_single_format_info.get(DB_FORMAT_2D_KEY)
        """Whether or not this format stores 2D structural information"""

        self.three_dim = d_single_format_info.get(DB_FORMAT_3D_KEY)
        """Whether or not this format stores 3D structural information"""

        self._lower_name: str = self.name.lower()
        """The format name all in lower-case"""

        self._disambiguated_name: str | None = None

    @property
    def disambiguated_name(self) -> str:
        """A unique name for this format which can be used to distinguish it from others which share the same extension,
        by appending the name of each with a unique index"""
        if self._disambiguated_name is None:
            l_formats_with_same_name = [x for x in self.parent.l_format_info
                                        if x and x._lower_name == self._lower_name]
            if len(l_formats_with_same_name) == 1:
                self._disambiguated_name = self._lower_name
            else:
                index_of_this = [i for i, x in enumerate(l_formats_with_same_name) if self is x][0]
                self._disambiguated_name = f"{self._lower_name}-{index_of_this}"
        return self._disambiguated_name

    def __str__(self):
        """When cast to string, convert to the name (extension) of the format"""
        return self.name

    def __int__(self):
        """When cast to int, return the ID of the format"""
        return self.id


@dataclass
class PropertyConversionInfo:
    """Class representing whether a given property is present in the input and/out output file formats, and a note on
    what its presence or absence means
    """
    key: str
    input_supported: bool | None
    output_supported: bool | None
    label: str = field(init=False)
    note: str = field(init=False)

    def __post_init__(self):
        """Set the label and note based on input/output status
        """
        self.label = FormatInfo.D_PROPERTY_ATTRS[self.key]

        if self.input_supported is None and self.output_supported is None:
            self.note = const.QUAL_NOTE_BOTH_UNKNOWN
        elif self.input_supported is None and self.output_supported is not None:
            self.note = const.QUAL_NOTE_IN_UNKNOWN
        elif self.input_supported is not None and self.output_supported is None:
            self.note = const.QUAL_NOTE_OUT_UNKNOWN
        elif self.input_supported == self.output_supported:
            self.note = ""
        elif self.input_supported:
            self.note = const.QUAL_NOTE_OUT_MISSING
        else:
            self.note = const.QUAL_NOTE_IN_MISSING

        if self.note:
            self.note = self.note.format(self.label)


@dataclass
class ConversionQualityInfo:
    """Class describing the quality of a conversion from one format to another with a given converter.
    """

    converter_name: str
    """The name of the converter"""

    in_format: str
    """The extension of the input file format"""

    out_format: str
    """The extension of the output file format"""

    qual_str: str
    """A string describing the quality of the conversion"""

    details: str
    """A string providing details on any possible issues with the conversion"""

    d_prop_conversion_info: dict[str, PropertyConversionInfo]
    """A dict of PropertyConversionInfo objects, which provide information on each property's support in the
    input and output file formats and a note on the implications
    """

    def __post_init__(self):
        """Regularize the converter name"""
        self.converter_name = regularize_name(self.converter_name)


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
        """

        self.parent = parent

        # Store references to needed data
        self._l_converts_to = l_converts_to

        # Build the conversion graphs - each format is a vertex, each conversion is an edge
        num_formats = len(parent.formats)

        l_supported_conversions = [x for x in l_converts_to if
                                   self.parent.get_converter_info(x[DB_CONV_ID_KEY]).name in L_SUPPORTED_CONVERTERS]
        l_registered_conversions = [x for x in l_supported_conversions if
                                    self.parent.get_converter_info(x[DB_CONV_ID_KEY]).name in L_REGISTERED_CONVERTERS]

        # We make separate graphs for all known conversions, all supported conversions, and all registered conversions
        self.graph: ig.Graph
        self.supported_graph: ig.Graph
        self.registered_graph: ig.Graph

        for support_type, l_conversions in (("", l_converts_to),
                                            ("supported_", l_supported_conversions),
                                            ("registered_", l_registered_conversions)):

            setattr(self, support_type+"graph",
                    ig.Graph(n=num_formats,
                             directed=True,
                             # Each vertex stores the disambiguated name of the format
                             vertex_attrs={DB_NAME_KEY: [x.disambiguated_name if x is not None else None
                                                         for x in parent.l_format_info]},
                             edges=[(x[DB_IN_ID_KEY], x[DB_OUT_ID_KEY]) for x in l_conversions],
                             # Each edge stores the id and name of the converter used for the conversion
                             edge_attrs={DB_CONV_ID_KEY: [x[DB_CONV_ID_KEY] for x in l_conversions],
                                         DB_NAME_KEY: [self.parent.get_converter_info(x[DB_CONV_ID_KEY]).name
                                                       for x in l_conversions]}))

    def _get_desired_graph(self,
                           only: Literal["all"] | Literal["supported"] | Literal["registered"] = "all") -> ig.Graph:
        if only == "all":
            return self.graph
        elif only == "supported":
            return self.supported_graph
        elif only == "registered":
            return self.registered_graph
        else:
            raise ValueError(f"Invalid value \"{only}\" for keyword argument `only`. Allowed values are \"all\" "
                             "(default), \"supported\", and \"registered\".")

    def _get_possible_converters(self, in_format_info: FormatInfo, out_format_info: FormatInfo,
                                 only: Literal["all"] | Literal["supported"] | Literal["registered"] = "all"):
        """Get a list of all converters which can convert from one format to another
        """
        graph = self._get_desired_graph(only)
        l_edges = graph.es.select(_source=in_format_info.id, _target=out_format_info.id)
        return [x[DB_NAME_KEY] for x in l_edges]

    @lru_cache(maxsize=None)
    def get_conversion_quality(self,
                               converter_name: str,
                               in_format: str | int,
                               out_format: str | int) -> ConversionQualityInfo | None:
        """Get an indication of the quality of a conversion from one format to another, or if it's not possible

        Parameters
        ----------
        converter_name : str
            The name of the converter
        in_format : str | int
            The extension or ID of the input file format
        out_format : str | int
            The extension or ID of the output file format

        Returns
        -------
        ConversionQualityInfo | None
            If the conversion is not possible, returns None. If the conversion is possible, returns a
            `ConversionQualityInfo` object with info on the conversion
        """

        # Check if this converter deals with ambiguous formats, so we know if we need to be strict about getting format
        # info
        if get_registered_converter_class(converter_name).supports_ambiguous_extensions:
            which_format = None
        else:
            which_format = 0

        # Get the full format info for each format
        in_format_info = self.parent.get_format_info(in_format, which_format)
        out_format_info: int = self.parent.get_format_info(out_format, which_format)

        # First check if the conversion is possible
        if converter_name not in self._get_possible_converters(in_format_info, out_format_info):
            return None

        # The conversion is possible. Now determine how many properties of the output format are not in the input
        # format and might end up being extrapolated
        num_out_props = 0
        num_new_props = 0
        any_unknown = False
        d_prop_conversion_info: dict[str, PropertyConversionInfo] = {}
        for prop in FormatInfo.D_PROPERTY_ATTRS:
            in_prop: bool | None = getattr(in_format_info, prop)
            out_prop: bool | None = getattr(out_format_info, prop)

            d_prop_conversion_info[prop] = PropertyConversionInfo(prop, in_prop, out_prop)

            # Check for None, indicating we don't have full information on both formats
            if in_prop is None or out_prop is None:
                any_unknown = True
            elif out_prop:
                num_out_props += 1
                if not in_prop:
                    num_new_props += 1

        # Determine the conversion quality
        if num_out_props > 0:
            qual_ratio = 1 - num_new_props/num_out_props
        else:
            qual_ratio = 1

        if any_unknown:
            qual_str = const.QUAL_UNKNOWN
        elif num_out_props == 0 or qual_ratio >= 0.8:
            qual_str = const.QUAL_VERYGOOD
        elif qual_ratio >= 0.6:
            qual_str = const.QUAL_GOOD
        elif qual_ratio >= 0.4:
            qual_str = const.QUAL_OKAY
        elif qual_ratio >= 0.2:
            qual_str = const.QUAL_POOR
        else:
            qual_str = const.QUAL_VERYPOOR

        # Construct the details string for info on possible issues with the conversion

        # Sort the keys by label alphabetically
        l_props: list[str] = list(d_prop_conversion_info.keys())
        l_props.sort(key=lambda x: d_prop_conversion_info[x].label)

        details = "\n".join([d_prop_conversion_info[x].note for x in l_props if d_prop_conversion_info[x].note])

        return ConversionQualityInfo(converter_name=converter_name,
                                     in_format=in_format,
                                     out_format=out_format,
                                     qual_str=qual_str,
                                     details=details,
                                     d_prop_conversion_info=d_prop_conversion_info)

    def get_possible_conversions(self,
                                 in_format: str | int,
                                 out_format: str | int,
                                 only: Literal["all"] | Literal["supported"] | Literal["registered"] = "all"
                                 ) -> list[tuple[ConverterInfo, FormatInfo, FormatInfo]]:
        """Get a list of converters which can perform a conversion from one format to another, disambiguating in the
        case of ambiguous formats and providing IDs for input/output formats for possible conversions

        Parameters
        ----------
        in_format : str | int
            The extension or ID of the input file format
        out_format : str | int
            The extension or ID of the output file format

        Returns
        -------
        list[tuple[ConverterInfo, FormatInfo, FormatInfo]]
            A list of tuples, where each tuple's first item is the ConverterInfo of a converter which can perform a
            matching conversion, the second is the info of the input format for this conversion, and the third is the
            info of the output format
        """
        l_in_format_infos = self.parent.get_format_info(in_format, which="all")
        l_out_format_infos = self.parent.get_format_info(out_format, which="all")

        # Start a list of all possible conversions
        l_possible_conversions = []

        # Iterate over all possible combinations of input and output formats
        for in_format_info, out_format_info in product(l_in_format_infos, l_out_format_infos):

            # Filter for converters which can perform this conversion
            l_converter_names = self._get_possible_converters(in_format_info, out_format_info, only=only)

            for converter_name in l_converter_names:
                l_possible_conversions.append((self.parent.get_converter_info(converter_name),
                                               in_format_info, out_format_info))

        return l_possible_conversions

    @lru_cache
    def _get_shared_attrs(self, source_format, target_format):
        """Get a list of attributes that both the source and target format feature
        """
        source_format_info = self.parent.get_format_info(source_format)
        target_format_info = self.parent.get_format_info(target_format)

        l_shared_attrs: list[str] = []

        for attr in FormatInfo.D_PROPERTY_ATTRS:
            if getattr(source_format_info, attr) and getattr(target_format_info, attr):
                l_shared_attrs.append(attr)

        return l_shared_attrs

    def _get_info_loss(self, path):
        """Get the number of attributes in both the first and last format which would be lost if a conversion path
        is traversed
        """
        l_shared_attrs = self._get_shared_attrs(path[0], path[-1])

        if len(l_shared_attrs) == 0:
            return 0

        l_kept_attrs = copy(l_shared_attrs)
        for i in range(len(path)-1):
            target_format_info = self.parent.get_format_info(i+1)

            # Check if each attr still in the shared list is kept here
            for attr in l_kept_attrs:
                if not getattr(target_format_info, attr):
                    l_kept_attrs.remove(attr)
                    if len(l_kept_attrs) == 0:
                        break

        num_lost_attrs = len(l_shared_attrs) - len(l_kept_attrs)

        return num_lost_attrs

    def get_conversion_pathway(self,
                               in_format: str | int | FormatInfo,
                               out_format: str | int | FormatInfo,
                               only: Literal["all"] | Literal["supported"] | Literal["registered"] = "all"
                               ) -> list[tuple[ConverterInfo, FormatInfo, FormatInfo]] | None:
        """Gets a pathway to convert from one format to another
        """

        in_format_info = self.parent.get_format_info(in_format)
        out_format_info = self.parent.get_format_info(out_format)

        # Check if the formats are the same
        if in_format_info is out_format_info:
            return None

        # First check if direct conversion is possible
        l_possible_direct_conversions = self.get_possible_conversions(in_format=in_format, out_format=out_format)
        if l_possible_direct_conversions:
            # TODO: When there's some better measure of conversion quality, use it to choose which converter to use
            return [l_possible_direct_conversions[0]]

        graph: ig.Graph = self._get_desired_graph(only)

        # Query the graph for the shortest paths to perform this conversion. If no conversions are possible, igraph
        # will print a warning, which we catch and suppress here
        with catch_warnings(record=True) as l_warnings:
            l_paths: list[list[int]] = graph.get_shortest_paths(in_format_info.id, to=out_format_info.id)
            for warning in l_warnings:
                if "Couldn't reach some vertices" not in str(warning.message):
                    print(warning, file=sys.stderr)

        # Check if any paths are possible
        if not l_paths or not l_paths[0]:
            return None

        # Check each path to find the first which doesn't lose any unnecessary info, or else the one which loses the
        # least
        best_path: list[int] | None = None
        best_info_loss: int | None = None
        for path in l_paths:
            info_loss = self._get_info_loss(path)
            if best_info_loss is None or info_loss < best_info_loss:
                best_path = path
                best_info_loss = info_loss
                if best_info_loss == 0:
                    break

        # Output the best path in the desired format
        l_steps: list[tuple[str, FormatInfo, FormatInfo]] = []
        for i in range(len(best_path)-1):
            source_id: int = best_path[i]
            target_id: int = best_path[i+1]
            converter_name: str = graph.es.select(_source=source_id, _target=target_id)[0][DB_NAME_KEY]
            l_steps.append((get_converter_info(converter_name),
                            self.parent.get_format_info(source_id),
                            self.parent.get_format_info(target_id)))

        return l_steps

    def get_possible_formats(self, converter_name: str) -> tuple[list[FormatInfo], list[FormatInfo]]:
        """Get a list of input and output formats that a given converter supports

        Parameters
        ----------
        converter_name : str
            The name of the converter

        Returns
        -------
        tuple[list[FormatInfo], list[FormatInfo]]
            A tuple of a list of the supported input formats and a list of the supported output formats
        """
        conv_id: int = self.parent.get_converter_info(converter_name).id

        l_conversion_edges = self.graph.es.select(**{DB_CONV_ID_KEY: conv_id})
        l_possible_in_format_ids = list({x.source for x in l_conversion_edges})
        l_possible_out_format_ids = list({x.target for x in l_conversion_edges})

        # Get the name for each format ID, and return lists of the names
        return ([self.parent.get_format_info(x) for x in l_possible_in_format_ids],
                [self.parent.get_format_info(x) for x in l_possible_out_format_ids])


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
        self.formats: list[dict[str, bool | int | str | None]] = d_data[DB_FORMATS_KEY]
        self.converters: list[dict[str, bool | int | str | None]] = d_data[DB_CONVERTERS_KEY]
        self.converts_to: list[dict[str, bool | int | str | None]] = d_data[DB_CONVERTS_TO_KEY]

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
                name: str = regularize_name(d_single_converter_info[DB_NAME_KEY])
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
            max_id: int = max([x[DB_ID_KEY] for x in self.converters])
            self._l_converter_info: list[ConverterInfo | None] = [None] * (max_id+1)

            # Fill the list with all converters in the dict
            for single_converter_info in self.d_converter_info.values():
                self._l_converter_info[single_converter_info.id] = single_converter_info

        return self._l_converter_info

    @property
    def d_format_info(self) -> dict[str, list[FormatInfo]]:
        """Generate the format info dict when needed
        """
        if self._d_format_info is None:
            self._init_formats_and_conversions()

        return self._d_format_info

    @property
    def l_format_info(self) -> list[FormatInfo | None]:
        """Generate the format info list (indexed by ID) when needed
        """
        if self._l_format_info is None:
            self._init_formats_and_conversions()

        return self._l_format_info

    @property
    def conversions_table(self) -> ConversionsTable:
        """Generates the conversions table when needed
        """

        if self._conversions_table is None:
            self._init_formats_and_conversions()

        return self._conversions_table

    def _init_formats_and_conversions(self):
        """Initializes the format list and dict and the conversions table"""

        # Start by initializing the list of conversions

        # Pre-size a list based on the maximum ID plus 1 (since IDs are 1-indexed)
        max_id: int = max([x[DB_ID_KEY] for x in self.formats])
        self._l_format_info: list[FormatInfo | None] = [None] * (max_id+1)

        for d_single_format_info in self.formats:
            lc_name: str = d_single_format_info[DB_FORMAT_EXT_KEY]

            format_info = FormatInfo(name=lc_name,
                                     parent=self,
                                     d_single_format_info=d_single_format_info)

            self._l_format_info[format_info.id] = format_info

        # Initialize the conversions table now
        self._conversions_table = ConversionsTable(l_converts_to=self.converts_to,
                                                   parent=self)

        # Use the conversions graph to prune any formats which have no valid conversions

        # Get a slice of the table which only includes supported converters
        supported_graph = self._conversions_table.supported_graph

        for format_id, format_info in enumerate(self._l_format_info):
            if not format_info:
                continue

            # Check if the format is supported as the input or output format for any conversion
            if supported_graph.degree(format_id) > 0:
                continue

            # If we get here, the format isn't supported for any conversions, so remove it from our list
            self._l_format_info[format_id] = None

        # Now create the formats dict, with only the pruned list of formats
        self._d_format_info: dict[str, list[FormatInfo]] = {}

        for format_info in self.l_format_info:

            if not format_info:
                continue

            lc_name = format_info.name.lower()

            # Each name may correspond to multiple formats, so we use a list for each entry to list all possible
            # formats for each name
            if lc_name not in self._d_format_info:
                self._d_format_info[lc_name] = []

            self._d_format_info[lc_name].append(format_info)

    def get_converter_info(self, converter_name_or_id: str | int) -> ConverterInfo:
        """Get a converter's info from either its name or ID
        """
        if isinstance(converter_name_or_id, str):
            try:
                return self.d_converter_info[converter_name_or_id]
            except KeyError:
                raise FileConverterDatabaseException(f"Converter name '{converter_name_or_id}' not recognised",
                                                     help=True)
        elif isinstance(converter_name_or_id, int):
            return self.l_converter_info[converter_name_or_id]
        else:
            raise FileConverterDatabaseException(f"Invalid key passed to `get_converter_info`: '{converter_name_or_id}'"
                                                 f" of type '{type(converter_name_or_id)}'. Type must be `str` or "
                                                 "`int`")

    @overload
    def get_format_info(self,
                        format_name_or_id: str | int | FormatInfo,
                        which: int | None = None) -> FormatInfo: ...

    @overload
    def get_format_info(self,
                        format_name_or_id: str | int | FormatInfo,
                        which: Literal["all"]) -> list[FormatInfo]: ...

    def get_format_info(self,
                        format_name_or_id: str | int | FormatInfo,
                        which: int | Literal["all"] | None = None) -> FormatInfo | list[FormatInfo]:
        """Gets the information on a given file format stored in the database

        Parameters
        ----------
        format_name_or_id : str | int | FormatInfo
            The name (extension) of the format, or its ID. In the case of ambiguous extensions which could apply to
            multiple formats, the ID must be used here or a FileConverterDatabaseException will be raised. This also
            allows passing a FormatInfo to this, in which case that object will be silently returned, to allow
            normalising the input to always be a FormatInfo when output from this
        which : int | None
            In the case that an extension string is provided which turns out to be ambiguous, which of the listed
            possibilities to use from the zero-indexed list. Default None, which raises an exception for an ambiguous
            format. 0 may be used to select the first in the database, which is often a good default choice. The literal
            string "all" may be used to request all possibilites, in which case this method will return a list (even if
            there are zero or one possibilities)

        Returns
        -------
        FormatInfo | list[FormatInfo]
        """

        if which == "all":
            return_as_list = True
        else:
            return_as_list = False

        if isinstance(format_name_or_id, str):
            # Silently strip leading period
            if format_name_or_id.startswith("."):
                format_name_or_id = format_name_or_id[1:]

            # Convert the format name to lower-case to handle it case-insensitively
            format_name_or_id = format_name_or_id.lower()

            # Check for a hyphen in the format, which indicates a preference from the user as to which, overriding the
            # `which` kwarg
            if "-" in format_name_or_id:
                l_name_segments = format_name_or_id.split("-")
                if len(l_name_segments) > 2:
                    raise FileConverterDatabaseException(f"Format name '{format_name_or_id} is improperly formatted - "
                                                         "It may contain at most one hyphen, separating the extension "
                                                         "from an index indicating which of the formats with that "
                                                         "extension to use, e.g. 'pdb-0', 'pdb-1', etc.",
                                                         help=True)
                format_name_or_id = l_name_segments[0]
                which = int(l_name_segments[1])

            l_possible_format_info = self.d_format_info.get(format_name_or_id, [])

            if which == "all":
                return l_possible_format_info

            elif len(l_possible_format_info) == 1:
                format_info = l_possible_format_info[0]

            elif len(l_possible_format_info) == 0:
                raise FileConverterDatabaseException(f"Format name '{format_name_or_id}' not recognised",
                                                     help=True)

            elif which is not None and which < len(l_possible_format_info):
                format_info = l_possible_format_info[which]

            else:
                msg = (f"Extension '{format_name_or_id}' is ambiguous and must be defined by ID. Possible formats "
                       "and their IDs are:")
                for possible_format_info in l_possible_format_info:
                    msg += (f"\n{possible_format_info.id}: {possible_format_info.disambiguated_name} "
                            f"({possible_format_info.note})")
                raise FileConverterDatabaseException(msg, help=True)

        elif isinstance(format_name_or_id, int):
            try:
                format_info = self.l_format_info[format_name_or_id]
            except IndexError:
                if return_as_list:
                    return []
                raise FileConverterDatabaseException(f"Format ID '{format_name_or_id}' not recognised",
                                                     help=True)

        elif isinstance(format_name_or_id, FormatInfo):
            # Silently return the FormatInfo if it was used as a key here
            format_info = format_name_or_id

        else:
            raise FileConverterDatabaseException(f"Invalid key passed to `get_format_info`: '{format_name_or_id}'"
                                                 f" of type '{type(format_name_or_id)}'. Type must be `str` or "
                                                 "`int`")
        if return_as_list:
            return [format_info]

        return format_info


# The database will be loaded on demand when `get_database()` is called
_database: DataConversionDatabase | None = None


def get_database_path() -> str:
    """Get the absolute path to the database file

    Returns
    -------
    str
    """

    # For an interactive shell, __file__ won't be defined for this module, so use the constants module instead
    reference_file = os.path.realpath(const.__file__)

    qualified_database_filename = os.path.join(os.path.dirname(reference_file), const.DATABASE_FILENAME)

    return qualified_database_filename


def load_database() -> DataConversionDatabase:
    """Load and return a new instance of the data conversion database from the JSON database file in this package. This
    function should not be called directly unless you specifically need a new instance of the database object and can't
    deepcopy the database returned by `get_database()`, as it's expensive to load it in.

    Returns
    -------
    DataConversionDatabase
    """

    # Find and load the database JSON file
    d_data: dict = json.load(open(get_database_path(), "r"))

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

    return get_database().d_converter_info[regularize_name(name)]


@overload
def get_format_info(format_name_or_id: str | int | FormatInfo,
                    which: int | None = None) -> FormatInfo: ...


@overload
def get_format_info(format_name_or_id: str | int | FormatInfo,
                    which: Literal["all"]) -> list[FormatInfo]: ...


def get_format_info(format_name_or_id: str | int | FormatInfo,
                    which: int | Literal["all"] | None = None) -> FormatInfo | list[FormatInfo]:
    """Gets the information on a given file format stored in the database

    Parameters
    ----------
    format_name_or_id : str | int | FormatInfo
        The name (extension) of the format, or its ID. In the case of ambiguous extensions which could apply to multiple
        formats, the ID must be used here or a FileConverterDatabaseException will be raised. This also allows passing a
        FormatInfo to this, in which case that object will be silently returned, to allow normalising the input to
        always be a FormatInfo when output from this
    which : int | None
        In the case that an extension string is provided which turns out to be ambiguous, which of the listed
        possibilities to use from the zero-indexed list. Default None, which raises an exception for an ambiguous
        format. 0 may be used to select the first in the database, which is often a good default choice. The literal
        string "all" may be used to request all possibilites, in which case this method will return a list (even if
        there are zero or one possibilities)

    Returns
    -------
    FormatInfo | list[FormatInfo]
    """

    return get_database().get_format_info(format_name_or_id, which)


def get_conversion_quality(converter_name: str,
                           in_format: str | int,
                           out_format: str | int) -> ConversionQualityInfo | None:
    """Get an indication of the quality of a conversion from one format to another, or if it's not possible

    Parameters
    ----------
    converter_name : str
        The name of the converter
    in_format : str | int
        The extension or ID of the input file format
    out_format : str | int
        The extension or ID of the output file format

    Returns
    -------
    ConversionQualityInfo | None
        If the conversion is not possible, returns None. If the conversion is possible, returns a
        `ConversionQualityInfo` object with info on the conversion
    """

    return get_database().conversions_table.get_conversion_quality(converter_name=regularize_name(converter_name),
                                                                   in_format=in_format,
                                                                   out_format=out_format)


def get_possible_conversions(in_format: str | int,
                             out_format: str | int) -> list[tuple[ConverterInfo, FormatInfo, FormatInfo]]:
    """Get a list of converters which can perform a conversion from one format to another and disambiguate in the case
    of ambiguous input/output formats

    Parameters
    ----------
    in_format : str | int
        The extension or ID of the input file format
    out_format : str | int
        The extension or ID of the output file format

    Returns
    -------
    list[tuple[ConverterInfo, FormatInfo, FormatInfo]]
        A list of tuples, where each tuple's first item is the ConverterInfo of a converter which can perform a matching
        conversion, the second is the info of the input format for this conversion, and the third is the info of the
        output format
    """

    return get_database().conversions_table.get_possible_conversions(in_format=in_format,
                                                                     out_format=out_format)


def get_conversion_pathway(in_format: str | int | FormatInfo,
                           out_format: str | int | FormatInfo,
                           only: Literal["all"] | Literal["supported"] | Literal["registered"] = "all"
                           ) -> list[tuple[ConverterInfo, FormatInfo, FormatInfo]] | None:
    """Get a list of conversions that can be performed to convert one format to another. This is primarily used when a
    direct conversion is not supported by any individual converter. Only one possible pathway will be returned,
    prioritising pathways which do not lose lose and then re-extrapolate any information stored by some formats and not
    others along the path.

    Parameters
    ----------
    in_format : str | int
        The input file format. For this function, the format must be defined uniquely, either by using a disambiguated
        extension, ID, or FormatInfo
    out_format : str | int
        The output file format. For this function, the format must be defined uniquely, either by using a disambiguated
        extension, ID, or FormatInfo
    only : Literal["all"] | Literal["supported"] | Literal["registered"], optional
        Which converters to limit the pathway search to:
        - "all" (default): All known converters
        - "supported": Only converters supported by this utility, even if not currently available (e.g. they don't work
          on your OS)
        - "registered": Only converters supported by this utility and currently available

    Returns
    -------
    list[tuple[ConverterInfo, FormatInfo, FormatInfo]] | None
        Will return `None` if no conversion pathway is possible or if the input and output formats are the same.
        Otherwise, will return a list of steps in the pathway, each being a tuple of:

        converter_info : ConverterInfo
            Info on the converter used to perform this step
        in_format : FormatInfo
            Input format for this step (if the first step, will be the input format to this function, otherwise will be
            the output format of the previous step)
        out_format : FormatInfo
            Output format from this step (if the last step, will be the output format for this function, otherwise will
            be the input format of the next step)
    """

    return get_database().conversions_table.get_conversion_pathway(in_format=in_format,
                                                                   out_format=out_format,
                                                                   only=only)


def disambiguate_formats(converter_name: str,
                         in_format: str | int | FormatInfo,
                         out_format: str | int | FormatInfo) -> tuple[FormatInfo, FormatInfo]:
    """Try to disambiguate formats by seeing if there's only one possible conversion between formats matching those
    provided.

    Parameters
    ----------
    converter_name : str
        The name of the converter
    in_format : str | int
        The extension or ID of the input file format
    out_format : str | int
        The extension or ID of the output file format

    Returns
    -------
    tuple[FormatInfo, FormatInfo]
        The input and output format for this conversion, if only one combination is possible

    Raises
    ------
    FileConverterDatabaseException
        If more than one format combination is possible for this conversion, or no conversion is possible
    """

    # Regularize the converter name so we don't worry about case/spacing mismatches
    converter_reg_name = regularize_name(converter_name)

    # Get all possible conversions, and see if we only have one for this converter
    l_possible_conversions = [x for x in get_possible_conversions(in_format, out_format)
                              if x[0].name == converter_reg_name]

    if len(l_possible_conversions) == 1:
        return l_possible_conversions[0][1], l_possible_conversions[0][2]
    elif len(l_possible_conversions) == 0:
        raise FileConverterDatabaseException(f"Conversion from {in_format} to {out_format} with converter "
                                             f"{converter_name} is not supported", help=True)
    else:
        msg = (f"Conversion from {in_format} to {out_format} with converter {converter_name} is ambiguous.\n"
               "Possible matching conversions are:\n")
        for _, possible_in_format, possible_out_format in l_possible_conversions:
            msg += (f"{possible_in_format.disambiguated_name} ({possible_in_format.note}) to "
                    f"{possible_out_format.disambiguated_name} ({possible_out_format.note})\n")
        # Trim the final newline from the message
        msg = msg[:-1]
        raise FileConverterDatabaseException(msg, help=True)


def get_possible_formats(converter_name: str) -> tuple[list[FormatInfo], list[FormatInfo]]:
    """Get a list of input and output formats that a given converter supports

    Parameters
    ----------
    converter_name : str
        The name of the converter

    Returns
    -------
    tuple[list[FormatInfo], list[FormatInfo]]
        A tuple of a list of the supported input formats and a list of the supported output formats
    """
    return get_database().conversions_table.get_possible_formats(converter_name=regularize_name(converter_name))


def _find_arg(tl_args: tuple[list[FlagInfo], list[OptionInfo]],
              arg: str) -> ArgInfo:
    """Find a specific flag or option in the lists
    """
    for l_args in tl_args:
        l_found = [x for x in l_args if x.flag == arg]
        if len(l_found) > 0:
            return l_found[0]
    # If we get here, it wasn't found in either list
    raise FileConverterDatabaseException(f"Argument '{arg}' was not found in the list of allowed arguments for this "
                                         "conversion")


def get_in_format_args(converter_name: str,
                       format_name: str,
                       arg: str | None = None) -> tuple[list[FlagInfo], list[OptionInfo]] | ArgInfo:
    """Get the input flags and options supported by a given converter for a given format (provided as its extension).
    Optionally will provide information on just a single flag or option if its value is provided as an optional argument

    Parameters
    ----------
    converter_name : str
        The converter name
    format_name : str
        The file format name (extension)
    arg : str | None
        If provided, only information on this flag or option will be provided

    Returns
    -------
    tuple[set[FlagInfo], set[OptionInfo]]
        A list of info for the allowed flags, and a set of info for the allowed options
    """

    converter_info = get_converter_info(converter_name)
    tl_args = converter_info.get_in_format_args(format_name)
    if not arg:
        return tl_args
    return _find_arg(tl_args, arg)


def get_out_format_args(converter_name: str,
                        format_name: str,
                        arg: str | None = None) -> tuple[list[FlagInfo], list[OptionInfo]]:
    """Get the output flags and options supported by a given converter for a given format (provided as its extension).
    Optionally will provide information on just a single flag or option if its value is provided as an optional argument

    Parameters
    ----------
    converter_name : str
        The converter name
    format_name : str
        The file format name (extension)
    arg : str | None
        If provided, only information on this flag or option will be provided

    Returns
    -------
    tuple[set[FlagInfo], set[OptionInfo]]
        A list of info for the allowed flags, and a set of info for the allowed options
    """

    converter_info = get_converter_info(converter_name)
    tl_args = converter_info.get_out_format_args(format_name)
    if not arg:
        return tl_args
    return _find_arg(tl_args, arg)
