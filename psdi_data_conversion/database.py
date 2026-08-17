"""@file psdi_data_conversion/database.py

Created 2025-02-03 by Bryan Gillis.

Python module provide utilities for accessing the converter database
"""

from __future__ import annotations

import json
import sys
import warnings
from copy import copy
from dataclasses import dataclass, field
from functools import lru_cache
from itertools import product
from logging import getLogger
from pathlib import Path
from typing import Any, Literal, NamedTuple, overload
from uuid import UUID
from warnings import catch_warnings

import igraph as ig

from psdi_data_conversion import constants as const
from psdi_data_conversion.converter import (L_REGISTERED_CONVERTERS, L_SUPPORTED_CONVERTERS,
                                            get_registered_converter_class)
from psdi_data_conversion.converters.base import FileConverter, FileConverterException
from psdi_data_conversion.file_io import get_package_path
from psdi_data_conversion.utils import regularize_name

# Database keys
# -------------

# Keys for converter-specific databases
DB_CONVERTER_KEY = "converter"
DB_KEY_PREFIX_KEY = "database_key_prefix"
DB_DESC_KEY = "desc"
DB_INFO_KEY = "info"
DB_SUPPORT_AMBIG_EXT_KEY = "supports_ambiguous_extensions"

DB_EXTRA_FORMATS_KEY = "extra_formats"

DB_SUPPORTED_FORMATS_KEY = "supported_formats"
DB_IN_ONLY_FORMATS_KEY = "in_only_formats"
DB_OUT_ONLY_FORMATS_KEY = "out_only_formats"

DB_SUPPORTED_CONVERSIONS_KEY = "supported_conversions"
DB_UNSUPPORTED_CONVERSIONS_KEY = "unsupported_conversions"

DB_IN_FLAGS_KEY = "in_flags"
DB_OUT_FLAGS_KEY = "out_flags"
DB_IN_OPTIONS_KEY = "in_options"
DB_OUT_OPTIONS_KEY = "out_options"

DB_FORMAT_ID_LIST_KEY = "format_ids"

# Keys for top-level and general items in the database
DB_FORMATS_KEY = "formats"
DB_CONVERTERS_KEY = "converters"
DB_CONVERTS_TO_KEY = "converts_to"
DB_ID_KEY = "id"
DB_NAME_KEY = "name"

# Keys for converter general info in the database
DB_DESCRIPTION_KEY = "description"
DB_FURTHER_INFO_KEY = "further_info"
DB_URL_KEY = "url"

# Keys for format general info in the database - some are duplicated here so they're also stored in the same format as
# other keys here
DB_FORMAT_EXT_KEY = "extension"
DB_FORMAT_C2X_KEY = "format"
DB_FORMAT_NOTE_KEY = "note"
DB_FORMAT_COMP_KEY = const.QUAL_COMP_KEY
DB_FORMAT_CONN_KEY = const.QUAL_CONN_KEY
DB_FORMAT_2D_KEY = const.QUAL_2D_KEY
DB_FORMAT_3D_KEY = const.QUAL_3D_KEY
DB_FORMAT_PRECISION_KEY = "precision"
DB_FORMAT_CONFIRMED_NEW_KEY = "confirmed_new"

# Keys for converts_to info in the database
DB_CONV_ID_KEY = "converters_id"
DB_IN_ID_KEY = "in_id"
DB_OUT_ID_KEY = "out_id"
DB_SUCCESS_KEY = "degree_of_success"
DB_WEIGHT_KEY = "weight"

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

# Chaining constants
# ------------------

# Each format property is assigned a weight with a different power of 2, plus a weight for taking any conversion step at
# all, to account for miscellaneous lossiness from a conversion that can't be quantified
D_PROP_BITS = {
    const.QUAL_COMP_KEY: 24,
    const.QUAL_CONN_KEY: 18,
    const.QUAL_2D_KEY: 12,
    const.QUAL_3D_KEY: 6
}
D_PROP_WEIGHTS = {key: 1 << bit for key, bit in D_PROP_BITS.items()}

STEP_BIT = 0
STEP_WEIGHT = 1 << STEP_BIT

# Number of bits the property weight section is offset within the full weight when everything is combined into a single
# 128-bit integer
PROP_WEIGHT_BIT_OFFSET = 64

# Minimum and maximum for digits of precision lost
PREC_MIN_DIGIT_LOSS = 0
PREC_MAX_DIGIT_LOSS = 12

# Number of bits separating weight bits for different levels of precision loss
PREC_GAP_BITS = 3

# Number of bits the precision weight section is offset within the full weight when everything is combined into a single
# 128-bit integer
PREC_WEIGHT_BIT_OFFSET = 16

# Number of bits the time weight section is offset within the full weight when everything is combined into a single
# 128-bit integer
TIME_WEIGHT_BIT_OFFSET = 0

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

    @property
    def uuid(self):
        """Returns the ID as a UUID object
        """
        return UUID(int=self.id)


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
            self.pretty_name = self.converter_class.meta.name
        except KeyError:
            self.converter_class = None
            self.pretty_name = name

        self.parent = parent
        """The parent database"""

        # Get info about the converter from the database
        self.id: int = d_single_converter_info.get(DB_ID_KEY, -1)
        """The converter's ID"""

        self.description: str = d_single_converter_info.get(DB_DESCRIPTION_KEY, "")
        """A description of the converter"""

        self.url: str = d_single_converter_info.get(DB_URL_KEY, "")
        """The official URL for the converter"""

        # Get necessary info about the converter from the class

        self.supported = self.name in L_SUPPORTED_CONVERTERS
        """Whether or not the converter is supported by this package. If a converter is supported but not registered,
        this usually means that a required binary is missing and must be supplied by the user"""

        self.registered = self.name in L_REGISTERED_CONVERTERS
        """Whether or not the converter is ready to be used by this package"""

        try:
            self._key_prefix = get_registered_converter_class(name).meta.database_key_prefix
        except KeyError:
            # We'll get a KeyError for converters in the database that don't yet have their own class, which we can
            # safely ignore
            self._key_prefix = None

        self._arg_info: dict[str, list[dict[str, int | str]]] = {}

        # Placeholders for members that are generated when needed
        self._d_in_flag_info: dict[int, FlagInfo] | None = None
        self._l_unsorted_in_flag_info: list[FlagInfo] | None = None
        self._d_out_flag_info: dict[int, FlagInfo] | None = None
        self._l_unsorted_out_flag_info: list[FlagInfo] | None = None
        self._d_in_option_info: dict[int, OptionInfo] | None = None
        self._l_unsorted_in_option_info: list[FlagInfo] | None = None
        self._d_out_option_info: dict[int, OptionInfo] | None = None
        self._l_unsorted_out_option_info: list[FlagInfo] | None = None

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

    @property
    def uuid(self):
        """Returns the ID as a UUID object
        """
        return UUID(int=self.id)

    def _create_d_arg_info(self, subclass: type[ArgInfo]):
        """Creates either the flag or option info dicts when needed
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
            raise FileConverterDatabaseException(f"Unrecognised subclass passed to `_create_d_arg_info`: {subclass}")

        for key_base, in_or_out in ((in_key_base, "in"),
                                    (out_key_base, "out")):

            d_arg_info: dict[int, ArgInfo] = {}

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
                                    description=d_single_arg_info[DB_DESCRIPTION_KEY],
                                    info=d_single_arg_info[DB_FURTHER_INFO_KEY],
                                    **optional_arg_info_kwargs)
                d_arg_info[arg_id] = arg_info

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

            if in_or_out == "in" and issubclass(subclass, FlagInfo):
                self._d_in_flag_info = d_arg_info
                self._l_unsorted_in_flag_info = list(d_arg_info.values())
            elif in_or_out == "out" and issubclass(subclass, FlagInfo):
                self._d_out_flag_info = d_arg_info
                self._l_unsorted_out_flag_info = list(d_arg_info.values())
            elif in_or_out == "in" and issubclass(subclass, OptionInfo):
                self._d_in_option_info = d_arg_info
                self._l_unsorted_in_option_info = list(d_arg_info.values())
            elif in_or_out == "out" and issubclass(subclass, OptionInfo):
                self._d_out_option_info = d_arg_info
                self._l_unsorted_out_option_info = list(d_arg_info.values())
            else:
                raise FileConverterDatabaseException(
                    f"Unrecognised subclass passed to `_create_d_arg_info`: {subclass}")

        return

    @property
    def d_in_flag_info(self) -> dict[int, FlagInfo] | None:
        """Generate the input flag info dict (indexed by ID) when needed. Returns None if the converter has no flag info
        in the database
        """
        if self._d_in_flag_info is None and self._key_prefix is not None:
            self._create_d_arg_info(FlagInfo)
        return self._d_in_flag_info

    @property
    def l_in_flag_info(self):
        """DEPRECATED: Generate the input flag info list (indexed by ID) when needed. Returns None if the converter has
        no flag info in the database
        """
        deprecation_msg = ("`l_in_flag_info` is deprecated as of version 0.4.0 and due to be removed in a future "
                           "release. To get an input FlagInfo from the format UUID, use `d_in_flag_info`. To "
                           "get an unsorted list of input FlagInfo, use `l_unsorted_in_flag_info`.")
        warnings.warn(deprecation_msg, DeprecationWarning)
        raise AttributeError(deprecation_msg)

    @property
    def l_unsorted_in_flag_info(self) -> list[FlagInfo] | None:
        """Generate the unsorted input flag info list when needed. Returns None if the converter has
        no flag info in the database
        """
        if self._l_unsorted_in_flag_info is None and self._key_prefix is not None:
            self._create_d_arg_info(FlagInfo)
        return self._l_unsorted_in_flag_info

    @property
    def d_out_flag_info(self) -> dict[int, FlagInfo] | None:
        """Generate the input flag info dict (indexed by ID) when needed. Returns None if the converter has no flag info
        in the database
        """
        if self._d_out_flag_info is None and self._key_prefix is not None:
            self._create_d_arg_info(FlagInfo)
        return self._d_out_flag_info

    @property
    def l_out_flag_info(self):
        """DEPRECATED: Generate the input flag info list (indexed by ID) when needed. Returns None if the converter has
        no flag info in the database
        """
        deprecation_msg = ("`l_out_flag_info` is deprecated as of version 0.4.0 and due to be removed in a future "
                           "release. To get an input FlagInfo from the format UUID, use `d_out_flag_info`. To "
                           "get an unsorted list of input FlagInfo, use `l_unsorted_out_flag_info`.")
        warnings.warn(deprecation_msg, DeprecationWarning)
        raise AttributeError(deprecation_msg)

    @property
    def l_unsorted_out_flag_info(self) -> list[FlagInfo] | None:
        """Generate the unsorted input flag info list when needed. Returns None if the converter has
        no flag info in the database
        """
        if self._l_unsorted_out_flag_info is None and self._key_prefix is not None:
            self._create_d_arg_info(FlagInfo)
        return self._l_unsorted_out_flag_info

    @property
    def d_in_option_info(self) -> dict[int, OptionInfo] | None:
        """Generate the input option info dict (indexed by ID) when needed. Returns None if the converter has no option
        info in the database
        """
        if self._d_in_option_info is None and self._key_prefix is not None:
            self._create_d_arg_info(OptionInfo)
        return self._d_in_option_info

    @property
    def l_in_option_info(self):
        """DEPRECATED: Generate the input option info list (indexed by ID) when needed. Returns None if the converter
        has no option info in the database
        """
        deprecation_msg = ("`l_in_option_info` is deprecated as of version 0.4.0 and due to be removed in a future "
                           "release. To get an input OptionInfo from the format UUID, use `d_in_option_info`. To "
                           "get an unsorted list of input OptionInfo, use `l_unsorted_in_option_info`.")
        warnings.warn(deprecation_msg, DeprecationWarning)
        raise AttributeError(deprecation_msg)

    @property
    def l_unsorted_in_option_info(self) -> list[OptionInfo] | None:
        """Generate the unsorted input option info list when needed. Returns None if the converter has
        no option info in the database
        """
        if self._l_unsorted_in_option_info is None and self._key_prefix is not None:
            self._create_d_arg_info(OptionInfo)
        return self._l_unsorted_in_option_info

    @property
    def d_out_option_info(self) -> dict[int, OptionInfo] | None:
        """Generate the input option info dict (indexed by ID) when needed. Returns None if the converter has no option
        info in the database
        """
        if self._d_out_option_info is None and self._key_prefix is not None:
            self._create_d_arg_info(OptionInfo)
        return self._d_out_option_info

    @property
    def l_out_option_info(self):
        """DEPRECATED: Generate the input option info list (indexed by ID) when needed. Returns None if the converter
        has no option info in the database
        """
        deprecation_msg = ("`l_out_option_info` is deprecated as of version 0.4.0 and due to be removed in a future "
                           "release. To get an input OptionInfo from the format UUID, use `d_out_option_info`. To "
                           "get an unsorted list of input OptionInfo, use `l_unsorted_out_option_info`.")
        warnings.warn(deprecation_msg, DeprecationWarning)
        raise AttributeError(deprecation_msg)

    @property
    def l_unsorted_out_option_info(self) -> list[OptionInfo] | None:
        """Generate the unsorted input option info list when needed. Returns None if the converter has
        no option info in the database
        """
        if self._l_unsorted_out_option_info is None and self._key_prefix is not None:
            self._create_d_arg_info(OptionInfo)
        return self._l_unsorted_out_option_info

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
            l_arg_info = self.l_unsorted_in_flag_info if in_or_out == "in" else self.l_unsorted_out_flag_info
        elif issubclass(subclass, OptionInfo):
            l_arg_info = self.l_unsorted_in_option_info if in_or_out == "in" else self.l_unsorted_out_option_info
        else:
            raise FileConverterDatabaseException(
                f"Unrecognised subclass passed to `_create_d_format_args`: {subclass}")

        d_format_args: dict[str | int, set[ArgInfo]] = {}
        d_parent_format_info_from_id = self.parent.d_format_info_from_id

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
            l_format_info = [d_parent_format_info_from_id[format_id] for format_id in s_formats]
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

    def get_in_format_args(self, in_format: str | int | UUID | FormatInfo) -> tuple[list[FlagInfo], list[OptionInfo]]:
        """Get the input flags and options supported for a given format (provided as its extension)

        Parameters
        ----------
        in_format : str
            The file format name (extension), ID, or FormatInfo

        Returns
        -------
        tuple[set[FlagInfo], set[OptionInfo]]
            A set of info for the allowed flags, and a set of info for the allowed options
        """

        l_in_format_infos = get_format_info(in_format, which="all")
        s_flag_ids = set()
        s_option_ids = set()

        for in_format_info in l_in_format_infos:
            in_format_id = in_format_info.id

            s_flag_ids.update(self.d_in_format_flags.get(in_format_id, set()))
            s_option_ids.update(self.d_in_format_options.get(in_format_id, set()))

        l_flag_ids = list(s_flag_ids)
        l_flag_ids.sort()
        l_flag_info = [self.d_in_flag_info[x] for x in l_flag_ids]

        l_option_ids = list(s_option_ids)
        l_option_ids.sort()
        l_option_info = [self.d_in_option_info[x] for x in l_option_ids]

        return l_flag_info, l_option_info

    def get_out_format_args(self, out_format: str | int | UUID | FormatInfo) -> tuple[list[FlagInfo], list[OptionInfo]]:
        """Get the output flags and options supported for a given format (provided as its extension)

        Parameters
        ----------
        out_format : str
            The file format name (extension), ID, or FormatInfo

        Returns
        -------
        tuple[set[FlagInfo], set[OptionInfo]]
            A set of info for the allowed flags, and a set of info for the allowed options
        """

        l_out_format_infos = get_format_info(out_format, which="all")
        s_flag_ids = set()
        s_option_ids = set()

        for out_format_info in l_out_format_infos:
            out_format_id = out_format_info.id

            s_flag_ids.update(self.d_out_format_flags.get(out_format_id, set()))
            s_option_ids.update(self.d_out_format_options.get(out_format_id, set()))

        l_flag_ids = list(s_flag_ids)
        l_flag_ids.sort()
        l_flag_info = [self.d_out_flag_info[x] for x in l_flag_ids]

        l_option_ids = list(s_option_ids)
        l_option_ids.sort()
        l_option_info = [self.d_out_option_info[x] for x in l_option_ids]

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

        self.composition: bool | None = d_single_format_info.get(DB_FORMAT_COMP_KEY)
        """Whether or not this format stores composition information"""

        self.connections: bool | None = d_single_format_info.get(DB_FORMAT_CONN_KEY)
        """Whether or not this format stores connections information"""

        self.two_dim: bool | None = d_single_format_info.get(DB_FORMAT_2D_KEY)
        """Whether or not this format stores 2D structural information"""

        self.three_dim: bool | None = d_single_format_info.get(DB_FORMAT_3D_KEY)
        """Whether or not this format stores 3D structural information"""

        self.precision: int | None = d_single_format_info.get(DB_FORMAT_PRECISION_KEY)
        """The precision of numeric information in the format, as the number of decimal places, or 0 if unknown"""

        self._lower_name: str = self.name.lower()
        """The format name all in lower-case"""

        self._disambiguated_name: str | None = None

    @property
    def uuid(self):
        """Returns the ID as a UUID object
        """
        return UUID(int=self.id)

    @property
    def disambiguated_name(self) -> str:
        """A unique name for this format which can be used to distinguish it from others which share the same extension,
        by appending the name of each with a unique index"""
        if self._disambiguated_name is None:
            l_formats_with_same_name = [x for x in self.parent.l_unsorted_format_info
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

    in_format: FormatInfo
    """The info for the input file format"""

    out_format: FormatInfo
    """The info for the output file format"""

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

        # Build the conversion graphs - each format is a vertex, each conversion is an edge

        num_formats = len(parent.formats)

        # igraph doesn't support int128s (used for UUIDs) for indices, so associate each format with an index
        self.d_indices_from_uuids: dict[int, int] = {}
        self.d_uuids_from_indices: dict[int, int] = {}

        for i, format in enumerate(parent.formats):
            uuid = format[DB_ID_KEY]
            self.d_indices_from_uuids[uuid] = i
            self.d_uuids_from_indices[i] = uuid

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
                                                         for x in parent.l_unsorted_format_info]},
                             edges=[(self.d_indices_from_uuids[x[DB_IN_ID_KEY]],
                                     self.d_indices_from_uuids[x[DB_OUT_ID_KEY]]) for x in l_conversions],
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
        l_edges = graph.es.select(_source=self.d_indices_from_uuids[in_format_info.id],
                                  _target=self.d_indices_from_uuids[out_format_info.id])
        return [x[DB_NAME_KEY] for x in l_edges]

    @lru_cache(maxsize=None)
    def get_conversion_quality(self,
                               converter_name: str,
                               in_format: str | int | UUID | FormatInfo,
                               out_format: str | int | UUID | FormatInfo) -> ConversionQualityInfo | None:
        """Get an indication of the quality of a conversion from one format to another, or if it's not possible

        Parameters
        ----------
        converter_name : str
            The name of the converter
        in_format : str | int | UUID | FormatInfo
            The extension or ID of the input file format
        out_format : str | int | UUID | FormatInfo
            The extension or ID of the output file format

        Returns
        -------
        ConversionQualityInfo | None
            If the conversion is not possible, returns None. If the conversion is possible, returns a
            `ConversionQualityInfo` object with info on the conversion
        """

        # Get all possible format infos for each format
        l_in_format_info = self.parent.get_format_info(in_format, "all")
        l_out_format_info = self.parent.get_format_info(out_format, "all")

        # First check if the conversion is possible for at least one combination
        l_found_combinations: list[tuple[FormatInfo, FormatInfo]] = []
        for in_format_info, out_format_info in product(l_in_format_info, l_out_format_info):
            if converter_name in self._get_possible_converters(in_format_info, out_format_info):
                l_found_combinations.append((in_format_info, out_format_info))
        if len(l_found_combinations) == 0:
            return None

        # Check if the conversion is ambiguous
        if len(l_found_combinations) > 1:
            msg = (f"Conversion from {in_format} to {out_format} with converter {converter_name} is ambiguous. Please "
                   "Use the ID or disambiguated name (listed below) of the desired conversion. Possible matching "
                   "conversions are:\n")
            for possible_in_format, possible_out_format in l_found_combinations:
                msg += (f"    {possible_in_format.id}: {possible_in_format.disambiguated_name} "
                        f"({possible_in_format.note}) to "
                        f"{possible_out_format.id}: {possible_out_format.disambiguated_name} "
                        f"({possible_out_format.note})\n")
            # Trim the final newline from the message
            msg = msg[:-1]
            raise FileConverterDatabaseException(msg, help=True)

        in_format_info, out_format_info = l_found_combinations[0]

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
                                     in_format=in_format_info,
                                     out_format=out_format_info,
                                     qual_str=qual_str,
                                     details=details,
                                     d_prop_conversion_info=d_prop_conversion_info)

    def get_possible_conversions(self,
                                 in_format: str | int | UUID | FormatInfo,
                                 out_format: str | int | UUID | FormatInfo,
                                 only: Literal["all"] | Literal["supported"] | Literal["registered"] = "all"
                                 ) -> list[tuple[ConverterInfo, FormatInfo, FormatInfo]]:
        """Get a list of converters which can perform a conversion from one format to another, disambiguating in the
        case of ambiguous formats and providing IDs for input/output formats for possible conversions

        Parameters
        ----------
        in_format : str | int | UUID | FormatInfo
            The extension or ID of the input file format
        out_format : str | int | UUID | FormatInfo
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
    def _get_shared_attrs(self, source_format_index, target_format_index):
        """Get a list of attributes that both the source and target format feature
        """
        source_format_info = self.parent.get_format_info(self.d_uuids_from_indices[source_format_index])
        target_format_info = self.parent.get_format_info(self.d_uuids_from_indices[target_format_index])

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
            target_format_info = self.parent.get_format_info(self.d_uuids_from_indices[path[i+1]])

            # Check if each attr still in the shared list is kept here
            for attr in l_kept_attrs:
                if not getattr(target_format_info, attr):
                    l_kept_attrs.remove(attr)
                    if len(l_kept_attrs) == 0:
                        break

        num_lost_attrs = len(l_shared_attrs) - len(l_kept_attrs)

        return num_lost_attrs

    def get_conversion_pathway(self,
                               in_format: str | int | UUID | FormatInfo,
                               out_format: str | int | UUID | FormatInfo,
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
            l_paths: list[list[int]] = graph.get_shortest_paths(self.d_indices_from_uuids[in_format_info.id],
                                                                to=self.d_indices_from_uuids[out_format_info.id])
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
            source_index = best_path[i]
            source_id: int = self.d_uuids_from_indices[source_index]
            target_index: int = best_path[i+1]
            target_id: int = self.d_uuids_from_indices[target_index]
            converter_name: str = graph.es.select(_source=source_index, _target=target_index)[0][DB_NAME_KEY]
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
        l_possible_in_format_ids = list({self.d_uuids_from_indices[x.source] for x in l_conversion_edges})
        l_possible_out_format_ids = list({self.d_uuids_from_indices[x.target] for x in l_conversion_edges})

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
        self._d_converter_info_from_name: dict[str, ConverterInfo] | None = None
        self._d_converter_info_from_id: dict[int, ConverterInfo] | None = None
        self._l_unsorted_converter_info: list[ConverterInfo] | None = None
        self._d_format_info_from_id: dict[int, FormatInfo] | None = None
        self._d_format_info_from_name: dict[str, FormatInfo] | None = None
        self._l_unsorted_format_info: list[FormatInfo] | None = None
        self._conversions_table: ConversionsTable | None = None

    def _init_converter_info(self):
        """Initialises the private dicts and lists for converter info
        """
        self._d_converter_info_from_name: dict[str, ConverterInfo] = {}
        self._d_converter_info_from_id: dict[int, ConverterInfo] = {}
        self._l_unsorted_converter_info: list[ConverterInfo] = []
        for d_single_converter_info in self.converters:
            name: str = regularize_name(d_single_converter_info[DB_NAME_KEY])
            if name in self._d_converter_info_from_name:
                logger.warning(f"Converter '{name}' appears more than once in the database. Only the first instance"
                               " will be used.")
                continue

            single_converter_info = ConverterInfo(name=name,
                                                  parent=self,
                                                  d_single_converter_info=d_single_converter_info,
                                                  d_data=self._d_data)
            self._d_converter_info_from_name[name] = single_converter_info
            self._d_converter_info_from_id[single_converter_info.id] = single_converter_info
            self._l_unsorted_converter_info.append(single_converter_info)

    @property
    def d_converter_info_from_name(self) -> dict[str, ConverterInfo]:
        """Generate the converter info dict (indexed by name) when needed
        """
        if self._d_converter_info_from_name is None:
            self._init_converter_info()

        return self._d_converter_info_from_name

    @property
    def d_converter_info_from_id(self) -> dict[int, ConverterInfo]:
        """Generate the converter info dict (indexed by ID) when needed
        """
        if self._d_converter_info_from_id is None:
            self._init_converter_info()

        return self._d_converter_info_from_id

    @property
    def d_converter_info(self) -> dict[str, ConverterInfo]:
        """DEPRECATED: Get a dict of converter info keyed by name"""
        warnings.warn("`d_converter_info` is deprecated as of version 0.4.0 and due to be removed in a future release. "
                      "To get a ConverterInfo from the converter name (the previous functionality of this), use "
                      "`d_converter_info_from_name`. To get a FormatInfo from the format UUID, use "
                      "`d_converter_info_from_id`.", DeprecationWarning)
        return self.d_converter_info_from_name

    @property
    def l_unsorted_converter_info(self) -> list[ConverterInfo]:
        """Generate the unsorted converter info list when needed
        """
        if self._l_unsorted_converter_info is None:
            self._init_converter_info()

        return self._l_unsorted_converter_info

    @property
    def l_converter_info(self):
        """DEPRECATED: Get a list of converter info keyed by ID"""
        deprecation_msg = ("`l_converter_info` is deprecated as of version 0.4.0 and due to be removed in a future "
                           "release. To get a ConverterInfo from the format UUID, use `d_converter_info_from_id`. To "
                           "get an unsorted list of ConverterInfos, use `l_unsorted_converter_info`.")
        warnings.warn(deprecation_msg, DeprecationWarning)
        raise AttributeError(deprecation_msg)

    @property
    def d_format_info_from_name(self) -> dict[str, list[FormatInfo]]:
        """Generate the format info from format name dict when needed
        """
        if self._d_format_info_from_name is None:
            self._init_formats_and_conversions()

        return self._d_format_info_from_name

    @property
    def d_format_info_from_id(self) -> dict[int, FormatInfo]:
        """Generate the format info from format ID dict when needed
        """
        if self._d_format_info_from_id is None:
            self._init_formats_and_conversions()

        return self._d_format_info_from_id

    @property
    def d_format_info(self) -> dict[str, list[FormatInfo]]:
        """DEPRECATED: Get a dict of format info keyed by format name"""
        warnings.warn("`d_format_info` is deprecated as of version 0.4.0 and due to be removed in a future release. "
                      "To get a FormatInfo from the format name (the previous functionality of this), use "
                      "`d_format_info_from_name`. To get a FormatInfo from the format UUID, use "
                      "`d_format_info_from_id`.", DeprecationWarning)
        return self.d_format_info_from_name

    @property
    def l_unsorted_format_info(self) -> list[FormatInfo]:
        """Generate the unsorted format info list when needed
        """
        if self._l_unsorted_format_info is None:
            self._init_formats_and_conversions()

        return self._l_unsorted_format_info

    @property
    def l_format_info(self):
        """DEPRECATED: Get a list of format info keyed by ID"""
        deprecation_msg = ("`l_format_info` is deprecated as of version 0.4.0 and due to be removed in a future "
                           "release. To get a FormatInfo from the format UUID, use `d_format_info_from_id`. To get an "
                           "unsorted list of FormatInfos, use `l_unsorted_format_info`.")
        warnings.warn(deprecation_msg, DeprecationWarning)
        raise AttributeError(deprecation_msg)

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

        # Make the dict of format info keyed by ID
        self._d_format_info_from_id: dict[int, FormatInfo] = {}

        for d_single_format_info in self.formats:
            lc_name: str = d_single_format_info[DB_FORMAT_EXT_KEY]

            format_info = FormatInfo(name=lc_name,
                                     parent=self,
                                     d_single_format_info=d_single_format_info)

            self._d_format_info_from_id[format_info.id] = format_info

        # Create a temporary version of the unsorted format info list. We'll create a pruned version later, but the
        # unpruned version is needed to create the conversions table, which is needed before we can prune it
        self._l_unsorted_format_info = list(self._d_format_info_from_id.values())

        # Initialize the conversions table now
        self._conversions_table = ConversionsTable(l_converts_to=self.converts_to,
                                                   parent=self)

        # Use the conversions graph to prune any formats which have no valid conversions

        # Get a slice of the table which only includes supported converters
        supported_graph = self._conversions_table.supported_graph
        d_indices_from_uuids: dict[int, int] = self._conversions_table.d_indices_from_uuids

        l_ids_to_remove: list[int] = []
        for format_id, format_info in self._d_format_info_from_id.items():
            if not format_info or supported_graph.degree(d_indices_from_uuids[format_id]) == 0:
                # The format isn't supported for any conversions, so mark it to be removed from the dict
                # (Can't remove while we're iterating over the dict)
                l_ids_to_remove.append(format_id)
        for id in l_ids_to_remove:
            del self._d_format_info_from_id[id]

        # Now create the formats from name dict
        self._d_format_info_from_name: dict[str, list[FormatInfo]] = {}

        for format_info in self._d_format_info_from_id.values():

            lc_name = format_info.name.lower()

            # Each name may correspond to multiple formats, so we use a list for each entry to list all possible
            # formats for each name
            if lc_name not in self._d_format_info_from_name:
                self._d_format_info_from_name[lc_name] = []

            self._d_format_info_from_name[lc_name].append(format_info)

        # Finally, create a list of format infos (with arbitrary index)
        self._l_unsorted_format_info = list(self._d_format_info_from_id.values())

    def _get_converter_list(self) -> str:
        return "\n".join([f"{x.pretty_name} (ID: {x.id} / {UUID(int=x.id)})" for x in self.l_unsorted_converter_info])

    @overload
    def get_converter_info(self, converter_name_or_id: str | int | UUID | ConverterInfo) -> ConverterInfo: ...

    @overload
    def get_converter_info(self, converter_name_or_id: None) -> list[ConverterInfo]: ...

    @overload
    def get_converter_info(self) -> list[ConverterInfo]: ...

    def get_converter_info(self, converter_name_or_id: str | int | UUID | ConverterInfo | None = None) -> (
            ConverterInfo | list[ConverterInfo]):
        """Gets the information on converters or a given converter stored in the database

        Parameters
        ----------
        converter_name_or_id : str | int | UUID | ConverterInfo | None
            The name or UUID of the converter to get info for. Default None, which results in a list being returned of
            the info for all converters in the database

        Returns
        -------
        ConverterInfo | list[ConverterInfo]
            If `converter_name_or_id` is provided, will return a single `ConverterInfo` (or raise an exception if the
            name is invalid). If not provided, a list of all `ConverterInfo` objects in the database will be returned

        Raises
        ------
        FileConverterDatabaseException
        If `name` is provided but does not match the name of a converter in the database
        """

        if isinstance(converter_name_or_id, str):
            try:
                return self.d_converter_info_from_name[regularize_name(converter_name_or_id)]
            except KeyError:
                try:
                    return self.d_converter_info_from_id[UUID(converter_name_or_id).int]
                except (KeyError, ValueError):
                    raise FileConverterDatabaseException(f"Converter '{converter_name_or_id}' not found as a name in "
                                                         "the database and/or was not recognised as a value UUID. "
                                                         "Known converters are:\n" +
                                                         self._get_converter_list(),
                                                         help=True)
        elif isinstance(converter_name_or_id, int):
            try:
                return self.d_converter_info_from_id[converter_name_or_id]
            except KeyError:
                raise FileConverterDatabaseException(f"Converter ID '{converter_name_or_id}' not found in the "
                                                     "database. Known converters are:\n" +
                                                     self._get_converter_list(),
                                                     help=True)
        elif isinstance(converter_name_or_id, UUID):
            try:
                return self.d_converter_info_from_id[converter_name_or_id.int]
            except KeyError:
                raise FileConverterDatabaseException(f"Converter ID '{converter_name_or_id}' not found in the "
                                                     "database. Known converters are:\n" +
                                                     self._get_converter_list(),
                                                     help=True)
        elif isinstance(converter_name_or_id, ConverterInfo):
            # Silently return if it's already a ConverterInfo
            return converter_name_or_id
        elif converter_name_or_id is None:
            return self.l_unsorted_converter_info
        else:
            raise FileConverterDatabaseException(f"Invalid key passed to `get_converter_info`: '{converter_name_or_id}'"
                                                 f" of type '{type(converter_name_or_id)}'. Type must be `str`, "
                                                 "`int`, or `UUID`")

    @overload
    def get_format_info(self,
                        format_name_or_id: str | int | UUID | FormatInfo,
                        which: int | None = None) -> FormatInfo: ...

    @overload
    def get_format_info(self,
                        format_name_or_id: str | int | UUID | FormatInfo,
                        which: Literal["all"]) -> list[FormatInfo]: ...

    def get_format_info(self,
                        format_name_or_id: str | int | UUID | FormatInfo,
                        which: int | Literal["all"] | None = None) -> FormatInfo | list[FormatInfo]:
        """Gets the information on a given file format stored in the database

        Parameters
        ----------
        format_name_or_id : str | int | UUID | FormatInfo
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
            # Check first if it's a UUID
            try:
                format_info = self.d_format_info_from_id[UUID(format_name_or_id).int]
                if which == "all":
                    return [format_info]
                return format_info
            except KeyError:
                raise FileConverterDatabaseException(f"Format ID '{format_name_or_id}' not recognised",
                                                     help=True)
            except ValueError:
                pass

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

            l_possible_format_info = self.d_format_info_from_name.get(format_name_or_id, [])

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
                msg = (f"Extension '{format_name_or_id}' is ambiguous and must be defined by disambiguated name or ID. "
                       "Possible formats are:")
                for possible_format_info in l_possible_format_info:
                    msg += (f"\n{possible_format_info.disambiguated_name} (ID: {possible_format_info.id}): " +
                            possible_format_info.note)
                raise FileConverterDatabaseException(msg, help=True)

        elif isinstance(format_name_or_id, int):
            try:
                format_info = self.d_format_info_from_id[format_name_or_id]
            except KeyError:
                if return_as_list:
                    return []
                raise FileConverterDatabaseException(f"Format ID '{format_name_or_id}' not recognised",
                                                     help=True)

        elif isinstance(format_name_or_id, UUID):
            try:
                format_info = self.d_format_info_from_id[format_name_or_id.int]
            except KeyError:
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


def get_database_path() -> Path:
    """Get the absolute path to the database file
    """
    qualified_database_filename = get_package_path() / const.DATABASE_FILENAME
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


@overload
def get_converter_info(converter_name_or_id: str | int | UUID | ConverterInfo) -> ConverterInfo: ...


@overload
def get_converter_info(converter_name_or_id: None) -> list[ConverterInfo]: ...


@overload
def get_converter_info() -> list[ConverterInfo]: ...


def get_converter_info(converter_name_or_id: str | int | UUID | ConverterInfo | None = None) -> (
        ConverterInfo | list[ConverterInfo]):
    """Gets the information on converters or a given converter stored in the database

    Parameters
    ----------
    converter_name_or_id : str | int | UUID | ConverterInfo | None
        The name or UUID of the converter to get info for. Default None, which results in a list being returned of the
        info for all converters in the database

    Returns
    -------
    ConverterInfo | list[ConverterInfo]
        If `converter_name_or_id` is provided, will return a single `ConverterInfo` (or raise an exception if the name
        is invalid). If not provided, a list of all `ConverterInfo` objects in the database will be returned

    Raises
    ------
    FileConverterDatabaseException
        If `name` is provided but does not match the name of a converter in the database
    """

    return get_database().get_converter_info(converter_name_or_id)


@overload
def get_format_info(format_name_or_id: str | int | UUID | FormatInfo,
                    which: int | None = None) -> FormatInfo: ...


@overload
def get_format_info(format_name_or_id: str | int | UUID | FormatInfo,
                    which: Literal["all"]) -> list[FormatInfo]: ...


def get_format_info(format_name_or_id: str | int | UUID | FormatInfo,
                    which: int | Literal["all"] | None = None) -> FormatInfo | list[FormatInfo]:
    """Gets the information on a given file format stored in the database

    Parameters
    ----------
    format_name_or_id : str | int | UUID | FormatInfo
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
                           in_format: str | int | UUID | FormatInfo,
                           out_format: str | int | UUID | FormatInfo) -> ConversionQualityInfo | None:
    """Get an indication of the quality of a conversion from one format to another, or if it's not possible

    Parameters
    ----------
    converter_name : str
        The name of the converter
    in_format : str | int | UUID | FormatInfo
        The extension or ID of the input file format
    out_format : str | int | UUID | FormatInfo
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


def get_possible_conversions(in_format: str | int | UUID | FormatInfo,
                             out_format: str | int | UUID | FormatInfo) -> list[
                                 tuple[ConverterInfo, FormatInfo, FormatInfo]]:
    """Get a list of converters which can perform a conversion from one format to another and disambiguate in the case
    of ambiguous input/output formats

    Parameters
    ----------
    in_format : str | int | UUID | FormatInfo
        The extension or ID of the input file format
    out_format : str | int | UUID | FormatInfo
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


def get_conversion_pathway(in_format: str | int | UUID | FormatInfo,
                           out_format: str | int | UUID | FormatInfo,
                           only: Literal["all"] | Literal["supported"] | Literal["registered"] = "all"
                           ) -> list[tuple[ConverterInfo, FormatInfo, FormatInfo]] | None:
    """Get a list of conversions that can be performed to convert one format to another. This is primarily used when a
    direct conversion is not supported by any individual converter. Only one possible pathway will be returned,
    prioritising pathways which do not lose lose and then re-extrapolate any information stored by some formats and not
    others along the path.

    Parameters
    ----------
    in_format : str | int | UUID | FormatInfo
        The input file format. For this function, the format must be defined uniquely, either by using a disambiguated
        extension, ID, or FormatInfo
    out_format : str | int | UUID | FormatInfo
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
                         in_format: str | int | UUID | FormatInfo,
                         out_format: str | int | UUID | FormatInfo) -> tuple[FormatInfo, FormatInfo]:
    """Try to disambiguate formats by seeing if there's only one possible conversion between formats matching those
    provided.

    Parameters
    ----------
    converter_name : str
        The name of the converter
    in_format : str | int | UUID | FormatInfo
        The extension or ID of the input file format
    out_format : str | int | UUID | FormatInfo
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
        msg = (f"Conversion from {in_format} to {out_format} with converter {converter_name} is ambiguous. Please "
               "Use the ID or disambiguated name (listed below) of the desired conversion. Possible matching "
               "conversions are:\n")
        for _, possible_in_format, possible_out_format in l_possible_conversions:
            msg += (f"    {possible_in_format.id}: {possible_in_format.disambiguated_name} "
                    f"({possible_in_format.note}) to "
                    f"{possible_out_format.id}: {possible_out_format.disambiguated_name} "
                    f"({possible_out_format.note})\n")
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


def get_in_format_args(converter_name: str | int | UUID | ConverterInfo,
                       format_name: str | int | UUID | FormatInfo,
                       arg: str | None = None) -> tuple[list[FlagInfo], list[OptionInfo]] | ArgInfo:
    """Get the input flags and options supported by a given converter for a given format (provided as its extension).
    Optionally will provide information on just a single flag or option if its value is provided as an optional argument

    Parameters
    ----------
    converter_name : str | int | UUID | ConverterInfo
        The converter name
    format_name : str | int | UUID | FormatInfo
        The file format name (extension)
    arg : str | None
        If provided, only information on this flag or option will be provided

    Returns
    -------
    tuple[list[FlagInfo], list[OptionInfo]]
        (if `arg` not provided) A list of info for the allowed flags, and a list of info for the allowed options. Each
        list is sorted by the ID of the flag or option.
    ArgInfo
        (if `arg` provided) Info on the provided flag or option
    """

    converter_info = get_converter_info(converter_name)
    tl_args = converter_info.get_in_format_args(format_name)
    if not arg:
        return tl_args
    return _find_arg(tl_args, arg)


def get_out_format_args(converter_name: str | int | UUID | ConverterInfo,
                        format_name: str | int | UUID | FormatInfo,
                        arg: str | None = None) -> tuple[list[FlagInfo], list[OptionInfo]] | ArgInfo:
    """Get the output flags and options supported by a given converter for a given format (provided as its extension).
    Optionally will provide information on just a single flag or option if its value is provided as an optional argument

    Parameters
    ----------
    converter_name : str | int | UUID | ConverterInfo
        The converter
    format_name : str | int | UUID | FormatInfo
        The file format (extension)
    arg : str | None
        If provided, only information on this flag or option will be provided

    Returns
    -------
    tuple[list[FlagInfo], list[OptionInfo]]
        (if `arg` not provided) A list of info for the allowed flags, and a list of info for the allowed options. Each
        list is sorted by the ID of the flag or option.
    ArgInfo
        (if `arg` provided) Info on the provided flag or option
    """

    converter_info = get_converter_info(converter_name)
    tl_args = converter_info.get_out_format_args(format_name)
    if not arg:
        return tl_args
    return _find_arg(tl_args, arg)


def get_conversion_prop_weight(in_format_info: FormatInfo, out_format_info: FormatInfo) -> int:
    """Get the property weight for a conversion from `in_format_info` to `out_format_info` (not including the offset
    applied to it when stored in the total weight).

    Parameters
    ----------
    in_format_info : FormatInfo
        The source format for the conversion
    out_format_info : FormatInfo
        The output format for the conversion

    Returns
    -------
    int
        64-bit bit weight, where bits set to 1 indicate the properties lost or potentially in this conversion
    """

    # Start the weight as the minimum weight for any conversion. We'll turn on bits for each property potentially lost
    prop_weight = STEP_WEIGHT

    for prop, bit in D_PROP_BITS.items():
        in_prop: bool | None = getattr(in_format_info, prop)
        out_prop: bool | None = getattr(out_format_info, prop)

        # Add a weight for this conversion if the input property status is True/Unknown and output is False/Unknown, to
        # be maximally conservative
        if (in_prop is True or in_prop is None) and not out_prop:
            prop_weight |= 1 << bit

    return prop_weight


def get_conversion_prec_weight(in_format_info: FormatInfo, out_format_info: FormatInfo) -> int:
    """Get the precision weight for a conversion from `in_format_info` to `out_format_info` (not including the offset
    applied to it when stored in the total weight).

    Parameters
    ----------
    in_format_info : FormatInfo
        The source format for the conversion
    out_format_info : FormatInfo
        The output format for the conversion

    Returns
    -------
    int
        64-bit bit weight, where the bit 3N is set to 1, with N being the number of decimal places of precision lost,
        bound to 0 <= N <= 12
    """

    # Calculate the precision loss, defaulting to the maximum if unknown
    prec_loss = PREC_MAX_DIGIT_LOSS
    if in_format_info.precision is not None and out_format_info.precision is not None:
        prec_loss = in_format_info.precision - out_format_info.precision
    prec_loss = min(max(prec_loss, PREC_MIN_DIGIT_LOSS), PREC_MAX_DIGIT_LOSS)

    return 1 << PREC_GAP_BITS*prec_loss


def get_conversion_time_weight(in_format_info: FormatInfo, out_format_info: FormatInfo) -> int:
    """Get the time weight for a conversion from `in_format_info` to `out_format_info` (not including the offset
    applied to it when stored in the total weight)

    TODO: Implement properly

    Parameters
    ----------
    in_format_info : FormatInfo
        The source format for the conversion
    out_format_info : FormatInfo
        The output format for the conversion

    Returns
    -------
    int
        64-bit bit weight, representing the weight based on the conversion time (implementation TBD)
    """
    return 0


def get_conversion_weight(in_format_info: FormatInfo, out_format_info: FormatInfo) -> int:
    """Get the combined weight for a conversion

    Parameters
    ----------
    in_format_info : FormatInfo
        The source format for the conversion
    out_format_info : FormatInfo
        The output format for the conversion

    Returns
    -------
    int
        128-bit weight
    """

    return calc_conversion_weight(get_conversion_prop_weight(in_format_info, out_format_info),
                                  get_conversion_prec_weight(in_format_info, out_format_info),
                                  get_conversion_time_weight(in_format_info, out_format_info))


def calc_conversion_weight(prop_weight: int, prec_weight: int, time_weight: int):
    """Calculate the combined weight for a conversion from its component weights. The weights must be in the provided
    range, or else the output will have undefined behaviour

    Parameters
    ----------
    prop_weight : int
        The conversion property weight, in the range 0 < prop_weight < 2**64
    prec_weight : int
        The conversion precision weight, in the range 0 < prop_weight < 2**48
    time_weight : int
        The conversion time weight, in the range 0 < prop_weight < 2**16

    Returns
    -------
    int
        The combined conversion weight
    """
    return ((prop_weight << PROP_WEIGHT_BIT_OFFSET) +
            (prec_weight << PREC_WEIGHT_BIT_OFFSET) +
            (time_weight << TIME_WEIGHT_BIT_OFFSET))


class ConversionWeightParts(NamedTuple):
    prop_weight: int
    prec_weight: int
    time_weight: int


def split_conversion_weight(conversion_weight: int):
    """Splits the total conversion weight into the parts for each component weight

    Parameters
    ----------
    conversion_weight : int
        The total conversion weight

    Returns
    -------
    ConversionWeightParts
        NamedTuple of prop_weight, prec_weight, and time_weight
    """

    prop_weight = conversion_weight >> PROP_WEIGHT_BIT_OFFSET
    conversion_weight -= prop_weight << PROP_WEIGHT_BIT_OFFSET

    prec_weight = conversion_weight >> PREC_WEIGHT_BIT_OFFSET
    conversion_weight -= prec_weight << PREC_WEIGHT_BIT_OFFSET

    time_weight = conversion_weight >> TIME_WEIGHT_BIT_OFFSET

    return ConversionWeightParts(prop_weight, prec_weight, time_weight)
