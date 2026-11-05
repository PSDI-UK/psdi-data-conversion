from http import HTTPStatus

from flask import jsonify, abort, send_file

from psdi_data_conversion.database import (FormatInfo, FlagInfo, OptionInfo, get_conversion_pathway, get_conversion_quality,
                                           get_converter_info, get_format_info, get_in_format_args,
                                           get_out_format_args, get_possible_conversions, get_possible_formats)

from psdi_data_conversion.database import (get_database, FileConverterDatabaseException)
from psdi_data_conversion.converter import run_converter

import uuid
import os

from tempfile import NamedTemporaryFile


def paginate_list(items, start, num):

    return items


def index_converter():

    all_converters = get_converter_info()

    result = []

    for converter in all_converters:

        result.append({
            "name": converter.name,
            "id": converter.uuid,
            "description": converter.description,
            "url": converter.url,
        })

    return jsonify(result)


def read_converter(converter_id):

    try:

        converter = get_database().get_converter_info(converter_id)

        return jsonify({
            "name": converter.name,
            "id": converter.uuid,
            "description": converter.description,
            "url": converter.url,
        })

    except FileConverterDatabaseException:

        return jsonify({"error": f"Converter {converter_id} was not found"}), HTTPStatus.NOT_FOUND


def read_converter_formats(converter_id):

    try:

        converter = get_database().get_converter_info(converter_id)

        all_formats = get_possible_formats(converter.name)

        input_formats = []
        output_formats = []

        input_format_names = all_formats[0]
        output_format_names = all_formats[1]

        for format_name in input_format_names:

            format_info = get_format_info(format_name)

            input_formats.append({
                "id": format_info.uuid,
                "label": format_info.note
            })

        for format_name in output_format_names:

            format_info = get_format_info(format_name)
            # print(f"get_out_format_args = ${get_out_format_args(converter, format_name, None)}")
            output_formats.append({
                "id": format_info.uuid,
                "label": format_info.note
            })

        return jsonify({
            "inputFormats": input_formats,
            "outputFormats": output_formats
        })


    except FileConverterDatabaseException:

        return jsonify({"error": f"Converter {converter_id} was not found"}), HTTPStatus.NOT_FOUND


def describe_flag(flag: FlagInfo):

    return {
        "id": flag.uuid,
        "flag": flag.flag,
        "description": flag.description,
        "info": flag.info
    }


def describe_option(option: OptionInfo):

    return {
        "id": option.uuid,
        "flag": option.flag,
        "description": option.description,
        "info": option.info,
        "brief": option.brief,
    }


def read_converter_format(converter_id, format_id):

    converter = get_database().get_converter_info(converter_id)

    if converter is None:
        return jsonify({"error": f"Converter {converter_id} was not found"}), HTTPStatus.NOT_FOUND

    format_info = get_format_info(format_id)

    if format_info is None:
        return jsonify({"error": f"Format {format_id} was not found"}), HTTPStatus.NOT_FOUND

    in_flags, in_options = get_in_format_args(converter, format_info, None)
    out_flags, out_options = get_out_format_args(converter, format_info, None)

    return jsonify({
        "inputFlags": [describe_flag(flag) for flag in in_flags],
        "inputOptions": [describe_option(option) for option in in_options],
        "outputFlags": [describe_flag(flag) for flag in out_flags],
        "outputOptions": [describe_option(option) for option in out_options],
    })


def create_format_entry(format_info):

    entry = { "id": format_info.uuid, "name": format_info.name }

    if (format_info.two_dim != None):
        entry |= { "twoDim": format_info.two_dim }

    if (format_info.three_dim != None):
        entry |= { "threeDim": format_info.three_dim }

    return entry


# def index_format():

#     all_formats = get_format_info()

#     result = []

#     for format in all_formats:

#         result.append({
#             "name": format.name,
#             "id": format.uuid
#         })

#     return jsonify(result)


def read_format(format_id):

    try:
        format_info = get_format_info(format_id)
    except FileConverterDatabaseException:
        return jsonify({"error": f"Format {format_id} was not found"}), HTTPStatus.NOT_FOUND

    if format_info is None:
        return jsonify({"error": f"Format {format_id} was not found"}), HTTPStatus.NOT_FOUND

    return jsonify(create_format_entry(format_info))


def get_conversion_information(converter_id, from_format, to_format):

    try:
        converter = get_database().get_converter_info(converter_id)
    except FileConverterDatabaseException:
        return jsonify({"error": f"Converter {converter_id} was not found"}), HTTPStatus.NOT_FOUND

    try:
        from_format_info = get_format_info(from_format)
    except FileConverterDatabaseException:
        return jsonify({"error": f"Format {from_format} was not found"}), HTTPStatus.NOT_FOUND

    try:
        to_format_info = get_format_info(to_format)
    except FileConverterDatabaseException:
        return jsonify({"error": f"Format {to_format} was not found"}), HTTPStatus.NOT_FOUND

    conversion_quality = get_conversion_quality(converter, from_format_info, to_format_info)

    return jsonify({
        "quality": conversion_quality.qual_str,
        "from": {
            "compositionSupported": conversion_quality.d_prop_conversion_info["composition"].input_supported,
            "connectionsSupported": conversion_quality.d_prop_conversion_info["connections"].input_supported,
            "twoDimSupported": conversion_quality.d_prop_conversion_info["two_dim"].input_supported,
            "threeDimSupported": conversion_quality.d_prop_conversion_info["three_dim"].input_supported
        },
        "to": {
            "compositionSupported": conversion_quality.d_prop_conversion_info["composition"].output_supported,
            "connectionsSupported": conversion_quality.d_prop_conversion_info["connections"].output_supported,
            "twoDimSupported": conversion_quality.d_prop_conversion_info["two_dim"].output_supported,
            "threeDimSupported": conversion_quality.d_prop_conversion_info["three_dim"].output_supported
        }
    })


def perform_conversion(converter_id, from_format, to_format, input_file, body, **kwargs):

    try:
        converter = get_database().get_converter_info(converter_id)
    except FileConverterDatabaseException:
        return (
            {"error": f"Converter {converter_id} was not found"},
            HTTPStatus.NOT_FOUND,
            { "Content-Type": "application/json" }
        )

    try:
        from_format_info = get_format_info(from_format)
    except FileConverterDatabaseException:
        return (
            {"error": f"Format {from_format} was not found"},
            HTTPStatus.NOT_FOUND,
            { "Content-Type": "application/json" }
        )

    try:
        to_format_info = get_format_info(to_format)
    except FileConverterDatabaseException:
        return (
            {"error": f"Format {to_format} was not found"},
            HTTPStatus.NOT_FOUND,
            { "Content-Type": "application/json" }
        )

    in_flags, in_options = get_in_format_args(converter, from_format_info, None)
    out_flags, out_options = get_out_format_args(converter, to_format_info, None)

    # Build strings for run_converter.

    from_arg_flags = ""
    from_args = ""
    from_flags = ""
    to_arg_flags = ""
    to_args = ""
    to_flags = ""

    for in_flag in in_flags:
        if (in_flag.flag in body["args"]) or (in_flag.uuid in body["args"]):
            to_flags += in_flag.flag

    for in_option in in_options:

        arg = None

        if in_option.flag in body["args"]:
            arg = body["args"][in_option.flag]

        if in_option.uuid in body["args"]:
            arg = body["args"][in_option.uuid]

        if arg != None and isinstance(arg, str):
            arg = arg.strip()
            if arg and "£" not in arg:
                from_arg_flags += in_option.flag
                from_args += f"{arg}£"

    for out_flag in out_flags:
        if (out_flag.flag in body["args"]) or (out_flag.uuid in body["args"]):
            to_flags += out_flag.flag

    for out_option in out_options:

        arg = None

        if out_option.flag in body["args"]:
            arg = body["args"][out_option.flag]

        if out_option.uuid in body["args"]:
            arg = body["args"][out_option.uuid]

        if arg != None and isinstance(arg, str):
            arg = arg.strip()
            if arg and "£" not in arg:
                to_arg_flags += out_option.flag
                to_args += f"{arg}£"

    extension = from_format_info.d_alias_exts[from_format_info.id]

    with NamedTemporaryFile(suffix=f".{extension}", delete_on_close=False) as infile:

        input_file.save(infile)
        infile.close()

        conversion_result = run_converter(
            filename = infile.name,
            name = converter.name,
            from_format = from_format_info.id,
            to_format = to_format_info.id,
            input_dir = "./psdi_data_conversion/static/uploads",
            output_dir = "./psdi_data_conversion/static/downloads",
            max_file_size = None,
            log_file = None,
            log_mode = "full",
            strict = True,
            permission_level = 2,
            archive_output = False,
            data = {
                # "coordinates": "neither",
                # "coordOption": "medium",
                "from_arg_flags": from_arg_flags,
                "from_args": from_args,
                "from_flags": from_flags,
                "to_arg_flags": to_arg_flags,
                "to_args": to_args,
                "to_flags": to_flags,
            },
            log_level = None,
            delete_input = True,
            abort_callback = abort
        )

        return send_file(str(os.path.abspath(conversion_result.output_filename)))
        os.unlink(conversion_result.output_filename)
