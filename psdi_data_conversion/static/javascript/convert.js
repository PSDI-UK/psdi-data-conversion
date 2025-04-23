/*
  convert.js
  Version 1.0, 24th January 2025

  This is the JavaScript which makes the convert.htm gui work.
*/


import { getInputFlags, getOutputFlags, getInputArgFlags, getOutputArgFlags } from "./data.js";
import {
    SAFE_CHAR_REGEX, commonConvertReady, convertFile, getExtCheck, splitArchiveExt, isArchiveExt
} from "./convert_common.js"

var token = "",
    in_ext = "",
    out_ext = "",
    in_str = "",
    out_str = "";

$(document).ready(function () {
    [token, in_str, in_ext, out_str, out_ext] = commonConvertReady("Open Babel");

    $('input[name="coordinates"]').change(coordOptionAvailability);
    $("#uploadButton").click(submitFile);
    $("#inFlags").click(blurBox);
    $("#outFlags").click(blurBox);

    getFlags("in", in_str);
    getFlags("out", out_str);
    getFlags("in_arg", in_str);
    getFlags("out_arg", out_str);
});

// On clicking on an option flag, the select box loses focus so that the light font and dark background (Chrome) appears immediately.
function blurBox(event) {
    this.blur();
}

// On ticking a checkbox, a text box for entry of an option flag argument appears next to it. On unticking, the text box disappears.
function enterArgument(event) {
    var arg_id = this.id.replace('check', 'text'),
        arg_label_id = this.id.replace('check', 'label');

    if ($('#' + this.id).is(':checked')) {
        // Show appropriate text box and its label
        $('#' + arg_id).show();
        $('#' + arg_label_id).show();
    }
    else {
        // Hide appropriate text box (empty) and its label
        $('#' + arg_id).val('').hide();
        $('#' + arg_label_id).hide();
    }
}


/**
 * Validate the input for a format option - if invalid, display the error message, if valid, hide it
 *
 * @param {*} event 
 */
function validateInput(event) {
    var err_id = this.id.replace('text', 'err')
    if (this.validity.patternMismatch) {
        $('#' + err_id).css({ display: "block" })
    } else {
        $('#' + err_id).css({ display: "none" })
    }
}

// Uploads a user-supplied file
function submitFile() {

    const file = $("#fileToUpload")[0].files[0],
        [fname, ext] = splitArchiveExt(file.name);


    var quality = sessionStorage.getItem("success"),
        start = quality.indexOf(':') + 2,
        finish = quality.lastIndexOf('(') - 1;

    quality = quality.substring(start, finish);

    const read_flags_text = $("#inFlags").find(":selected").text(),
        read_flags = extractFlags(read_flags_text);

    const write_flags_text = $("#outFlags").find(":selected").text(),
        write_flags = extractFlags(write_flags_text);

    var read_arg_flags = '',
        write_arg_flags = '',
        read_args = '',
        write_args = '',
        all_args_entered = true;

    const checked_in = $('input[name=in_arg_check]:checked'),
        checked_out = $('input[name=out_arg_check]:checked');

    let security_passed = true;

    checked_in.each(function () {
        read_arg_flags += $("#" + this.id).val()[0];
        const e = $("#in_arg_text" + this.id.substring(this.id.length - 1, this.id.length));
        const arg = e.val();

        if (e[0].validity.patternMismatch) {
            security_passed = false;
        }

        if (/\S/.test(arg)) {
            read_args += arg.trim() + '£';
        }
        else {
            all_args_entered = false;
        }
    })

    checked_out.each(function () {
        write_arg_flags += $("#" + this.id).val()[0];
        const e = $("#out_arg_text" + this.id.substring(this.id.length - 1, this.id.length));
        const arg = e.val();

        if (e[0].validity.patternMismatch) {
            security_passed = false;
        }

        if (/\S/.test(arg)) {
            write_args += arg.trim() + '£';
        }
        else {
            all_args_entered = false;
        }
    })

    let alert_msg = '';

    if (!security_passed) {
        alert_msg += 'ERROR: One or more ticked options contains invalid characters. They must match the regex /' +
            SAFE_CHAR_REGEX + '/ .\n';
    }

    if (!all_args_entered) {
        alert_msg += 'ERROR: All ticked options need additional information to be entered into the associated ' +
            'text box.\n';
    }

    if (alert_msg) {
        alert(alert_msg);
        return
    }

    const coordinates = $('input[name="coordinates"]:checked').val(),
        coordOption = $('input[name="coordOptions"]:checked').val();

    // If the uploaded file was an archive, we'll expect to download one too. Otherwise we'll expect to download
    // something with the output extension
    const download_fname = isArchiveExt(ext) ? fname + "-" + out_ext + "." + ext : fname + "." + out_ext;

    var form_data = new FormData();

    form_data.append("token", token);
    form_data.append("converter", 'Open Babel');
    form_data.append("from", in_ext);
    form_data.append("to", out_ext);
    form_data.append("from_full", sessionStorage.getItem("in_str"));
    form_data.append("to_full", sessionStorage.getItem("out_str"));
    form_data.append("success", quality);
    form_data.append("from_flags", read_flags);
    form_data.append("to_flags", write_flags);
    form_data.append("from_arg_flags", read_arg_flags);
    form_data.append("from_args", read_args);
    form_data.append("to_arg_flags", write_arg_flags);
    form_data.append("to_args", write_args);
    form_data.append("coordinates", coordinates);
    form_data.append("coordOption", coordOption);
    form_data.append("fileToUpload", file);
    form_data.append("upload_file", true);
    form_data.append("check_ext", getExtCheck());

    convertFile(form_data, download_fname, fname);
}

// Retrieves option flags from selected text
function extractFlags(flags_text) {
    var flags = "",
        regex = /: /g,
        match = "";

    while ((match = regex.exec(flags_text)) != null) {
        flags += flags_text[match.index - 1];
    }

    return flags;
}

// Retrieves read or write option flags associated with a file format
function getFlags(type, str) {

    try {
        const [ext, note] = str.split(": ");

        if (type == "in") {
            getInputFlags(ext, note).then((flags) => {
                populateFlagBox(flags, type);
            });
        }
        else if (type === "out") {
            getOutputFlags(ext, note).then((flags) => {
                populateFlagBox(flags, type);
            });
        }
        else if (type == "in_arg") {

            const in_arg_str_array = str.split(": "),
                in_arg_ext = in_arg_str_array[0],          // e.g. "ins"
                in_arg_note = in_arg_str_array[1];         // e.g. "ShelX"

            getInputArgFlags(in_arg_ext, in_arg_note).then((argFlags) => {
                addCheckboxes(argFlags, "in_arg");
            });
        }
        else if (type == "out_arg") {

            const out_arg_str_array = str.split(": "),
                out_arg_ext = out_arg_str_array[0],          // e.g. "ins"
                out_arg_note = out_arg_str_array[1];         // e.g. "ShelX"

            getOutputArgFlags(out_arg_ext, out_arg_note).then((argFlags) => {
                addCheckboxes(argFlags, "out_arg");
            });
        }

        return true;
    }
    catch (e) {
        return false;
    }
}

// Adds checkboxes for read or write option flags requiring an argument
function addCheckboxes(argFlags, type) {

    var container = $(`#${type}Flags`),
        flagCount = 0;

    if (argFlags.length > 0) {

        $(`#${type}Label`).show();

        for (const argFlag of argFlags) {

            const flag = argFlag.flag;
            const brief = argFlag.brief.replace(/^N\/A$/, "");
            const description = argFlag.description.replace(/^N\/A$/, "");
            const furtherInfo = argFlag.further_info.replace(/^N\/A$/, "");

            container.append(`
                <tr>
                    <td><input type='checkbox' id="${type}_check${flagCount}" name=${type}_check value="${flag}"></input></td>
                    <td><label for="${type}_check${flagCount}">${flag} [${brief}]: ${description}<label></td>
                    <td><input type='text' id="${type}_text${flagCount}" placeholder='-- type info. here --'
                         pattern='` + SAFE_CHAR_REGEX + `'></input>
                         <p class="init-hidden" id="${type}_err${flagCount}"><strong>ERROR:</strong> Input contains
                         invalid characters; it must match the regex 
                         <code class="secondary">/` + SAFE_CHAR_REGEX + `/</code></p>
                    </td>
                    <td><span id= ${type}_label${flagCount}>${furtherInfo}</span></td>
                </tr>`);

            $(`#${type}_text${flagCount}`).hide();
            $(`#${type}_text${flagCount}`).on('input', validateInput);
            $(`#${type}_label${flagCount}`).hide();
            $(`#${type}_check${flagCount}`).change(enterArgument);

            flagCount++;
        }

        container.append(`<br>`);
    }
    else {
        $(`#${type}Label`).hide();

        if (type == 'in_arg') {
            $("#flag_break").hide();
        }
    }
}

// Populates a read or write option flag box
function populateFlagBox(entries, type) {

    const el = $("#" + type + "Flags");
    const disp = $("#" + type + "FlagList");
    const flagInfo = $("#" + type + "FlagInfo");

    let infoLines = [];

    if (entries.length != 0) {

        disp.css({ display: "inline" });

        for (const entry of entries) {

            el.append(new Option(`${entry.flag}: ${entry.description}`));

            const info = `${entry.flag}: ${entry.further_info}`;

            if (!info.match(/.: N\/A/)) {
                infoLines.push(info);
            }
        }

    } else {

        $("#" + type + "_label").hide();
        $("#" + type + "_flag_break").hide();
        el.hide();
    }

    el.append(new Option(""));

    for (const infoLine of infoLines) {

        const div = $("<div>");

        div.text(infoLine);

        flagInfo.append(div);
    }
}

// Disable coordinate options if calculation type is 'neither,' otherwise enable
function coordOptionAvailability(event) {
    const calcType = $('input[name="coordinates"]:checked').val();

    if (calcType == 'neither') {
        $('input[name="coordOptions"]').prop({ disabled: true });
    }
    else {
        $('input[name="coordOptions"]').prop({ disabled: false });
    }
}

