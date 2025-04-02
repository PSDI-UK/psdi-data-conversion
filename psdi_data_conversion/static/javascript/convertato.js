/*
  convertato.js
  Version 1.0, 8th November 2024

  This is the JavaScript which makes the convertato.htm gui work.
*/

import { commonConvertReady, convertFile, getExtCheck, splitArchiveExt, isArchiveExt } from "./convert_common.js"

var token = "",
    in_ext = "",
    out_ext = "",
    in_str = "",
    out_str = "";

$(document).ready(function () {
    [token, in_str, in_ext, out_str, out_ext] = commonConvertReady("Atomsk");
    $("#uploadButton").click(submitFile);
});

// Uploads a user-supplied file
function submitFile() {
    const file = $("#fileToUpload")[0].files[0],
        [fname, ext] = splitArchiveExt(file.name);

    var quality = sessionStorage.getItem("success"),
        start = quality.indexOf(':') + 2,
        finish = quality.lastIndexOf('(') - 1;

    quality = quality.substring(start, finish);

    const read_flags_text = $("#inFlags").find(":selected").text(),
        read_flags = '';

    const write_flags_text = $("#outFlags").find(":selected").text(),
        write_flags = '';

    var count = 0,
        read_arg_flags = '',
        write_arg_flags = '',
        read_args = '',
        write_args = '',
        all_args_entered = true;

    const checked_in = $('input[name=in_arg_check]:checked'),
        checked_out = $('input[name=out_arg_check]:checked');

    checked_in.each(function () {
        read_arg_flags += $("#" + this.id).val()[0];
        const arg = $("#in_arg_text" + this.id.substring(this.id.length - 1, this.id.length)).val();

        if (/\S/.test(arg)) {
            read_args += arg.trim() + '£';
        }
        else {
            all_args_entered = false;
        }
    })

    checked_out.each(function () {
        write_arg_flags += $("#" + this.id).val()[0];
        const arg = $("#out_arg_text" + this.id.substring(this.id.length - 1, this.id.length)).val();

        if (/\S/.test(arg)) {
            write_args += arg.trim() + '£';
        }
        else {
            all_args_entered = false;
        }
    })

    if (!all_args_entered) {
        alert('All ticked option flags need additional information to be entered into the associated text box.');
        return;
    }

    const coordinates = 'neither', //$('input[name="coordinates"]:checked').val(),
        coordOption = 'medium'; //$('input[name="coordOptions"]:checked').val(),

    // If the uploaded file was an archive, we'll expect to download one too. Otherwise we'll expect to download
    // something with the output extension
    const download_fname = isArchiveExt(ext) ? fname + "-" + out_ext + "." + ext : fname + "." + out_ext;

    var form_data = new FormData();

    form_data.append("token", token);
    form_data.append("converter", 'Atomsk');
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