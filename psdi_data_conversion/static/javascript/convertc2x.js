/*
  convertato.js
  Version 1.0, 10th January 2025

  This is the JavaScript which makes the convertc2x.htm gui work.
*/

import { commonConvertReady, convertFile, getExtCheck } from "./convert_common.js"

var token = "",
    max_file_size = 0,
    in_ext = "",
    out_ext = "",
    in_str = "",
    out_str = "";

$(document).ready(function () {
    [token, max_file_size, in_str, in_ext, out_str, out_ext] = commonConvertReady("c2x");
    $("#uploadButton").click(submitFile);
});

// Uploads a user-supplied file
function submitFile() {
    const file = $("#fileToUpload")[0].files[0],
        fname = file.name.split(".")[0],
        extension = file.name.split(".")[1];

    var quality = sessionStorage.getItem("success"),
        start = quality.indexOf(':') + 2,
        finish = quality.lastIndexOf('(') - 1;

    quality = quality.substring(start, finish);

    if (extension != in_ext) {
        alert("The file extension is not " + in_ext + ": please select another file or change the 'from' format on the 'Home' page.");
        $("#uploadButton").css({ "background-color": "var(--psdi-bg-color-secondary)", "color": "gray" });
        return;
    }

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
        coordOption = 'medium', //$('input[name="coordOptions"]:checked').val(),
        download_fname = file.name.split(".")[0] + "." + out_ext;

    var form_data = new FormData();

    form_data.append("token", token);
    form_data.append("converter", 'c2x');
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
