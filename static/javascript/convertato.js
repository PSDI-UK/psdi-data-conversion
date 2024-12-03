/*
  convertato.js
  Version 1.0, 8th November 2024

  This is the JavaScript which makes the convertato.htm gui work.
*/

import {checkFile} from "./convert.js"

const fromList = new Array(),
      toList = new Array();

var token = "",
    last_select = "",
    in_ext = "",
    out_ext = "";

$(document).ready(function() {
    token = sessionStorage.getItem("token");

    const in_str = sessionStorage.getItem("in_str"),
          out_str = sessionStorage.getItem("out_str");

    // $$$$$ FUNCTION FOR THIS? $$$$$
    const in_str_array = in_str.split(": ");
    in_ext = in_str_array[0];                // e.g. "ins"
    const in_note = in_str_array[1];         // e.g. "ShelX"

    const out_str_array = out_str.split(": ");
    out_ext = out_str_array[0];
    const out_note = out_str_array[1];

    $("#heading").html("Convert from \'" + in_ext + "\' (" + in_note + ") to \'" + out_ext + "\' (" + out_note + ") using Atomsk");

    $("#fileToUpload").change(checkFile);
    $("#uploadButton").click(submitFile);
});

// $$$$$$$$$$ Retained in case of future need to write to a log file $$$$$$$$$$
// Writes user input to a server-side file
function writeLog(message) {
    var jqXHR = $.get(`/data/`, {
            'token': token,
            'data': message
        })
        .fail(function(e) {
            // For debugging
            console.log("Error writing to log");
            console.log(e.status);
            console.log(e.responseText);
        })
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
        $("#uploadButton").css({"background-color": "var(--psdi-bg-color-secondary)", "color": "gray"});
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

    checked_in.each(function() {
        read_arg_flags += $("#" + this.id).val()[0];
        const arg = $("#in_arg_text" + this.id.substring(this.id.length - 1, this.id.length)).val();

        if (/\S/.test(arg)) {
            read_args += arg.trim() + '£';
        }
        else {
            all_args_entered = false;
        }
    })

    checked_out.each(function() {
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

    convertFile(form_data, download_fname, fname);
}

// Converts user-supplied file to another format and downloads the resulting file
function convertFile(form_data, download_fname, fname) {
    var jqXHR = $.ajax({
            url: `/convert/`,
            type: "POST",
            data: form_data,
            processData: false,
            contentType: false,
            success: async function() {
                const delay = ms => new Promise(response => setTimeout(response, ms));

                downloadFile(`../downloads/${fname}.log.txt`, fname + '.log.txt')
                await delay(300);
                downloadFile(`../downloads/${download_fname}`, download_fname)
                await delay(300);

                var fdata = new FormData();

                fdata.append("filename", download_fname);
                fdata.append("logname", fname + '.log.txt');

                $.ajax({
                    url: `/delete/`,
                    type: "POST",
                    data: fdata,
                    processData: false,
                    contentType: false
                })
                .fail(function(e) {
                    // For debugging
                    console.log("Error deleting remote files after download");
                    console.log(e.status);
                    console.log(e.responseText);
                })
            },
            error: function(data) {
                //alert("ajax error, FormData: " + data);
                }
            })
            .done(response => {
                alert("To the best of our knowledge, this conversion has worked. Your output file should download automatically " +
                      "when you close this alert. Please report any problems by clicking on 'Contact' in the navigation bar.");
            })
            .fail(function(e) {
                alert("This conversion has failed. Please provide feedback on the conversion " +
                      "that you were attempting by clicking on 'Contact' in the navigation bar.");

                // For debugging
                console.log("Error converting file");
                console.log(e.status);
                console.log(e.responseText);
            })
}

// A link is created, clicked and removed, resulting in the download of a file
function downloadFile(path, filename) {
    const a = $("<a>")
          .attr("href", path)
          .attr("download", filename)
          .appendTo("body");
    a[0].click();
    a.remove();
}
