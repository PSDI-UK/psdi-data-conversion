/*
  convert.js
  Version 1.0, 8th November 2024

  This is the JavaScript which makes the convert.htm gui work.
*/

import { getInputFlags, getOutputFlags, getInputArgFlags, getOutputArgFlags } from "./data.js";

const MEGABYTE = 1024*1024;
const MAX_FILESIZE = 1*MEGABYTE;

const fromList = new Array(),
      toList = new Array();

var token = "",
    last_select = "",
    in_ext = "",
    out_ext = "",
    in_str = "",
    out_str = "";

export function commonConvertReady(converter) {
    token = sessionStorage.getItem("token");

    in_str = sessionStorage.getItem("in_str");
    out_str = sessionStorage.getItem("out_str");

    // $$$$$ FUNCTION FOR THIS? $$$$$
    const in_str_array = in_str.split(": ");
    in_ext = in_str_array[0];                // e.g. "ins"
    const in_note = in_str_array[1];         // e.g. "ShelX"

    const out_str_array = out_str.split(": ");
    out_ext = out_str_array[0];
    const out_note = out_str_array[1];

    $("#heading").html("Convert from \'" + in_ext + "\' (" + in_note + ") to \'" + out_ext + "\' (" + out_note + 
        ") using " + converter);

    $("#fileToUpload").change(checkFile);

    return [token, in_str, in_ext, out_str, out_ext];
}

$(document).ready(function() {
    [token, in_str, in_ext, out_str, out_ext] = commonConvertReady("Open Babel");

    $('input[name="coordinates"]').change(coordOptionAvailability);
    $("#uploadButton").click(submitFile);

    getFlags("in", in_str);
    getFlags("out", out_str);
    getFlags("in_arg", in_str);
    getFlags("out_arg", out_str);
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
    var //flags_text = $('#' + this.id).val(),
        arg_id = this.id.replace('check', 'text'),
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
          read_flags = extractFlags(read_flags_text);

    const write_flags_text = $("#outFlags").find(":selected").text(),
          write_flags = extractFlags(write_flags_text);

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

    const coordinates = $('input[name="coordinates"]:checked').val(),
          coordOption = $('input[name="coordOptions"]:checked').val(),
          download_fname = file.name.split(".")[0] + "." + out_ext;

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

// Retrieves option flags requiring arguments from selected text
function extractArgFlags(flags_text) {
    var flags = "",
        regex = /[.]/g,
        match = "";

    if (flags_text.length > 0)
        flags += flags_text[0];

    while ((match = regex.exec(flags_text)) != null && match.index < flags_text.length - 1) {
        flags += flags_text[match.index + 1];
    }

    return flags;
}

// Converts user-supplied file to another format and downloads the resulting file
export function convertFile(form_data, download_fname, fname) {
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
                let errLog = '/static/downloads/' + fname + '-' + download_fname + ".err";

                fetch(errLog, {cache: "no-store"})
                .then(function (response) {
                    if (response.status==404) {
                        return "An unknown error occurred, which produced no error log. Please provide feedback on " +
                            "the conversion that you were attempting by clicking on 'Contact' in the navigation bar.";
                    } else {
                        return response.text();
                    }
                })
                .then(function (text) {
                    if (text!="")
                        alert(text);
                })

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

// Retrieves read or write option flags associated with a file format
function getFlags (type, str) {

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
                    <td><input type='text' id=${type}_text${flagCount} placeholder='-- type info. here --'></input></td>
                    <td><span id= ${type}_label${flagCount}>${furtherInfo}</span></td>
                </tr>`);

            $(`#${type}_text${flagCount}`).hide();
            $(`#${type}_label${flagCount}`).hide();
            $(`#${type}_check${flagCount}`).change(enterArgument);

            flagCount++;
        }
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

        disp.css({display: "inline"});

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

        const p = $("<p>");

        p.text(infoLine);

        flagInfo.append(p);
    }
}

// Disable coordinate options if calculation type is 'neither,' otherwise enable
function coordOptionAvailability(event) {
    const calcType = $('input[name="coordinates"]:checked').val();

    if (calcType == 'neither') {
        $('input[name="coordOptions"]').prop({disabled: true}); 
    }       
    else {
        $('input[name="coordOptions"]').prop({disabled: false}); 
    }
}

// Check that the file meets requirements for upload
function checkFile(event) {

    let allGood = true;
    let file = this.files[0];
    let message = "";
    
    // Check file has the proper extension
    const file_name = file.name;
    const file_name_array = file_name.split(".");
    const extension = file_name_array[1];
    if (extension != in_ext) {
        message += "The file extension is not " + in_ext +
            ": please select another file or change the 'from' format on the 'Home' page.";
        allGood = false;
    }

    // Check file does not exceed maximum size
    if (file.size > MAX_FILESIZE) {
        if (message!=="")
            message += "\n\n";
        message += "The file exceeds the maximum size limit of " + (MAX_FILESIZE/MEGABYTE).toFixed(2) +
          " MB; its size is " + (file.size/MEGABYTE).toFixed(2) + " MB.";
        allGood = false;
    }

    if(allGood) {
        $("#uploadButton").css({"background-color": "var(--ifm-color-primary)", "color": "var(--ifm-hero-text-color)"});
        $("#uploadButton").prop({disabled: false});
    } else {
        $("#uploadButton").css({"background-color": "var(--psdi-bg-color-secondary)", "color": "gray"});
        $("#uploadButton").prop({disabled: true});
        alert(message);
    }
}

