/*
  convert.js
  Version 1.0, 9th September 2024

  This is the JavaScript which makes the convert.htm gui work.
*/

var token = "",
    fromList = new Array(),
    toList = new Array(),
    last_select = "",
    in_ext = "",
    out_ext = "";

$(document).ready(function() {
    token = sessionStorage.getItem("token");

    const font = sessionStorage.getItem("font"),
          size = sessionStorage.getItem("size"),
          weight = sessionStorage.getItem("weight"),
          letter = sessionStorage.getItem("letter"),
          line = sessionStorage.getItem("line"),
          colour = sessionStorage.getItem("colour"),
          back = sessionStorage.getItem("back");

    if (font != null) {
        $(".normalText, .middle, h3, #uploadButton").css({
            fontFamily: font,
            fontSize: size,
            fontWeight: weight,
            letterSpacing: letter
        });

        $(".normalText, .middle").css({lineHeight: line});
        $(".normalText").css({color: colour});
        $("h1, h2").css({letterSpacing: letter});
        $("h1, h3").css({color: colour});
        $("h3").css({fontSize: Number(size.substring(0, 2)) + 4 + 'px'});
        $("form, #upper, #inFlags, #outFlags").css({background: back});
    }

    const in_str = sessionStorage.getItem("in_str"),
          out_str = sessionStorage.getItem("out_str");

    // $$$$$ FUNCTION FOR THIS? $$$$$
    const in_str_array = in_str.split(": ");
    in_ext = in_str_array[0];                // e.g. "ins"
    const in_note = in_str_array[1];         // e.g. "ShelX"

    const out_str_array = out_str.split(": ");
    out_ext = out_str_array[0];
    const out_note = out_str_array[1];

    $("#heading").html("Convert from \'" + in_ext + "\' (" + in_note + ") to \'" + out_ext + "\' (" + out_note + ") using Open Babel");

    getFlags("in", in_str);
    getFlags("out", out_str);

    $('input[name="coordinates"]').change(coordOptionAvailability);
    $("#fileToUpload").change(checkExtension);
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

// $$$$$$$$$$ write separate function for debugging $$$$$$$$$$$

// Queries the PostgreSQL database
function queryDatabase(query, sel, callback) {
    var jqXHR = $.post(`/query/`, {
            'token': token,
            'data': query
        })
        .done(response => {
            callback(response, sel);
        })
        .fail(function(e) {
            // For debugging
            console.log("Error querying database");
            console.log(e.status);
            console.log(e.responseText);
        })
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
        $("#uploadButton").css({"background-color": "#e5e1e6", "color": "gray"});
        return;
    }

    const read_flags_text = $("#inFlags").find(":selected").text(),
          read_flags = extractFlags(read_flags_text);

    const write_flags_text = $("#outFlags").find(":selected").text(),
          write_flags = extractFlags(write_flags_text);

    const coordinates = $('input[name="coordinates"]:checked').val();
    const coordOption = $('input[name="coordOptions"]:checked').val();

    const download_fname = file.name.split(".")[0] + "." + out_ext;

    var form_data = new FormData();

    form_data.append("token", token);
    form_data.append("from", in_ext);
    form_data.append("to", out_ext);
    form_data.append("from_full", sessionStorage.getItem("in_str"));
    form_data.append("to_full", sessionStorage.getItem("out_str"));
    form_data.append("success", quality);
    form_data.append("from_flags", read_flags);
    form_data.append("to_flags", write_flags);
    form_data.append("coordinates", coordinates);
    form_data.append("coordOption", coordOption);
    form_data.append("fileToUpload", file);
    form_data.append("upload_file", true);

    convertFile(form_data, download_fname, fname);
}

// Retrieves option flags from selected text
function extractFlags(flags_text) {
    var flags = "",
        regex = /:/g,
        match = "";

    while ((match = regex.exec(flags_text)) != null) {
        flags += flags_text[match.index - 1];
    }

    return flags;
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

// Retrieves read or write option flags associated with a file format
function getFlags (type, str) {
    var query = ``;

    if (type == "in") {
        // $$$$$$$$$ MAKE A FUNCTION FOR THIS? $$$$$$$$$
        const in_str_array = str.split(": "),
              in_ext = in_str_array[0],          // e.g. "ins"
              in_note = in_str_array[1];         // e.g. "ShelX"

        query = `SELECT DISTINCT Flag, Description, Further_info FROM OBFlags_in
                 WHERE ID IN (SELECT DISTINCT OBFlags_in_ID FROM OBFormat_to_Flags_in
                              WHERE Formats_ID=(SELECT ID FROM Formats WHERE Extension = '${in_ext}' AND Note = '${in_note}'))`;
    }
    else {
        const out_str_array = str.split(": "),
              out_ext = out_str_array[0],          // e.g. "ins"
              out_note = out_str_array[1];         // e.g. "ShelX"

        query = `SELECT DISTINCT Flag, Description, Further_info FROM OBFlags_out
                 WHERE ID IN (SELECT DISTINCT OBFlags_out_ID FROM OBFormat_to_Flags_out
                              WHERE Formats_ID=(SELECT ID FROM Formats WHERE Extension = '${out_ext}' AND Note = '${out_note}'))`;
    }

    try {
        queryDatabase(query, type, populateFlagBox);
        return true;
    }
    catch (e) {
        return false;
    }
}

// Populates a read or write option flag box
function populateFlagBox(response, type) {
    var el = $("#" + type + "Flags"),
        disp = $("#" + type + "FlagList"),
        flag = '',
        info = '<p>',
        count = 0;

    if (response.length != 0) {
        disp.css({display: "inline"});
        response += '$';

        for (var i = 1; i < response.length; i++) {
            if (response[i] == '$') {                 // End of all information about a particular option flag
                const original_length = info.length;
                info = info.replace(/.: N\/A/, '');   // Ensures 'N/A' not displayed

                if (info.length == original_length) {
                    info += '</p><p>';
                }

                count = -1;
            }
            else if (count == 2) {                    // End of what appears in the flag box for a particular option flag
                el.append(new Option(flag));
                info += response[i];                  // Start of further information about a particular option flag
                flag = '';
                count++;
            }
            else if (count == 3) {                    // Adds to further information
                info += response[i];                  
            }
            else if (response[i] == '£' && response[i - 1] != '$' && count == 0) {
                flag += ': ';                         // Break between flag and its description
                info += ': ';
                count++;
            }
            else if (response[i] == '£') {            // Acknowledges change from one piece of data to another
                count++;
            }
            else if (response[i] != '£') {            // Adds to what appears in the flag box
                flag += response[i];

                if (count == 0) {
                    info += response[i];              // Adds flag to start of further information
                }
            }
        }
    }

    el.append(new Option(""));
    $("#" + type + "FlagInfo").html(info);
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

// File upload is allowed only if its extension matches the 'from' format
function checkExtension(event) {
    const file_name = this.files[0].name;
    const file_name_array = file_name.split(".");
    const extension = file_name_array[1];

    if (extension != in_ext) {
        $("#uploadButton").css({"background-color": "#e5e1e6", "color": "gray"});
        $("#uploadButton").prop({disabled: true});
        alert("The file extension is not " + in_ext + ": please select another file or change the 'from' format on the 'Home' page.");
    }
    else {
        $("#uploadButton").css("background-color", "#011e41"); // TODO: COMBINE TWO LINES
        $("#uploadButton").css("color", "#e5e1e6");
        $("#uploadButton").prop({disabled: false});
    }
}

