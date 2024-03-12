/*
  format.js
  Version 1.0, 11th March 2024

  This is the JavaScript which makes the gui work.
*/

import format_db_byte_array from './byte.js';

var fromList = new Array();
var toList = new Array();
var last_select = "";
var db;

$(document).ready(function() {
    initSqlJs({
        locateFile: filename => 'static/javascript/node_modules/sql.js/dist/sql-wasm.wasm'}).then(SQL => {
        db = new SQL.Database(format_db_byte_array);

        // Populates the "Convert from" selection list
        var query = `SELECT DISTINCT Form.Extension, Form.Note FROM Formats Form, Converts_to Conv
                     WHERE Form.ID=Conv.in_ID UNION SELECT DISTINCT Extension, Note
                     FROM OBFormats WHERE Input="true" ORDER BY Extension, Note ASC`
        populateList(query, "from");

        // Populates the "Convert to" selection list
        var query = `SELECT DISTINCT Form.Extension, Form.Note FROM Formats Form, Converts_to Conv
                     WHERE Form.ID=Conv.out_ID UNION SELECT DISTINCT Extension, Note
                     FROM OBFormats WHERE Output="true" ORDER BY Extension, Note ASC`
        populateList(query, "to");
    });

    $("#fromList").click(populateConversionSuccess);
    $("#toList").click(populateConversionSuccess);
    $("#searchTo").keyup(filterOptions);
    $("#searchFrom").keyup(filterOptions);
    $("#enter").click(submitUserInput);
    $("#cancel").click(hideTextInput);
    $("#noteButton").click(toggleNotes);
    $("#radioYes").click(showFileConversionDetails);
    $("#radioNo").click(hideFileConversionDetails);
    $("#fileToUpload").change(checkExtension);
    $("#uploadButton").click(submitFile);
    $("#success").click(showConverterDetails);
});

// Database must be closed to avoid memory leak
window.onbeforeunload = function() {
    db.close();
};

// Finds ID from table Formats given Extension and Note
function findFormatID(sel) {
    const str = $("#" + sel).val();    // e.g. "ins: ShelX"
    const str_array = str.split(": ");
    const ext = str_array[0];          // e.g. "ins"
    const note = str_array[1];         // e.g. "ShelX"

    try {
        var result = db.exec(`SELECT ID FROM Formats
                              WHERE Extension="${ext}" AND Note="${note}"`);
    }
    catch (e) {
        console.log(e);
        return "";
    }

    return result[0].values;
}

// Finds ID from table OBFormats given Extension
function findOBFormatID(sel) {
    const str = $("#" + sel).val();   // e.g. "ins: ShelX"
    const str_array = str.split(": ");
    const ext = str_array[0];         // e.g. "ins"

    try {
        var result = db.exec(`SELECT ID FROM OBFormats WHERE Extension="${ext}"`);
    }
    catch (e) {
        console.log(e);
        return "";
    }

    return result[0].values;
}

// Checks for the existence of a selected file extension in table OBFormats
function findExtension(sel) {
    alert("Top of findExtension");
    const ext = getFormat($("#" + sel).val()); // e.g. "ins: ShelX" --> "ins"
    alert("ext = " + ext);
    try {
        var result = db.exec(`SELECT Extension FROM OBFormats WHERE Extension="${ext}"`);
        if (result.toString() == "") {
            return "";
        }
        else {
            return result[0].values;
        }
    }
    catch (e) {
        console.log(e);
        return "";
    }
}

// Selects a file format; populates the "Conversion success" selection list given input and output IDs;
// and removes converter details and text input (if showing)
function populateConversionSuccess(event) {
    const selectedText = getSelectedText(this);

    if (this.id == "fromList") {
        $("#searchFrom").val(selectedText);
    }
    else {
        $("#searchTo").val(selectedText);
    }

    const from_text = $("#searchFrom").val();
    const to_text = $("#searchTo").val();

    if (!(from_text == "File format not found" || to_text == "File format not found")) {
        hideConverterDetails();
        hideTextInput();
    }
    else if (from_text == "File format not found") {
        formatNotFound($("#searchFrom"));
    }
    else {
        formatNotFound($("#searchTo"));
    }

    this.selectionStart = -1;
    this.selectionEnd = -1;
    this.blur();
    emptySuccess();
    hideOffer();

    if (!(to_text == "-- select --" || to_text == "File format not found" ||
          from_text == "-- select --" || from_text == "File format not found")) {
        try {
            const ID_a = findFormatID("searchFrom");
            const ID_b = findFormatID("searchTo");

            const query = `SELECT C.Name, C_to.Degree_of_success FROM Converters C, Converts_to C_to
                           WHERE C_to.in_ID=${ID_a} AND C_to.out_ID=${ID_b} AND C.ID=C_to.Converters_ID ORDER BY C.Name ASC`

            populateList(query, "success");
        }
        catch (e) {
            const ID_a = getFormat($("#searchFrom").val());
            const ID_b = getFormat($("#searchTo").val());

            if (ID_a.toString() != ID_b.toString() && ID_a != "" && ID_b != "") {
                $("#success").html("-- select --");
                conversionSuccessEmpty();

                var query = `SELECT DISTINCT Extension, Note
                             FROM OBFormats WHERE Input="true"`

                var inputs = db.exec(query)[0].values;

                query = `SELECT DISTINCT Extension, Note
                         FROM OBFormats WHERE Output="true"`

                var outputs = db.exec(query)[0].values;

                for (var i = 0; i < inputs.length; i++) {
                    if (from_text == inputs[i].join(": ")) {
                        for (var j = 0; j < outputs.length; j++) {
                            if (to_text == outputs[j].join(": ")) {
                                showOffer();
                                break;
                            }
                        }
                    }
                }
            }
            else {
                console.log(e + "\nNot necessarily a problem - unsupported conversions are not in the database.");
            }
        }
    }
}

// Retrieve selected text from the "Conversion success" textarea
function getSelectedText(el) {
    const text = el.value;
    const before = text.substring(0, el.selectionStart);
    const after = text.substring(el.selectionEnd, text.length);

    el.selectionStart = before.lastIndexOf("\n") >= 0 ? before.lastIndexOf("\n") + 1 : 0;
    el.selectionEnd = after.indexOf("\n") >= 0 ? el.selectionEnd + after.indexOf("\n") : text.length;

    return el.value.substring(el.selectionStart, el.selectionEnd);
}

// Hides converter details
function hideConverterDetails() {
    $("#converter").css({display: "none"});
    $("h2").css({display: "none"});
}

// Prompts the user for feedback if "File format not found" is selected
function formatNotFound(element) {
    emptySuccess();
    hideConverterDetails();
    showTextInput();
    last_select = element.id;

    $("#message").html("");
    $("#message1").html("Missing file format? If so, please enter the format, the format(s) " +
                        (element.attr('id') == "searchFrom" ? "to" : "from") + " which it should be converted and a reason.");

    $("#enter").css({"background-color": "#993366",
                     "color": "white",
                     "font-family": "Lato",
                     "font-size": "14px"});
    $("#cancel").css({"background-color": "#993366",
                      "color": "white",
                      "font-family": "Lato",
                      "font-size": "14px"});
}

// Prompts the user for feedback if the "Conversion success" selection list is empty
function conversionSuccessEmpty() {
    $("#message").html("Missing conversion? If you believe that this conversion should be supported, please enter a reason.");
    $("#message1").html("The displayed 'from' and 'to' formats will be automatically added to your message.");

    $("#success").disabled = true;
    showTextInput();
    last_select = "success";
}

// Shows the text input element and associated buttons and prompt
function showTextInput() {
    $("#userInput").css({display: "block"});
}

// Submits user input
function submitUserInput() {
    const from = $("#searchFrom").val();
    const to = $("#searchTo").val();
    const reason = $("#in").val();
    var input = '';

    if (reason != "") {
        if (last_select == "success") {
            input = 'From ' + from + ' to ' + to + '   ' + reason;
        }
        else {
            input = reason;
        }

        writeLog(input);
        hideTextInput();
    }
}

// Show Open Babel conversion offer
function showOffer() {
    const from_format = getFormat($("#searchFrom").val());
    const to_format = getFormat($("#searchTo").val());

    $("#question").html("Would you like to convert a file from " + from_format + " to " + to_format + " using Open Babel?");
    $("#question").css({display: "inline"});
    $("#radioButtons").css({display: "inline"});
    $("#radioNo").prop("checked", true);
}

// Hide Open Babel conversion offer
function hideOffer() {
    $("#converter").css({display: "none"});
    $("#question").css({display: "none"});
    $("#radioButtons").css({display: "none"});
    $("#convertFile").css({display: "none"});
    $("#inFlags").empty();
    $("#outFlags").empty();
    $("#flags").css({display: "none"});
}

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

// Uploads a user-supplied file
function submitFile() {
    const file = $("#fileToUpload")[0].files[0],
          extension = file.name.split(".")[1];
    
    const from = $("#searchFrom").val(),    // e.g. "ins: ShelX"
          from_format = from.split(": ")[0];

    if (extension != from_format) {
        alert("The file extension is not " + from_format + ": please select another file or change the 'from' format.");
        return;
    }

    const to = $("#searchTo").val(),
          to_format = to.split(": ")[0];

    const read_flags_text = $("#inFlags").find(":selected").text(),
          read_flags = extractFlags(read_flags_text);

    const write_flags_text = $("#outFlags").find(":selected").text(),
          write_flags = extractFlags(write_flags_text);

    const download_fname = file.name.split(".")[0] + "." + to_format;

    var form_data = new FormData();

    form_data.append("token", token);
    form_data.append("from", from_format);
    form_data.append("to", to_format);
    form_data.append("from_flags", read_flags);
    form_data.append("to_flags", write_flags);
    form_data.append("fileToUpload", file);
    form_data.append("upload_file", true);

    convertFile(form_data, download_fname);
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
function convertFile(form_data, download_fname) {
    var jqXHR = $.ajax({
            url: `/convert/`,
            type: "POST",
            data: form_data,
            processData: false,
            contentType: false,
            success: function() {
                const a = $("<a>")
                      .attr("href", `static/downloads/${download_fname}`)
                      .attr("download", download_fname)
                      .appendTo("body");
                a[0].click();
                a.remove();
                },
            error: function(data) {
                alert("ajax error, FormData: " + data);
                }
            })
            .fail(function(e) {
                // For debugging
                console.log("Error converting file");
                console.log(e.status);
                console.log(e.responseText);
            })
}

// Hides the text input element and associated buttons and prompt
function hideTextInput() {
    $("#in").val("");
    $("#userInput").css({display: "none"});

    if (!($("#searchFrom").val() == "File format not found" ||
          $("#searchTo").val() == "File format not found")) {
        $("#success").prop({disabled: false});
    }
}

// Shows notes if hidden; hides notes if shown. Button text changed as appropriate.
function toggleNotes(event) {
    if (this.value == "Show notes") {
        this.value = "Hide notes";
        $("#notes").css({visibility: "visible"});
    }
    else {
        this.value = "Show notes";
        $("#notes").css({visibility: "hidden"});
    }
}

// Displays converter details given its name and offers an Open Babel conversion if available
function showConverterDetails(event) {
    var selectedText = getSelectedText(this);

    const text = this.value,
          str_array = selectedText.split(": ", 1),
          conv_name = str_array[0],                                                    // e.g. "Open Babel"
          result = db.exec(`SELECT * FROM Converters WHERE Name="${conv_name}"`);

    try {
        $("#name").html("" + result[0].values[0][1]);
        $("#description").html("" + result[0].values[0][2]);
        $("#url").html("" + result[0].values[0][3]);

        var visit = $("#visit");
        visit.attr("href", "" + result[0].values[0][3]);

        $("#converter").css({display: "block"});
        $("h2").css({display: "block"});
    }
    catch (e) {
        // No need for error message if no options shown
    }

    // Search textarea for "Open Babel"
    this.selectionStart = 0;
    this.selectionEnd = 0;

    while (this.selectionStart < text.length) {
        selectedText = getSelectedText(this);
        const name = selectedText.split(": ")[0];

        if (name == "Open Babel") {
            showOffer();
            break;
        }

        this.selectionEnd += 1;
        this.selectionStart = this.selectionEnd;
    }

    this.selectionStart = -1;
    this.selectionEnd = -1;
    this.blur();
}

// Only options having user filter input as a substring (case insensitive) are included in the selection list
function filterOptions(event) {
    const str = this.value.toLowerCase();
    var box;
    var list;
    var text = "-- select --\n";

    if (this.id == "searchFrom") {
        box = document.getElementById("fromList");
        list = fromList;
    }
    else {
        box = document.getElementById("toList");
        list = toList;
    }

    for (var i = 1; i < list.length - 1; i++) {
        if (list[i].toLowerCase().includes(str)) {
            text += list[i];
        }
    }

    text += "File format not found";
    box.value = text;
    emptySuccess();
    hideConverterDetails();
    hideOffer();
}

// Empties the "Conversion success" textarea
function emptySuccess() {
    $("#success").html("");
}

// Retrieves a file format from a string (e.g. "ins: ShelX") from a selection list
function getFormat(str) {
    const str_array = str.split(": ");
    return str_array[0];               // e.g. "ins"
}

// Shows file selection and conversion elements
function showFileConversionDetails(event) {
    $("#convertFile").css({display: "block"});
    $("#uploadButton").css({"background-color": "#993366",
                            "color": "white",
                            "font-family": "Lato",
                            "font-size": "14px"});

    var in_found = getFlags("in"),
        out_found = getFlags("out");

    if (in_found || out_found) {
        $("#flags").css({display: "grid"});
    }
}

// Retrieves read or write option flags associated with a file format
function getFlags (type) {
    var id = "",
        query = ``,
        el = $("#" + type + "Flags");

    if (type == "in") {
        id = findOBFormatID("searchFrom");

        query = `SELECT DISTINCT Flag, Description FROM OBFlags_in
                 WHERE ID IN (SELECT DISTINCT OBFlags_in_ID
                              FROM OBFormat_to_Flags_in WHERE OBFormats_ID=${id})`;
    }
    else {
        id = findOBFormatID("searchTo");

        query = `SELECT DISTINCT Flag, Description FROM OBFlags_out
                 WHERE ID IN (SELECT DISTINCT OBFlags_out_ID
                              FROM OBFormat_to_Flags_out WHERE OBFormats_ID=${id})`;
    }

    try {
        var result = db.exec(query)[0].values;
        el.append(new Option("-- select if required --"));

        for (var i = 0; i < result.length; i++) {
            var row = result[i].join(": ");
            el.append(new Option(row));
        }

        el.append(new Option(""));
        return true;
    }
    catch (e) {
        el.append(new Option("-- no option flags available --"));
        el.append(new Option(""));
        return false;
    }
}

// Hides file selection and conversion elements
function hideFileConversionDetails(event) {
    $("#convertFile").css({display: "none"});
    $("#inFlags").empty();
    $("#outFlags").empty();
    $("#flags").css({display: "none"});
}

// File upload is allowed only if its extension matches the 'from' format
function checkExtension(event) {
    const file_name = this.files[0].name;
    const file_name_array = file_name.split(".");
    const extension = file_name_array[1];
    const from_format = getFormat($("#searchFrom").val());

    if (extension != from_format) {
        $("#uploadButton").prop({disabled: true});
        alert("The file extension is not " + from_format + ": please select another file or change the 'from' format.");
    }
    else {
        $("#uploadButton").prop({disabled: false});
    }
}

// Populates a selection list
function populateList(query, sel) {
    var el = $("#" + sel),
        successText = "-- select --",
        text = "-- select --\n",
        rows = [];

    const result = db.exec(query)[0].values;

    while (result.length > 0) {
        rows.push(result.pop().join(": "));
    };

    rows.sort(function(a, b) {
        return a.toLowerCase().localeCompare(b.toLowerCase());
    });

    fromList[0] = "-- select --\n";
    toList[0] = "-- select --\n";

    for (var i = 0; i < rows.length; i++) {

        if ( sel == "success") {
            successText += "\n" + rows[i];
        }

        text += rows[i] + "\n";

        // This is for creating permanent arrays of options
        if (sel == "from") {  
            fromList[i + 1] = rows[i] + "\n";
        }
        else if (sel == "to") {
            toList[i + 1] = rows[i] + "\n";
        }
    }

    text += "File format not found";

    if (sel == "from") {
        $("#fromList").html(text);
        fromList[fromList.length] = "File format not found";
    }
    else if (sel == "to") {
        $("#toList").html(text);
        toList[toList.length] = "File format not found";
    }
    else {
        el.html(successText);
    }
}

