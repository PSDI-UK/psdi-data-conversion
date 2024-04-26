/*
  format.js
  Version 1.0, 26th April 2024

  This is the JavaScript which makes the Format and Converter Selection gui work.
*/

var fromList = new Array(),
    toList = new Array(),
    last_select = "";

$(document).ready(function() {
    // Populates the "Convert from" and "Convert to" selection lists
    var query = `SELECT DISTINCT Extension, Note FROM Formats ORDER BY Extension, Note ASC`

    queryDatabase(query, "from", populateList);
    queryDatabase(query, "to", populateList);

    sessionStorage.setItem("token", token);
    sessionStorage.setItem("in_str", "");
    sessionStorage.setItem("out_str", "");
//    sessionStorage.setItem("message", "");
  //  sessionStorage.setItem("message1", "");

    $("#fromList").click(populateConversionSuccess);
    $("#toList").click(populateConversionSuccess);
    $("#searchTo").keyup(filterOptions);
    $("#searchFrom").keyup(filterOptions);
    $("#yesButton").click(goToConversionPage);
    $("#success").click(showConverterDetails);
    $("#resetButton").click(resetAll);
});

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

    sessionStorage.setItem("in_str", from_text);
    sessionStorage.setItem("out_str", to_text);

    this.selectionStart = -1;
    this.selectionEnd = -1;
    this.blur();
    emptySuccess();
    hideConverterDetails();
    hideOffer();

    try {
        const in_str = $("#searchFrom").val(), // e.g. "ins: ShelX"
              in_str_array = in_str.split(": "),
              in_ext = in_str_array[0],           // e.g. "ins"
              in_note = in_str_array[1];          // e.g. "ShelX"

        const out_str = $("#searchTo").val(),
              out_str_array = out_str.split(": "),
              out_ext = out_str_array[0],
              out_note = out_str_array[1];

        const query = `SELECT C.Name, C_to.Degree_of_success FROM Converters C, Converts_to C_to
                       WHERE C_to.in_ID=(SELECT ID FROM Formats WHERE Extension = '${in_ext}' AND Note = '${in_note}')
                       AND C_to.out_ID=(SELECT ID FROM Formats WHERE Extension = '${out_ext}' AND Note = '${out_note}')
                       AND C.ID=C_to.Converters_ID ORDER BY C.Name ASC`

        queryDatabase(query, "success", populateList);
    }
    catch (e) {
        const ID_a = getFormat($("#searchFrom").val()),
              ID_b = getFormat($("#searchTo").val());

        if (ID_a.toString() != ID_b.toString() && ID_a != "" && ID_b != "") {
                conversionSuccessEmpty();
        }
    }
}

// Retrieve selected text from the "Conversion success" textarea
function getSelectedText(el) {
    const text = el.value,
          before = text.substring(0, el.selectionStart),
          after = text.substring(el.selectionEnd, text.length);

    el.selectionStart = before.lastIndexOf("\n") >= 0 ? before.lastIndexOf("\n") + 1 : 0;
    el.selectionEnd = after.indexOf("\n") >= 0 ? el.selectionEnd + after.indexOf("\n") : text.length;

    return el.value.substring(el.selectionStart, el.selectionEnd);
}

// Hides converter details
function hideConverterDetails() {
    $("#converter").css({display: "none"});
    $("h3").css({display: "none"});
}

// Prompts the user for feedback if "File format not found" is selected
//function formatNotFound(element) {
  //  emptySuccess();
    //hideConverterDetails();
//    showTextInput();
  //  last_select = element.id;

//    sessionStorage.setItem("message", "");
  //  sessionStorage.setItem("message1", "Missing file format? If so, please enter the format, the format(s) " +
    //                       (element.attr('id') == "searchFrom" ? "to" : "from") + " which it should be converted and a reason.");
//}

// Prompts the user for feedback if the "Conversion success" selection list is empty
function conversionSuccessEmpty() {
//    sessionStorage.setItem("message", "Missing conversion? If you believe that this conversion should be supported, please enter a reason.");
  //  sessionStorage.setItem("message1", "The displayed 'from' and 'to' formats will be automatically added to your message.");

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
    const from_format = getFormat($("#searchFrom").val()),
          to_format = getFormat($("#searchTo").val()),
          quest = "ould you like to convert a file from '" + from_format + "' to '" + to_format + "' on this site using Open Babel?";

    if ($("#name").html() == "Open Babel") {
        $("#info").html("");
        $("#visit").html("visit website");
        $("#question").html("W" + quest);
    }
    else {
        $("#info").html("This converter is not currently supported on our website; however, you can find out how to use it at");
        $("#visit").html("this website.");
        $("#question").html("As an alternative, w" + quest);
    }

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

//    if (!($("#searchFrom").val() == "File format not found" ||
  //        $("#searchTo").val() == "File format not found")) {
    $("#success").prop({disabled: false});
    //}
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

// Displays converter details given its name
function showConverterDetails(event) {
    var selectedText = getSelectedText(this);

    const str_array = selectedText.split(": ", 1),
          conv_name = str_array[0];                                     // e.g. "Open Babel"

    const query = `SELECT * FROM Converters WHERE Name='${conv_name}'`;

    queryDatabase(query, this, displayConverterDetails);
}

// Displays converter details and offers an Open Babel conversion if available 
function displayConverterDetails(response, el) {
    var count = 0,
        str = '';

    response += '£';

    try {
        for (var i = 1; i < response.length; i++) {
            if (response[i] == '£') {
                if (count == 1) {
                    $("#name").html(str);
                }
                else if (count == 2) {
                    $("#description").html(str);
                }
                else if (count == 3) {
                    $("#url").html(str);
                    break;
                }

                count += 1;
                str = '';
            }
            else {
                str += response[i];
            }
        }

        var visit = $("#visit");
        visit.attr("href", str);

        $("#converter").css({display: "block"});
        $("h3").css({display: "block"});
    }
    catch (e) {
        // No need for error message if no options shown
    }

    // Search textarea for "Open Babel"
    const text = el.value;

    el.selectionStart = 0;
    el.selectionEnd = 0;

    while (el.selectionStart < text.length) {
        const selectedText = getSelectedText(el),
              name = selectedText.split(": ")[0];

        if (name == "Open Babel") {
            showOffer();
            break;
        }
        else {
            $("#info").html("This converter is not currently supported on our website; however, you can find out how to use it at");
        }

        el.selectionEnd += 1;
        el.selectionStart = el.selectionEnd;
    }

    el.selectionStart = -1;
    el.selectionEnd = -1;
    el.blur();
}

// Only options having user filter input as a substring (case insensitive) are included in the selection list
function filterOptions(event) {
    const str = this.value.toLowerCase();
    var box, list,
        count = 0,
        text = "";

    if (this.id == "searchFrom") {
        box = $("#fromList");
        list = fromList;
    }
    else {
        box = $("#toList");
        list = toList;
    }

    box.children().remove();

    for (var i = 0; i < list.length; i++) {
        if (list[i].toLowerCase().includes(str)) {
            box.append($('<option>', { text: list[i] }));
            count += 1;
        }
    }

    if (this.id == "searchFrom") {
        $("#fromLabel").html("Convert from (" + count + "):");
    }
    else {
        $("#toLabel").html("Convert to (" + count + "):");
    }

    $("#success").prop({disabled: true});
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

// Stores chosen formats and switches to the Conversion page
function goToConversionPage(event) {
    const a = $("<a>")
          .attr("href", `static/content/convert.htm`)
          .appendTo("body");

    a[0].click();
    a.remove();
}

// Shows file selection and conversion elements
function showFileConversionDetails(event) {
    $("#convertFile").css({display: "block"});

    var in_found = getFlags("in"),
        out_found = getFlags("out");

    if (in_found || out_found) {
        $("#flags").css({display: "grid"});
    }
}

// Retrieves read or write option flags associated with a file format
function getFlags (type) {
    var query = ``;

    if (type == "in") {
        // $$$$$$$$$ MAKE A FUNCTION FOR THIS? $$$$$$$$$
        const in_str = $("#searchFrom").val(),   // e.g. "ins: ShelX"
              in_str_array = in_str.split(": "),
              in_ext = in_str_array[0],          // e.g. "ins"
              in_note = in_str_array[1];         // e.g. "ShelX"

        query = `SELECT DISTINCT Flag, Description FROM OBFlags_in
                 WHERE ID IN (SELECT DISTINCT OBFlags_in_ID FROM OBFormat_to_Flags_in
                              WHERE Formats_ID=(SELECT ID FROM Formats WHERE Extension = '${in_ext}' AND Note = '${in_note}'))`;
    }
    else {
        const out_str = $("#searchTo").val(),      // e.g. "ins: ShelX"
              out_str_array = out_str.split(": "),
              out_ext = out_str_array[0],          // e.g. "ins"
              out_note = out_str_array[1];         // e.g. "ShelX"

        query = `SELECT DISTINCT Flag, Description FROM OBFlags_out
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

// $$$$$$$$$$ COMMENT NEEDED HERE $$$$$$$$$$
function populateFlagBox(response, type) {
    var el = $("#" + type + "Flags"),
        flag = '';

    for (var i = 1; i < response.length; i++) {
        if (response[i] == '£' && response[i - 1] != '$') {
            flag += ': ';
        }
        else if (response[i] == '$') {
            el.append(new Option(flag));
            flag = '';
        }
        else if (i == response.length - 1) { // $$$$$$$$ PUT A '$' AT THE END OF 'response' instead? $$$$$$$$
            flag += response[i];
            el.append(new Option(flag));
        }
        else if (response[i] != '£') {
            flag += response[i];
        }
    }

    el.append(new Option(""));
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
function populateList(response, sel) {
    var el = $("#" + sel),
        successText = "",
        format = '',
        rows = [];

    for (var i = 1; i < response.length; i++) {
        if (response[i] == '£' && response[i - 1] != '$') {
            format += ': ';
        }
        else if (response[i] == '$') {
            rows.push(format);
            format = '';
        }
        else if (i == response.length - 1) {
            format += response[i];
            rows.push(format);
        }
        else if (response[i] != '£') {
            format += response[i];
        }
    }

    rows.sort(function(a, b) {
        return a.toLowerCase().localeCompare(b.toLowerCase());
    });

    $("#success").prop({disabled: true});

    for (var i = 0; i < rows.length; i++) {
        const support = rows[i].substring(0, 10) == "Open Babel" ? " (supported)" : " (unsupported)";

        if ( sel == "success") {
            if (rows.length > 0) {
                $("#success").prop({disabled: false});
            }

            $("#success").append($('<option>', { text: "" + rows[i] + support }));
        }

        if (sel == "from") {
            $("#fromList").append($('<option>', { text: rows[i] }));
            fromList[i] = rows[i] + "\n";
        }
        else if (sel == "to") {
            $("#toList").append($('<option>', { text: rows[i] }));
            toList[i] = rows[i] + "\n";
        }
    }

    if (sel != "success") {
        $("#fromLabel").html("Convert from (" + fromList.length + "):");
        $("#toLabel").html("Convert to (" + toList.length + "):");
    }
    else {
        const ID_a = getFormat($("#searchFrom").val()),
              ID_b = getFormat($("#searchTo").val());

        if (ID_a.toString() != ID_b.toString() && ID_a != "" && ID_b != "") {
            conversionSuccessEmpty();
        }
    }
}

// Resets the filtering, format list and converter list boxes
function resetAll() {
    $("#searchFrom").val("");
    $("#searchFrom").keyup();

    $("#searchTo").val("");
    $("#searchTo").keyup();
}
