/*
  format.js
  Version 1.0, 16th May 2024

  This is the JavaScript which makes the Format and Converter Selection gui work.
*/

var fromList = new Array(),
    toList = new Array();

$(document).ready(function() {
    // Populates the "Convert from" and "Convert to" selection lists
    var query = `SELECT DISTINCT Extension, Note FROM Formats ORDER BY Extension, Note ASC`

    queryDatabase(query, "from", populateList);
    queryDatabase(query, "to", populateList);

    sessionStorage.setItem("token", token);
    sessionStorage.setItem("in_str", "");
    sessionStorage.setItem("out_str", "");
    sessionStorage.setItem("success", "");

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
        // Can do without an error message if the 'Conversion options' box remains empty;
        // however, consider a greyed-out message inside the box (using some of the commented out code below).
    }
}

// Retrieve selected text from the "Conversion success" textarea
// $$$$$$$$$$ Can delete this PROVIDED the mobile 'phone select box issue can be solved without it? BUT NEED TO DO IT ANOTHER WAY!! $$$$$$$$$$
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
    $("#radioButtons").css({display: "inline"}); // $$$$$$$$$$ TODO: 'radioButtons' is no longer an appropriate id $$$$$$$$$$
}

// Hide Open Babel conversion offer
function hideOffer() {
    $("#converter").css({display: "none"});
    $("#question").css({display: "none"});
    $("#radioButtons").css({display: "none"});
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

// Displays converter details given its name
function showConverterDetails(event) {
    var selectedText = getSelectedText(this);

    sessionStorage.setItem("success", selectedText);

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

    // Search textarea for "Open Babel"     $$$$$ textarea? $$$$$
    const text = el.value;

    el.selectionStart = 0;
    el.selectionEnd = 0;

    while (el.selectionStart < text.length) {
        const selectedText = getSelectedText(el),
              name = selectedText.split(": ")[0];

        showOffer();

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

// Empties the "Conversion success" textarea     $$$$$ textarea? $$$$$
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
}

// Resets the filtering, format list and converter list boxes
function resetAll() {
    $("#searchFrom").val("");
    $("#searchFrom").keyup();

    $("#searchTo").val("");
    $("#searchTo").keyup();
}
