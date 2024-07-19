/*
  report.js
  Version 1.0, 26th June 2024

  This is the JavaScript which makes the report.htm gui work.
*/

var token = "",
    fromList = new Array(),
    toList = new Array(),
    formatList = new Array();

$(document).ready(function() {
    token = sessionStorage.getItem("token");

    $("#success").css({display: "none"});

    // Populates the "Convert from" and "Convert to" selection lists
    var query = `SELECT DISTINCT Extension, Note FROM Formats ORDER BY Extension, Note ASC`

    queryDatabase(query, "from", populateList);
    queryDatabase(query, "to", populateList);
    queryDatabase(query, "format", populateList);

    const font = sessionStorage.getItem("font"),
          size = sessionStorage.getItem("size"),
          weight = sessionStorage.getItem("weight"),
          letter = sessionStorage.getItem("letter"),
          line = sessionStorage.getItem("line"),
          colour = sessionStorage.getItem("colour"),
          back = sessionStorage.getItem("back");

    if (font != null) {
        $(".normalText, .middle, #resetButton, #resetButton2, #reportButton").css({
            fontFamily: font,
            fontSize: size,
            fontWeight: weight,
            letterSpacing: letter
        });

        $(".normalText, .middle").css({lineHeight: line});
        $(".normalText, h1").css({color: colour});
        $("h1, h2").css({letterSpacing: letter});
        $("form, select, #upper, #missingFormat, #searchFrom, #searchTo, #searchFormats, #in").css({background: back});
    }

//    $("#searchFrom").val(sessionStorage.getItem("in_str")); // $$$$$ CHECK if can remove these items from entire website $$$$$
  //  $("#searchTo").val(sessionStorage.getItem("out_str"));

    $("#reason").change(display);
    $("#fromList").click(populateConversionSuccess);
    $("#toList").click(populateConversionSuccess);
    $("#formatList").click(populateConversionSuccess);
    $("#searchTo").keyup(filterOptions);
    $("#searchFrom").keyup(filterOptions);
    $("#searchFormats").keyup(filterOptions);
    $("#resetButton").click(resetAll);
    $("#resetButton2").click(resetAll);
    $("#reportButton").click(submitUserInput);
});

// Included in this file for convenience. When the 'Report' button is clicked, a user's missing conversion report
// is only sent if the undisplayed conversion success box is empty (i.e., the conversion really is missing)
function populateConversionSuccess(event) {
    const selectedText = getSelectedText(this);

    if (this.id == "fromList") {
        $("#searchFrom").val(selectedText);
    }
    else if (this.id == "toList") {
        $("#searchTo").val(selectedText);
    }
    else {
        $("#searchFormats").val(selectedText);
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

//        const ID_a = getFormat($("#searchFrom").val()),
  //            ID_b = getFormat($("#searchTo").val());

    //    if (ID_a.toString() != ID_b.toString() && ID_a != "" && ID_b != "") {
      //          conversionSuccessEmpty();
        //}
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
    $("h3").css({display: "none"});
}

// Submits user input
function submitUserInput() {
    const from = $("#searchFrom").val(),
          to = $("#searchTo").val();

    var reason = $("#in").val(),
        missing = $("#missingFormat").val(),
        date = new Date(),
        now = "";

    now += date.toDateString();
    now += " " + date.toTimeString();

    reason = reason.replaceAll("'", "`");
    missing = missing.replaceAll("'", "`");

    if (reason.length > 9 && reason.length < 501) {
        if ($("#reason").val() == "format") {
            if (missing.length > 1 && missing.length < 101) {
                const query = `INSERT INTO Missing_Format (ID, Date, Format, Reason, Reviewed)
                               VALUES ((SELECT COALESCE(MAX(ID), 0) FROM Missing_Format) + 1, '${now}', '${missing}', '${reason}', 'no')`;

                updateDatabase(query);
            }
            else {
                alert("Please enter the missing format (2 to 100 characters).");
            }
        }
        else {
            if (from != "" && to != "") {
                if ($("#success option").length == 0) {
                    const query = `INSERT INTO Missing_Conversion (ID, Date, Convert_from, Convert_to, Reason, Reviewed)
                               VALUES ((SELECT COALESCE(MAX(ID), 0) FROM Missing_Conversion) + 1, '${now}', '${from}', '${to}', '${reason}', 'no')`;

                    updateDatabase(query);
                }
                else {
                    alert("At least one converter is capable of carrying out this conversion, therefore your report has not been sent. If you wish to send feedback about this conversion, please click on 'Contact' in the navigation bar.");
                }
            }
            else if (to != "") {
                alert("Please select 'from' format.");
            }
            else if (from != "") {
                alert("Please select 'to' format.");
            }
            else {
                alert("Please select 'to' and 'from' formats.");
            }
        }
    }
    else {
        alert("Please enter a reason, etc. (10 to 500 characters).");
    }
}

// Hide Open Babel conversion offer (not required to do anything on this page)
function hideOffer() {}

// Writes user input to a server-side file
// $$$$$$$$$$ Retain for now in case logging to file is required for some other purpose $$$$$$$$$$
//function writeLog(message) {
  //  var jqXHR = $.get(`/data/`, {
    //        'token': token,
      //      'data': message
        //})
//        .done(response => {
  //          alert("Report received!");
    //    })
      //  .fail(function(e) {
        //    alert("Reporting failed. Please provide feedback by clicking on 'Contact' in the navigation bar.");

          //  // For debugging
            //console.log("Error writing to log");
//            console.log(e.status);
  //          console.log(e.responseText);
    //    })
//}

// $$$$$$$$$$ TODO: write separate function for debugging $$$$$$$$$$$

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

// Updates the PostgreSQL database
function updateDatabase(query) {
    var jqXHR = $.post(`/query/`, {
        'token': token,
        'data': query
    })
    .done(response => {
        alert("Report received!");
    })
    .fail(function(e) {
        alert("Reporting failed. Please provide feedback by clicking on 'Contact' in the navigation bar.");

        // For debugging
        console.log("Error writing to log");
        console.log(e.status);
        console.log(e.responseText);
    })
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
    else if (this.id == "searchTo") {
        box = $("#toList");
        list = toList;
    }
    else {
        box = $("#formatList");
        list = formatList;
    }

    box.children().remove();

    for (var i = 0; i < list.length; i++) {
        if (list[i].toLowerCase().includes(str)) {
            box.append($('<option>', { text: list[i] }));
            count += 1;
        }
    }

    if (this.id == "searchFrom") {
        $("#fromLabel").html("Select format to convert from (" + count + "):");
    }
    else if (this.id == "searchTo") {
        $("#toLabel").html("Select format to convert to (" + count + "):");
    }
    else {
        $("#formatLabel").html("Check that the format is not present in the list. If it is, consider reporting a missing conversion. (" + count + ")");
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
        else if (sel == "format") {
            $("#formatList").append($('<option>', { text: rows[i] }));
            formatList[i] = rows[i] + "\n";
        }
    }

    if (sel != "success") {
        $("#fromLabel").html("Select format to convert from (" + fromList.length + "):");
        $("#toLabel").html("Select format to convert to (" + toList.length + "):");
        $("#formatLabel").html("Check that the format is not present in the list. If it is, consider reporting a missing conversion. (" + toList.length + ")");
    }
}

// Resets the filtering and format list boxes
function resetAll() {
    $("#searchFrom").val("");
    $("#searchFrom").keyup();

    $("#searchTo").val("");
    $("#searchTo").keyup();

    $("#searchFormats").val("");
    $("#searchFormats").keyup();
}

// Displays format or conversion related content as appropriate
function display(event) {
    const selectedText = getSelectedText(this);

    $("#in").val("");
    $("#missingFormat").val("");

    if (selectedText == "conversion") {
        $("#in_out_formats").css({display: "block"});
        $("#formats").css({display: "none"});
        $("#missing").css({display: "none"});
        $("#message").css({display: "inline"});
        $("#message").html("Explain why the conversion is required and provide a link to appropriate documentation if possible [max 500 characters].");
        $("#message1").html("The displayed 'from' and 'to' formats will be automatically submitted with your message.");
        $("#userInput").css({display: "block"});
    }
    else if (selectedText == "format") {
        $("#formats").css({display: "block"});
        $("#in_out_formats").css({display: "none"});
        $("#missing").css({display: "block"});
        $("#message").css({display: "none"});
        $("#message1").html("Enter details of the file conversions expected for this format and provide a link to appropriate documentation if possible [max 500 characters].");
        $("#userInput").css({display: "block"});
    }
    else {
        $("#formats").css({display: "none"});
        $("#missing").css({display: "none"});
        $("#userInput").css({display: "none"});
    }
}

