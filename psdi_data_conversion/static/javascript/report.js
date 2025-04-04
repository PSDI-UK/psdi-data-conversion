/*
  report.js
  Version 1.0, 26th June 2024

  This is the JavaScript which makes the report.htm gui work.
*/

import { disableDirtyForms, initDirtyForms } from "./common.js";
import { getAllFormats, getConverters } from "./data.js";

var token = "",
    fromList = new Array(),
    toList = new Array(),
    formatList = new Array();

const pageURL = new URL(window.location.toLocaleString());

// If coming to this page from a local version of the app, warn if the user clicks a link in the header, which will
// take them to another page on the public app
const originUrlParam = pageURL.searchParams.get('origin');
const localOrigin = (originUrlParam != null && originUrlParam.toLowerCase() == "local") ? true : false;
let originAlertDisplayed = false;

if (localOrigin) {
    $("#psdi-header").click(function (e) {
        if (originAlertDisplayed) return;
        alert("WARNING: This is the online, public version of the app, which you arrived at from the local " +
            "app. Press the back button on your browser to return to the local version, or click \"OK\" on this " +
            "alert to continue to the public app.");
        originAlertDisplayed = true;
    });
}

$(document).ready(function () {

    token = sessionStorage.getItem("token");

    $("#success").css({ display: "none" });

    // Populates the "Convert from" and "Convert to" selection lists
    getAllFormats().then((allFormats) => {
        populateList(allFormats, "from");
        populateList(allFormats, "to");
        populateList(allFormats, "format");
    });

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

    initDirtyForms();
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

        getConverters(in_ext, in_note, out_ext, out_note).then((converters) => {
            populateList(converters, "success");
        });
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
    $("#converter").css({ display: "none" });
    $("h3").css({ display: "none" });
}

// Submits user input
function submitUserInput() {

    const from = $("#searchFrom").val(),
        to = $("#searchTo").val();

    var reason = $("#in").val(),
        missing = $("#missingFormat").val()

    if (reason.length > 9 && reason.length < 501) {
        if ($("#reason").val() == "format") {
            if (missing.length > 1 && missing.length < 101) {
                submitFeedback({
                    type: "missingFormat",
                    missing: missing,
                    reason: reason
                });
            }
            else {
                alert("Please enter the missing format (2 to 100 characters).");
            }
        }
        else {
            if (from != "" && to != "") {
                if ($("#success option").length == 0) {
                    submitFeedback({
                        type: "missingConversion",
                        from: from,
                        to: to,
                        reason: reason
                    });
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
function hideOffer() { }

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

// Submit feedback
function submitFeedback(data) {
    disableDirtyForms();

    $.post(`/feedback/`, {
        'token': token,
        'data': JSON.stringify(data)
    })
        .done(() => {
            alert("Report received!");
        })
        .fail(function (e) {
            alert("Reporting failed. Please provide feedback by clicking on 'Contact' in the navigation bar.");

            // For debugging
            console.error("Error submitting feedback", e.status, e.responseText);
        });
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

    $("#success").prop({ disabled: true });
    emptySuccess();
    hideConverterDetails();
    hideOffer();
}

// Empties the "Conversion success" textarea
function emptySuccess() {
    $("#success").html("");
}

// Populates a selection list
function populateList(entries, sel) {

    let rows = [];

    if ((sel === "from") || (sel === "to") || (sel === "format")) {

        rows = entries.map(entry => `${entry.extension}: ${entry.note}`);

    } else if (sel === "success") {

        rows = entries.map(entry => `${entry.name}: ${entry.degree_of_success}`);
    }

    rows.sort(function (a, b) {
        return a.toLowerCase().localeCompare(b.toLowerCase());
    });

    $("#success").prop({ disabled: true });

    for (var i = 0; i < rows.length; i++) {
        const support = rows[i].substring(0, 10) == "Open Babel" ? " (supported)" : " (unsupported)";

        if (sel == "success") {
            if (rows.length > 0) {
                $("#success").prop({ disabled: false });
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
        $("#in_out_formats").css({ display: "block" });
        $("#message").css({ display: "inline" });

        $("#userInput").css({ display: "block" });

        $("#formats").css({ display: "none" });
        $("#missing").css({ display: "none" });

        $("#message").html("Explain why the conversion is required and provide a link to appropriate documentation if possible [max 500 characters].");
        $("#message1").html("The displayed 'from' and 'to' formats will be automatically submitted with your message.");
    }
    else if (selectedText == "format") {
        $("#formats").css({ display: "block" });
        $("#missing").css({ display: "block" });

        $("#userInput").css({ display: "block" });

        $("#in_out_formats").css({ display: "none" });
        $("#message").css({ display: "none" });

        $("#message1").html("Enter details of the file conversions expected for this format and provide a link to appropriate documentation if possible [max 500 characters].");
    }
    else {
        $("#in_out_formats").css({ display: "none" });
        $("#message").css({ display: "none" });

        $("#formats").css({ display: "none" });
        $("#missing").css({ display: "none" });

        $("#userInput").css({ display: "none" });
    }
}

