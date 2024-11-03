/*
  format.js
  Version 1.0, 26th June 2024

  This is the JavaScript which makes the Format and Converter Selection gui work.
*/

import { loadAccessibilitySettings } from './accessibility.js';
import { getInputFormats, getOutputFormats, getOutputFormatsForInputFormat,
    getInputFormatsForOutputFormat, getConverters, getConverterByName } from "./data.js";

var fromList = new Array(),
    toList = new Array();

$(document).ready(function() {

    loadAccessibilitySettings();

    // Populates the "Convert from" selection list
    getInputFormats().then((formats) => {
        populateList(formats, "from");
    });

    // Populates the "Convert to" selection list
    getOutputFormats().then((formats) => {
        populateList(formats, "to");
    });

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

// Remove the loading cover when everything is loaded
$(window).on('load', function() {
    $("#cover").hide();
});

// Selects a file format; populates the "Conversion success" selection list given input and output IDs;
// filters the output format list when an input format is selected, and vice versa (formats not convertible
// to/from the selected format are removed); and removes converter details and text input (if showing)
function populateConversionSuccess(event) {
    const selectedText = getSelectedText(this);
    let filterQuery = ``;

    if (this.id == "fromList") {
        $("#searchFrom").val(selectedText);
    }
    else {
        $("#searchTo").val(selectedText);
    }

    const from_text = $("#searchFrom").val(),
          to_text = $("#searchTo").val();

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

        if (this.id == "fromList") { // && !isOption(out_str, "toList")) {
            toList = [];
            $("#toList").children().remove();
            getOutputFormatsForInputFormat(in_ext, in_note).then(formats => populateList(formats, "to"));
        }
        else if (this.id == "toList") { // && !isOption(in_str, "fromList")) {
            fromList = [];
            $("#fromList").children().remove();
            getInputFormatsForOutputFormat(out_ext, out_note).then(formats => populateList(formats, "from"));
        }

        getConverters(in_ext, in_note, out_ext, out_note).then((converters) => {
            populateList(converters, "success");
        });
    }
    catch (e) {
        // Can do without an error message if the 'Conversion options' box remains empty;
        // however, consider a greyed-out message inside the box (using some of the commented out code below).
    }
}

// Returns 'true' if textbox text is exactly the same as a select box option
function isOption(str, boxId) {
    var isOption = false;

    $("#" + boxId + " > option").each(function() {
        if (str == this.text) {
            isOption = true;
        }
    });

    return isOption;
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
        $("#info").html("This converter is not supported on our website; however, you can find out how to use it at");
        $("#visit").html("this website.");
        $("#question").html("As an alternative, w" + quest);
    }

    $("#question").css({display: "inline"});
    $("#offer").css({display: "inline"});
}

// Hide Open Babel conversion offer
function hideOffer() {
    $("#converter").css({display: "none"});
    $("#question").css({display: "none"});
    $("#offer").css({display: "none"});
}

// Displays converter details given its name
function showConverterDetails(event) {
    var selectedText = getSelectedText(this);

    if (selectedText != "") {
        sessionStorage.setItem("success", selectedText);

        const str_array = selectedText.split(": ", 1),
              conv_name = str_array[0];                                     // e.g. "Open Babel"

        getConverterByName(conv_name).then((converter) => {
            if (converter) {
                $("#name").html(converter.name);
                $("#description").html(converter.description);
                $("#url").html(converter.url);
                $("#visit").attr("href", converter.url);
                $("#converter").css({display: "block"});
                $("h3").css({display: "block"});
            }
        });

        const el = this;

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
}

// Only options having user filter input as a substring (case insensitive) are included in the selection list $$$$$$$$$$ REVISE $$$$$$$$$$
function filterOptions(event) {
    const str = event.target.value.toLowerCase();
    var box, list,
        count = 0,
        text = "";

    if (event.target.id == "searchFrom") {
        toList = [];
        $("#toList").children().remove();
        getOutputFormats().then(formats => populateList(formats, "to"));
        box = $("#fromList");
        list = fromList;
    }
    else {
        fromList = [];
        $("#fromList").children().remove();
        getInputFormats().then(formats => populateList(formats, "from"));
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

    if (event.target.id == "searchFrom") {
        $("#fromLabel").html("Select format to convert from (" + count + "):");
    }
    else {
        $("#toLabel").html("Select format to convert to (" + count + "):");
    }

    emptySuccess();
    hideConverterDetails();
    hideOffer();
}

// Only options having user filter input as a substring (case insensitive) are included in the slection list
function filter(id) {
    try {
        const str = $("#" + id).val().toLowerCase();
        var box, list,
            count = 0,
            text = "";

        if (id == "searchFrom") {
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

        if (id == "searchFrom") {
            $("#fromLabel").html("Select format to convert from (" + count + "):");
        }
        else {
            $("#toLabel").html("Select format to convert to (" + count + "):");
        }

        emptySuccess();
        hideConverterDetails();
        hideOffer();
    }
    catch (e) {
        // No need for an error message here. No need to filter if text box is empty.
    }
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
function populateList(entries, sel) {
    const in_str = $("#searchFrom").val(), // e.g. "ins: ShelX"
          out_str = $("#searchTo").val();

    let rows;

    if ((sel === "from") || (sel === "to")) {

        rows = entries.map(entry => `${entry.extension}: ${entry.note}`);

    } else if (sel === "success") {

        rows = entries.map(entry => `${entry.name}: ${entry.degree_of_success}`);
    }

    rows.sort(function(a, b) {
        return a.toLowerCase().localeCompare(b.toLowerCase());
    });

    for (var i = 0; i < rows.length; i++) {
        const support = rows[i].substring(0, 10) == "Open Babel" ? " (supported)" : " (unsupported)";

        if ( sel == "success") {
            $("#success").append($('<option>', { text: rows[i] + support }));
        }

        if (sel == "from") {
            $("#fromList").append($('<option>', { value: rows[i], text: rows[i] }));
            fromList[i] = rows[i] + "\n";
        }
        else if (sel == "to") {
            $("#toList").append($('<option>', { value: rows[i], text: rows[i] }));
            toList[i] = rows[i] + "\n";
        }
    }

    if (sel == "from" && !isOption(out_str, "toList")) {
        filter("searchFrom");
    }
    else if (sel == "to" && !isOption(in_str, "fromList")) {
        filter("searchTo");
    }

    if (sel == "from") {
        $("#fromLabel").html("Select format to convert from (" + $("#fromList").children('option').length + "):");
        $("#fromList option[value='" + in_str + "']").prop('selected', 'selected');
        $("#fromList").hide().show();
    }
    else if (sel == "to") {
        $("#toLabel").html("Select format to convert to (" + $("#toList").children('option').length + "):");
        $("#toList option[value='" + out_str + "']").prop('selected', 'selected');
        $("#toList").hide().show();
    }
}

// Resets the filtering, format list and converter list boxes
function resetAll() {
    $("#fromList").children().remove();
    $("#toList").children().remove();

    $("#searchFrom").val("");
    $("#searchTo").val("");

    // Populates the "Convert from" selection list
    getInputFormats().then(formats => populateList(formats, "from"));

    // Populates the "Convert to" selection list
    getOutputFormats().then(formats => populateList(formats, "to"));
}
