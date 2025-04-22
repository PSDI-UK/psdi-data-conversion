/*
  format.js
  Version 1.0, 27th November 2024

  This is the JavaScript which makes the Format and Converter Selection gui work.
*/

import { disableDirtyForms, cleanDirtyForms, initDirtyForms, loadServiceMode, loadProductionMode } from "./common.js";
import {
    databaseLoaded, getInputFormats, getOutputFormats, getOutputFormatsForInputFormat,
    getInputFormatsForOutputFormat, getConverters, getConverterByName, getLevelChemInfo
} from "./data.js";

var fromList = new Array(),
    toList = new Array(),
    qualityCriteriaCount = 0,
    qualityMeasureSum = 0;

$(document).ready(function () {

    // Populate the available formats if the database successfully loaded, otherwise display an error message
    if (databaseLoaded()) {
        // Populates the "Convert from" selection list
        getInputFormats().then((formats) => {
            populateList(formats, "from");
        });

        // Populates the "Convert to" selection list
        getOutputFormats().then((formats) => {
            populateList(formats, "to");
        });
    } else {
        displayDatabaseLoadError();;
    }

    sessionStorage.setItem("token", token);
    sessionStorage.setItem("max_file_size", max_file_size);
    sessionStorage.setItem("max_file_size_ob", max_file_size_ob);
    sessionStorage.setItem("service_mode", service_mode);
    sessionStorage.setItem("production_mode", production_mode);
    sessionStorage.setItem("in_str", "");
    sessionStorage.setItem("out_str", "");
    sessionStorage.setItem("success", "");

    loadServiceMode();
    loadProductionMode();

    $("#fromList").click(populateConversionSuccess);
    $("#toList").click(populateConversionSuccess);
    $("#searchTo").keyup(filterOptions);
    $("#searchFrom").keyup(filterOptions);
    $("#yesButton").click(goToConversionPage);
    $("#success").click(showConverterDetails);
    $("#resetButton").click(resetAll);
    $("#showButton").click(showQualityDetails);

    initDirtyForms();
});

/**
 * Hides the format and converter selection and displays an error message that the database could not be loaded.
 */
function displayDatabaseLoadError() {
    $("#format-selection").css({ display: "none" });
    $("#database-error").css({ display: "block" });
}

/**
 * Gets the input and output extensions and their notes from what's currently selected in the search boxes.
 * @returns {Array<str>} Input extension, Input note, Output extension, Output note
 */
function getExtAndNotes() {
    const in_str = $("#searchFrom").val(), // e.g. "ins: ShelX"
        in_str_array = in_str.split(": "),
        in_ext = in_str_array[0],           // e.g. "ins"
        in_note = in_str_array[1];          // e.g. "ShelX"

    const out_str = $("#searchTo").val(),
        out_str_array = out_str.split(": "),
        out_ext = out_str_array[0],
        out_note = out_str_array[1];

    return [in_ext, in_note, out_ext, out_note];
}

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
        const [in_ext, in_note, out_ext, out_note] = getExtAndNotes();

        if (this.id == "fromList") {
            toList = [];
            $("#toList").children().remove();
            getOutputFormatsForInputFormat(in_ext, in_note).then(formats => populateList(formats, "to"));
        }
        else if (this.id == "toList") {
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

// Shows how the conversion quality was determined for the selected conversion/converter combination
function showQualityDetails(event) {
    const [in_ext, in_note, out_ext, out_note] = getExtAndNotes();
    const converter = sessionStorage.getItem("success").split(": ")[0]

    getLevelChemInfo(in_ext, in_note, out_ext, out_note).then(entries => displayLevelChemInfo(entries));
}

/**
 * Compiles details on quality issues for a conversion between two formats
 * @param {*} entries 
 * @param {bool} format 
 * @returns {Array<str>} A long description of quality issues, the percent rating of the conversion, a short description
 *                       of the quality
 */
function getQualityDetails(entries, format = false) {
    var composition_in = entries[0].composition,
        composition_out = entries[1].composition,
        connections_in = entries[0].connections,
        connections_out = entries[1].connections,
        two_dim_in = entries[0].two_dim,
        two_dim_out = entries[1].two_dim,
        three_dim_in = entries[0].three_dim,
        three_dim_out = entries[1].three_dim;

    qualityCriteriaCount = 0;
    qualityMeasureSum = 0;

    var quality = qualityDetail(composition_in, composition_out, 'Composition') +
        qualityDetail(connections_in, connections_out, 'Connections') +
        qualityDetail(two_dim_in, two_dim_out, '2D coordinates') +
        qualityDetail(three_dim_in, three_dim_out, '3D coordinates');
    let percent, qualityText;

    if (qualityCriteriaCount == 0) {
        quality = 'Conversion quality details not available for this converter/conversion combination.';
        percent = "N/A";
        qualityText = "N/A";
    }
    else {
        qualityText = 'very poor';
        percent = qualityMeasureSum * 20 / qualityCriteriaCount;

        percent = (Math.round(percent * 100) / 100);

        if (percent >= 80.0) {
            qualityText = 'very good';
        }
        else if (percent >= 60.0) {
            qualityText = 'good';
        }
        else if (percent >= 40.0) {
            qualityText = 'okay';
        }
        else if (percent >= 20.0) {
            qualityText = 'poor';
        }

        if (quality != '') {
            quality = "<strong>WARNING:</strong> Potential data loss or extrapolation in conversion due to a " +
                "mismatch in representation between input and output formats:\n<ul>" + quality + "</ul>";
        }
    }

    // Strip the HTML formatting from the output if it isn't desired
    if (!format) {
        quality = quality.replaceAll(/<.*?>/g, "");
    }

    return [quality, percent, qualityText];
}

// Assemble quality information based on level of chemical information and display it.
function displayLevelChemInfo(entries) {
    var [quality, percent, qualityText] = getQualityDetails(entries);

    if (percent != "N/A") {

        quality += '-----------------------------\n' +
            'Total score: ' + percent + '%\n' +
            'Conversion quality: ' + qualityText + '\n' +
            '-----------------------------\n';

    }

    alert(quality);
}

// Determine a statement of quality based on level of chemical information and populate the converter select box.
function getQuality(entries, rows) {
    var composition_in = entries[0].composition,
        composition_out = entries[1].composition,
        connections_in = entries[0].connections,
        connections_out = entries[1].connections,
        two_dim_in = entries[0].two_dim,
        two_dim_out = entries[1].two_dim,
        three_dim_in = entries[0].three_dim,
        three_dim_out = entries[1].three_dim,
        qualityText = ' Conversion quality is very poor.';

    qualityCriteriaCount = 0;
    qualityMeasureSum = 0;

    qualityDetail(composition_in, composition_out, 'Composition') +
        qualityDetail(connections_in, connections_out, 'Connections') +
        qualityDetail(two_dim_in, two_dim_out, '2D coordinates') +
        qualityDetail(three_dim_in, three_dim_out, '3D coordinates');

    if (qualityCriteriaCount == 0) {
        qualityText = ' Conversion quality not tested.';
    }
    else {
        var percent = qualityMeasureSum * 20 / qualityCriteriaCount;

        percent = (Math.round(percent * 100) / 100);

        if (percent >= 80.0) {
            qualityText = ' Conversion quality is very good.';
        }
        else if (percent >= 60.0) {
            qualityText = ' Conversion quality is good.';
        }
        else if (percent >= 40.0) {
            qualityText = ' Conversion quality is okay.';
        }
        else if (percent >= 20.0) {
            qualityText = ' Conversion quality is poor.';
        }
    }

    rows.sort(function (a, b) {
        return a.toLowerCase().localeCompare(b.toLowerCase());
    });

    for (var i = 0; i < rows.length; i++) {
        const support = rows[i].substring(0, 10) == "Open Babel" ||
            rows[i].substring(0, 6) == "Atomsk" ||
            rows[i].substring(0, 3) == "c2x" ? " supported on this site." : " not supported on this site." //" (supported)" : " (unsupported)";

        $("#success").append($('<option>', { text: rows[i] + support + qualityText }));
    }
}

// Returns 'true' if textbox text is exactly the same as a select box option
function isOption(str, boxId) {
    var isOption = false;

    $("#" + boxId + " > option").each(function () {
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
    $("#converter").css({ display: "none" });
    $("h3").css({ display: "none" });
}

// Show conversion offer
function showOffer() {

    const [in_ext, in_note, out_ext, out_note] = getExtAndNotes();

    const quest = " like to convert a file from '" + in_ext + "' to '" + out_ext + "' on this site";

    if ($("#name").html() == "Open Babel" || $("#name").html() == "Atomsk" || $("#name").html() == "c2x") {
        $("#info").html("");
        $("#visit").html("visit website");
        $("#question").html("Would you" + quest + " using " + $("#name").html() + "?");
        $("#yesButton").css({ display: "inline" });

        // Check for any warnings about this conversion
        getLevelChemInfo(in_ext, in_note, out_ext, out_note).then(function (entries) {
            const [quality, _percent, _qualityText] = getQualityDetails(entries, true);
            if (quality != "") {
                $("#formatWarning").html(quality);
                $("#formatWarning").css({ display: "block" });
            }
        });

    }
    else {
        $("#info").html("This converter is not supported on our website; however, you can find out how to use it at");
        $("#visit").html("this website.");
        $("#question").html("As an alternative, if you would" + quest + ", please select a supported converter above (if available).");
        $("#yesButton").css({ display: "none" });
    }

    $("#question").css({ display: "inline" });
    $("#offer").css({ display: "inline" });
}

// Hide conversion offer
function hideOffer() {
    $("#converter").css({ display: "none" });
    $("#question").css({ display: "none" });
    $("#offer").css({ display: "none" });
    $("#formatWarning").css({ display: "none" });
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
                $("#converter").css({ display: "block" });
                $("h3").css({ display: "block" });
            }

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
        });
    }
}

// Create content for conversion quality details based on level of chemical information
function qualityDetail(input, output, type) {

    let label;
    if (type.includes("2D")) {
        label = "2D atomic coordinates are"
    } else if (type.includes("3D")) {
        label = "3D atomic coordinates are"
    } else if (type.includes("Composition")) {
        label = "Atomic composition is"
    } else if (type.includes("Connections")) {
        label = "Atomic connections are"
    } else {
        label = "The '" + type + "' property is"
    }

    if (input == true && output == true) {
        qualityCriteriaCount += 1;
        qualityMeasureSum += 5;
        return '';
    }
    else if (input == true && output == false) {
        return '<li><strong>Potential data loss:</strong> ' + label +
            ' represented in the input format but not the output format.</li>\n';
    }
    else if (input == false && output == true) {
        // We penalize the quality if the output format contains information that the input doesn't, as this might
        // result in info being extrapolated
        qualityCriteriaCount += 1;
        return '<li><strong>Potential data extrapolation:</strong> ' + label +
            ' represented in the output format but not the input format.</li>\n';
    }
    else {
        return '';
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

    disableDirtyForms();

    var path = ``;

    if ($("#name").html() == "Open Babel") {
        path = `static/content/convert.htm`;
    }
    else if ($("#name").html() == "Atomsk") {
        path = `static/content/convertato.htm`;
    }
    else if ($("#name").html() == "c2x") {
        path = `static/content/convertc2x.htm`;
    }

    const a = $("<a>")
        .attr("href", path)
        .appendTo("body");

    a[0].click();
    a.remove();
}

// Populates a selection list
function populateList(entries, sel) {
    const [in_ext, in_note, out_ext, out_note] = getExtAndNotes();

    let rows;

    if ((sel === "from") || (sel === "to")) {

        rows = entries.map(entry => `${entry.extension}: ${entry.note}`);

    } else if (sel === "success") {

        rows = entries.map(entry => `${entry.name}:`);

        if (in_ext != "" && out_ext != "") {
            const quality = getLevelChemInfo(in_ext, in_note, out_ext, out_note).then(formats => getQuality(formats, rows));
        }
    }

    if (sel !== "success") {
        rows.sort(function (a, b) {
            return a.toLowerCase().localeCompare(b.toLowerCase());
        });

        for (var i = 0; i < rows.length; i++) {
            const support = rows[i].substring(0, 10) == "Open Babel" || rows[i].substring(0, 6) == "Atomsk" ? " (supported)" : " (unsupported)";

            if (sel == "from") {
                $("#fromList").append($('<option>', { value: rows[i], text: rows[i] }));
                fromList[i] = rows[i] + "\n";
            }
            else if (sel == "to") {
                $("#toList").append($('<option>', { value: rows[i], text: rows[i] }));
                toList[i] = rows[i] + "\n";
            }
        }
    }

    const in_str = in_ext + ": " + in_note;
    const out_str = out_ext + ": " + out_note;
    if (sel == "from" && !isOption(in_str, "toList")) {
        filter("searchFrom");
    }
    else if (sel == "to" && !isOption(out_str, "fromList")) {
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

    cleanDirtyForms();
}
