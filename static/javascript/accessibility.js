/*
  accessibility.js
  Version 1.0, 31st May 2024

  This is the JavaScript which makes the Accessibility gui work.
*/

var fromList = new Array(),
    toList = new Array();

$(document).ready(function() {
    const font = sessionStorage.getItem("font"),
          fontOpt = sessionStorage.getItem("fontOpt"),
          size = sessionStorage.getItem("size"),
          sizeOpt = sessionStorage.getItem("sizeOpt"),
          weight = sessionStorage.getItem("weight"),
          weightOpt = sessionStorage.getItem("weightOpt"),
          letter = sessionStorage.getItem("letter"),
          letterOpt = sessionStorage.getItem("letterOpt"),
          line = sessionStorage.getItem("line"),
          lineOpt = sessionStorage.getItem("lineOpt"),
          colour = sessionStorage.getItem("colour"),
          colourOpt = sessionStorage.getItem("colourOpt"),
          back = sessionStorage.getItem("back"),
          backOpt = sessionStorage.getItem("backOpt");

    if (font != null) {
        $(".normalText, .middle, #applyButton").css({
            fontFamily: font,
            fontSize: size,
            fontWeight: weight,
            letterSpacing: letter
        });

        $(".normalText").css({
            lineHeight: line,
            color: colour
        });

        $(".middle").css({lineHeight: line});
        $("h1, h2").css({letterSpacing: letter});
        $("h1").css({color: colour});
        $("form, select, #upper").css({background: back});

        $("#font").val(fontOpt).change();
        $("#size").val(sizeOpt).change();
        $("#weight").val(weightOpt).change();
        $("#letter").val(letterOpt).change();
        $("#line").val(lineOpt).change();
        $("#colour").val(colourOpt).change();
        $("#background").val(backOpt).change();
    }

    $("#font").change(changeFont);
    $("#size").change(changeFontSize);
    $("#weight").change(changeFontWeight);
    $("#letter").change(changeLetterSpacing);
    $("#line").change(changeLineSpacing);
    $("#colour").change(changeFontColour);
    $("#background").change(changeBackground);
    $("#applyButton").click(applySettings);
/*    $("#fromList").click(populateConversionSuccess);
    $("#toList").click(populateConversionSuccess);
    $("#searchTo").keyup(filterOptions);
    $("#searchFrom").keyup(filterOptions);
    $("#yesButton").click(goToConversionPage);
    $("#success").click(showConverterDetails);
    $("#resetButton").click(resetAll);         */
});

// Changes the font for accessibility purposes and ensures that the correct default line spacing is applied.
function changeFont(event) {
    const font = $("#font").find(":selected").text(),
          line = $("#line").find(":selected").text();

    var text = $(".normalText, .middle, #applyButton");

    switch (font) {
        case 'Arial':
            if (line == 'Default') {
                text.css({lineHeight: 1.145});
            }

            text.css({fontFamily: 'Arial, sans-serif'});
            break;

        case 'Comic Sans':
            if (line == 'Default') {
                text.css({lineHeight: 1.4});
            }

            text.css({fontFamily: 'Comic Sans MS, Comic Sans, sans-serif'});
            break;

        case 'Lexend':
            if (line == 'Default') {
                text.css({lineHeight: 1.3});
            }

            text.css({fontFamily: 'Lexend, sans-serif'});
            break;

        case 'Open Sans':
            if (line == 'Default') {
                text.css({lineHeight: 1.4});
            }

            text.css({fontFamily: 'Open Sans, sans-serif'});
            break;

        case 'Tahoma':
            if (line == 'Default') {
                text.css({lineHeight: 1.25});
            }

            text.css({fontFamily: 'Tahoma, sans-serif'});
            break;

        case 'Trebuchet':
            if (line == 'Default') {
                text.css({lineHeight: 1.2});
            }

            text.css({fontFamily: 'Trebuchet MS, Trebuchet, sans-serif'});
            break;

        case 'Verdana':
            if (line == 'Default') {
                text.css({lineHeight: 1.25});
            }

            text.css({fontFamily: 'Verdana, sans-serif'});
            break;

        default:
            if (line == 'Default') {
                text.css({lineHeight: 1.218});
            }

            text.css({fontFamily: 'Lato, sans-serif'});
            break;
    }
}

// Changes the letter spacing for accessibility purposes.
function changeLetterSpacing(event) {
    const space = $("#letter").find(":selected").text();
    var text = $(".normalText, .middle, #applyButton, h1, h2");

    switch (space) {
        case '0.5':
            text.css({letterSpacing: '0.5px'});
            break;

        case '1.0':
            text.css({letterSpacing: '1px'});
            break;

        case '1.5':
            text.css({letterSpacing: '1.5px'});
            break;

        case '2.0':
            text.css({letterSpacing: '2px'});
            break;

        case '2.5':
            text.css({letterSpacing: '2.5px'});
            break;

        case '3.0':
            text.css({letterSpacing: '3px'});
            break;

        case '3.5':
            text.css({letterSpacing: '3.5px'});
            break;

        default:
            text.css({letterSpacing: '0px'});
            break;
    }
}

// Changes the line spacing for accessibility purposes.
function changeLineSpacing(event) {
    const space = $("#line").find(":selected").text();
    var text = $(".normalText, .middle");

    switch (space) {
        case '1.1':
            text.css({lineHeight: 1.1});
            break;

        case '1.2':
            text.css({lineHeight: 1.2});
            break;

        case '1.3':
            text.css({lineHeight: 1.3});
            break;

        case '1.4':
            text.css({lineHeight: 1.4});
            break;

        case '1.5':
            text.css({lineHeight: 1.5});
            break;

        case '1.6':
            text.css({lineHeight: 1.6});
            break;

        case '1.7':
            text.css({lineHeight: 1.7});
            break;

        // Ensures that the correct default line spacing is applied to the current font.
        default:
            const font = $("#font").find(":selected").text();

            switch (font) {
                case 'Arial':
                    text.css({lineHeight: 1.145});
                    break;

                case 'Comic Sans':
                case 'Open Sans':
                    text.css({lineHeight: 1.4});
                    break;

                case 'Lexend':
                    text.css({lineHeight: 1.3});
                    break;

                case 'Tahoma':
                case 'Verdana':
                    text.css({lineHeight: 1.25});
                    break;

                case 'Trebuchet':
                    text.css({lineHeight: 1.2});
                    break;

                default:
                    text.css({lineHeight: 1.218});
                    break;
            }

            break;
    }
}

// Changes the font size for accessibility purposes.
function changeFontSize(event) {
    const size = $("#size").find(":selected").text();
    var text = $(".normalText, .middle, #applyButton");

    switch (size) {
        case '15':
            text.css({fontSize: '15px'});
            break;

        case '16':
            text.css({fontSize: '16px'});
            break;

        case '17':
            text.css({fontSize: '17px'});
            break;

        case '18':
            text.css({fontSize: '18px'});
            break;

        case '19':
            text.css({fontSize: '19px'});
            break;

        case '20':
            text.css({fontSize: '20px'});
            break;

        case '21':
            text.css({fontSize: '21px'});
            break;

        default:
            text.css({fontSize: '14px'});
            break;
    }
}

// Changes the font weight for accessibility purposes.
function changeFontWeight(event) {
    const weight = $("#weight").find(":selected").text();
    var text = $(".normalText, .middle, #applyButton");

    switch (weight) {
        case 'Bold':
            text.css({fontWeight: 'bold'});
            break;

        default:
            text.css({fontWeight: 'normal'});
            break;
    }
}

// Changes the font colour for accessibility purposes.
function changeFontColour(event) {
    const colour = $("#colour").find(":selected").text();
    var text = $(".normalText, h1");

    switch (colour) {
        case 'Black':
            text.css({color: 'black'});
            break;

        case 'Red':
            text.css({color: 'red'});
            break;

        case 'Orange':
            text.css({color: 'orange'});
            break;

        case 'Green':
            text.css({color: 'green'});
            break;

        case 'Purple':
            text.css({color: 'purple'});
            break;

        case 'Brown':
            text.css({color: 'brown'});
            break;

        default:
            text.css({color: '#011e41'});
            break;
    }
}

// Changes the background colour for accessibility purposes. Font colour is black when background colour is not white.
function changeBackground(event) {
    const colour = $("#background").find(":selected").text();
    var text = $(".normalText");

    switch (colour) {
        case 'Mustard':
            $("form, select, #upper").css({background: '#eddd6e'});
            break;

        case 'Peach':
            $("form, select, #upper").css({background: '#edd1b0'});
            break;

        case 'Lemon':
            $("form, select, #upper").css({background: '#f8fd89'});
            break;

        default:
            $("form, select, #upper").css({background: 'white'});
            break;
    }
}

// Applies accessibility settings to the entire website.
function applySettings(event) {
    sessionStorage.setItem("font", $(".normalText").css('fontFamily'));
    sessionStorage.setItem("size", $(".normalText").css('fontSize'));
    sessionStorage.setItem("weight", $(".normalText").css('fontWeight'));
    sessionStorage.setItem("letter", $(".normalText").css('letterSpacing'));
    sessionStorage.setItem("line", $(".normalText").css('lineHeight'));
    sessionStorage.setItem("colour", $(".normalText").css('color'));
    sessionStorage.setItem("back", $("form").css('background'));

    sessionStorage.setItem("fontOpt", $("#font").find(":selected").text());
    sessionStorage.setItem("sizeOpt", $("#size").find(":selected").text());
    sessionStorage.setItem("weightOpt", $("#weight").find(":selected").text());
    sessionStorage.setItem("letterOpt", $("#letter").find(":selected").text());
    sessionStorage.setItem("lineOpt", $("#line").find(":selected").text());
    sessionStorage.setItem("colourOpt", $("#colour").find(":selected").text());
    sessionStorage.setItem("backOpt", $("#background").find(":selected").text());
}

// Selects a file format; populates the "Conversion success" selection list given input and output IDs;
// and removes converter details and text input (if showing)
/*function populateConversionSuccess(event) {
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
}     */

// Retrieve selected text from the "Conversion success" textarea
// $$$$$$$$$$ Can delete this PROVIDED the mobile 'phone select box issue can be solved without it? BUT NEED TO DO IT ANOTHER WAY!! $$$$$$$$$$
/*function getSelectedText(el) {
    const text = el.value,
          before = text.substring(0, el.selectionStart),
          after = text.substring(el.selectionEnd, text.length);

    el.selectionStart = before.lastIndexOf("\n") >= 0 ? before.lastIndexOf("\n") + 1 : 0;
    el.selectionEnd = after.indexOf("\n") >= 0 ? el.selectionEnd + after.indexOf("\n") : text.length;

    return el.value.substring(el.selectionStart, el.selectionEnd);
} */

// Hides converter details
/*function hideConverterDetails() {
    $("#converter").css({display: "none"});
    $("h3").css({display: "none"});
} */

// Show Open Babel conversion offer
/*function showOffer() {
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
} */

// Hide Open Babel conversion offer
/*function hideOffer() {
    $("#converter").css({display: "none"});
    $("#question").css({display: "none"});
    $("#radioButtons").css({display: "none"});
} */

// $$$$$$$$$$ write separate function for debugging $$$$$$$$$$$

// Queries the PostgreSQL database
/*function queryDatabase(query, sel, callback) {
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
} */

// Displays converter details given its name
/*function showConverterDetails(event) {
    var selectedText = getSelectedText(this);

    sessionStorage.setItem("success", selectedText);

    const str_array = selectedText.split(": ", 1),
          conv_name = str_array[0];                                     // e.g. "Open Babel"

    const query = `SELECT * FROM Converters WHERE Name='${conv_name}'`;

    queryDatabase(query, this, displayConverterDetails);
} */

// Displays converter details and offers an Open Babel conversion if available 
/*function displayConverterDetails(response, el) {
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
} */

// Only options having user filter input as a substring (case insensitive) are included in the selection list
/*function filterOptions(event) {
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
} */

// Empties the "Conversion success" textarea     $$$$$ textarea? $$$$$
/*function emptySuccess() {
    $("#success").html("");
} */

// Retrieves a file format from a string (e.g. "ins: ShelX") from a selection list
/*function getFormat(str) {
    const str_array = str.split(": ");
    return str_array[0];               // e.g. "ins"
} */

// Stores chosen formats and switches to the Conversion page
/*function goToConversionPage(event) {
    const a = $("<a>")
          .attr("href", `static/content/convert.htm`)
          .appendTo("body");

    a[0].click();
    a.remove();
} */

// Populates a selection list
/*function populateList(response, sel) {
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
} */

// Resets the filtering, format list and converter list boxes
/*function resetAll() {
    $("#searchFrom").val("");
    $("#searchFrom").keyup();

    $("#searchTo").val("");
    $("#searchTo").keyup();
} */
