/*
  accessibility.js
  Version 1.0, 7th June 2024

  This is the JavaScript which makes the Accessibility gui work.
*/

let r = document.querySelector(':root');

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
          darkColour = sessionStorage.getItem("darkColour"),
          darkColourOpt = sessionStorage.getItem("darkColourOpt"),
          lightColour = sessionStorage.getItem("lightColour"),
          lightColourOpt = sessionStorage.getItem("lightColourOpt"),
          back = sessionStorage.getItem("back"),
          backOpt = sessionStorage.getItem("backOpt");

    if (font != null) {
        $(".normalText, .middle, #resetButton, #applyButton").css({
            fontFamily: font,
            fontSize: size,
            fontWeight: weight,
            letterSpacing: letter
        });

        if (darkColour !== "default ") {
            r.style.setProperty('--psdi-dark-text-color-body', darkColour);
            r.style.setProperty('--psdi-dark-text-color-heading', darkColour);
        }

        if (lightColour !== "default ") {
            r.style.setProperty('--psdi-dark-text-color-body', lightColour);
            r.style.setProperty('--psdi-dark-text-color-heading', lightColour);
        }

        $(".middle").css({lineHeight: line});
        $("h1, h2").css({letterSpacing: letter});
        $("form, select, #upper").css({background: back});

        $("#font").val(fontOpt).change();
        $("#size").val(sizeOpt).change();
        $("#weight").val(weightOpt).change();
        $("#letter").val(letterOpt).change();
        $("#line").val(lineOpt).change();
        $("#darkColour").val(darkColourOpt).change();
        $("#lightColour").val(lightColourOpt).change();
        $("#background").val(backOpt).change();
    }

    $("#font").change(changeFont);
    $("#size").change(changeFontSize);
    $("#weight").change(changeFontWeight);
    $("#letter").change(changeLetterSpacing);
    $("#line").change(changeLineSpacing);
    $("#dark-colour").change(changeFontColourDark);
    $("#light-colour").change(changeFontColourLight);
    $("#background").change(changeBackground);
    $("#resetButton").click(resetSelections);
    $("#applyButton").click(applySettings);
});

// Changes the font for accessibility purposes and ensures that the correct default line spacing is applied.
function changeFont(event) {
    const font = $("#font").find(":selected").text(),
          line = $("#line").find(":selected").text();

    var text = $(".normalText, .middle, #resetButton, #applyButton");

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
    var text = $(".normalText, .middle, #resetButton, #applyButton, h1, h2");

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
    var text = $(".normalText, .middle, #resetButton, #applyButton");

    switch (size) {
        case '14':
            text.css({fontSize: '14px'});
            break;

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
            text.css({fontSize: '1rem'});
            break;
    }
}

// Changes the font weight for accessibility purposes.
function changeFontWeight(event) {
    const weight = $("#weight").find(":selected").text();
    var text = $(".normalText, .middle, #resetButton, #applyButton");

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

function changeFontColourDark(event) {
    return changeFontColour(event, "dark");
}

function changeFontColourLight(event) {
    return changeFontColour(event, "light");
}

function changeFontColour(event, lightOrDark="dark") {
    const colour = $("#"+lightOrDark+"-colour").find(":selected").text();
    var text = $(".normalText");
    let new_colour;

    switch (colour) {
        case 'Black':
            new_colour = 'black';
            break;

        case 'White':
            new_colour = 'white';
            break;

        case 'Red':
            new_colour = 'red';
            break;

        case 'Orange':
            new_colour = 'orange';
            break;

        case 'Green':
            new_colour = 'green';
            break;

        case 'Purple':
            new_colour = 'purple';
            break;

        case 'Brown':
            new_colour = 'brown';
            break;

        default:
            new_colour = 'default';
            break;
    }

    if (new_colour==='default') {
        if (lightOrDark=="dark") {
            r.style.setProperty('--psdi-'+lightOrDark+'-text-color-body', '#000000');
            r.style.setProperty('--psdi-'+lightOrDark+'-text-color-heading', '#041E42');
        } else {
            r.style.setProperty('--psdi-'+lightOrDark+'-text-color-body', '#ffffff');
            r.style.setProperty('--psdi-'+lightOrDark+'-text-color-heading', '#ffffff');
        }
    } else {
        r.style.setProperty('--psdi-'+lightOrDark+'-text-color-body', new_colour);
        r.style.setProperty('--psdi-'+lightOrDark+'-text-color-heading', new_colour);
    }
}

// Changes the background colour for accessibility purposes.
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

// Reverts all select boxes to 'Default'
function resetSelections(event) {
    $("#font, #size, #weight, #letter, #line, #colour, #background").val('Default').change();
}

// Applies accessibility settings to the entire website.
function applySettings(event) {
    sessionStorage.setItem("font", $(".normalText").css('fontFamily'));
    sessionStorage.setItem("size", $(".normalText").css('fontSize'));
    sessionStorage.setItem("weight", $(".normalText").css('fontWeight'));
    sessionStorage.setItem("letter", $(".normalText").css('letterSpacing'));
    sessionStorage.setItem("line", $(".normalText").css('lineHeight'));
    sessionStorage.setItem("darkColour", r.style.getPropertyValue('--psdi-dark-text-color-body'));
    sessionStorage.setItem("lightColour", r.style.getPropertyValue('--psdi-light-text-color-body'));
    sessionStorage.setItem("back", $("form").css('background'));

    sessionStorage.setItem("fontOpt", $("#font").find(":selected").text());
    sessionStorage.setItem("sizeOpt", $("#size").find(":selected").text());
    sessionStorage.setItem("weightOpt", $("#weight").find(":selected").text());
    sessionStorage.setItem("letterOpt", $("#letter").find(":selected").text());
    sessionStorage.setItem("lineOpt", $("#line").find(":selected").text());
    sessionStorage.setItem("darkColourOpt", $("#dark-colour").find(":selected").text());
    sessionStorage.setItem("lightColourOpt", $("#light-colour").find(":selected").text());
    sessionStorage.setItem("backOpt", $("#background").find(":selected").text());

    alert("The settings have been applied to the entire website.");
}

