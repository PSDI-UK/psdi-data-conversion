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

        r.style.setProperty('--ifm-font-family-base', font);
        r.style.setProperty('--ifm-heading-font-family', font);

        r.style.setProperty('--ifm-font-size-base', size);

        r.style.setProperty('--ifm-font-size-base', size);

        r.style.setProperty('--ifm-font-weight-base', weight);

        r.style.setProperty('--psdi-letter-spacing-base', letter);

        r.style.setProperty('--psdi-dark-text-color-body', darkColour);
        r.style.setProperty('--psdi-dark-text-color-heading', darkColour);
        r.style.setProperty('--psdi-light-text-color-body', lightColour);
        r.style.setProperty('--psdi-light-text-color-heading', lightColour);

        r.style.setProperty('--ifm-line-height-base', line);

        r.style.setProperty('--ifm-background-color', back);

        $("#font").val(fontOpt).change();
        $("#size").val(sizeOpt).change();
        $("#weight").val(weightOpt).change();
        $("#letter").val(letterOpt).change();
        $("#line").val(lineOpt).change();
        $("#dark-colour").val(darkColourOpt).change();
        console.log("\""+lightColourOpt+"\"");
        $("#light-colour").val(lightColourOpt).change();
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
            r.style.setProperty('--ifm-font-family-base', 'Arial, sans-serif');
            r.style.setProperty('--ifm-heading-font-family', 'Arial, sans-serif');
            break;

        case 'Comic Sans':
            r.style.setProperty('--ifm-font-family-base', 'Comic Sans MS, Comic Sans, sans-serif');
            r.style.setProperty('--ifm-heading-font-family', 'Comic Sans MS, Comic Sans, sans-serif');
            break;

        case 'Lexend':
            r.style.setProperty('--ifm-font-family-base', 'Lexend, sans-serif');
            r.style.setProperty('--ifm-heading-font-family', 'Lexend, sans-serif');
            break;

        case 'Open Sans':
            r.style.setProperty('--ifm-font-family-base', 'Open Sans, sans-serif');
            r.style.setProperty('--ifm-heading-font-family', 'Open Sans, sans-serif');
            break;

        case 'Tahoma':
            r.style.setProperty('--ifm-font-family-base', 'Tahoma, sans-serif');
            r.style.setProperty('--ifm-heading-font-family', 'Tahoma, sans-serif');
            break;

        case 'Trebuchet':
            r.style.setProperty('--ifm-font-family-base', 'Trebuchet MS, Trebuchet, sans-serif');
            r.style.setProperty('--ifm-heading-font-family', 'Trebuchet MS, Trebuchet, sans-serif');
            break;

        case 'Verdana':
            r.style.setProperty('--ifm-font-family-base', 'Verdana, sans-serif');
            r.style.setProperty('--ifm-heading-font-family', 'Verdana, sans-serif');
            text.css({fontFamily: 'Verdana, sans-serif'});
            break;

        default:
            r.style.setProperty('--ifm-font-family-base', 'Lato, Helvetica, Arial, Lucida, sans-serif');
            r.style.setProperty('--ifm-heading-font-family', 'Oswald, Helvetica, Arial, Lucida, sans-serif');
            break;
    }
}

// Changes the letter spacing for accessibility purposes.
function changeLetterSpacing(event) {
    const space = $("#letter").find(":selected").text();

    if (space == "Default") {
        r.style.setProperty('--psdi-letter-spacing-base', "normal");
    } else {
        r.style.setProperty('--psdi-letter-spacing-base', space+"px");
    }
}

// Changes the line spacing for accessibility purposes.
function changeLineSpacing(event) {
    const space = $("#line").find(":selected").text();
    var text = $(".normalText, .middle");
    
    if (space=="Default") {
        r.style.setProperty('--ifm-line-height-base', "1.5");
    } else {
        r.style.setProperty('--ifm-line-height-base', space);
    }
}

// Changes the font size for accessibility purposes.
function changeFontSize(event) {
    const size = $("#size").find(":selected").text();

    if (size=="Default") {
        r.style.setProperty('--ifm-font-size-base', "1rem");
    } else {
        r.style.setProperty('--ifm-font-size-base', size+"px");
    }
}

// Changes the font weight for accessibility purposes.
function changeFontWeight(event) {
    const weight = $("#weight").find(":selected").text();

    if (weight=="Default") {
        r.style.setProperty('--ifm-font-weight-base', "normal");
    } else {
        r.style.setProperty('--ifm-font-weight-base', weight.toLowerCase());
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
    let new_colour;

    new_colour = colour.toLowerCase();

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

    switch (colour) {
        case 'Grey':
            r.style.setProperty('--ifm-background-color', "lightgrey");
            break;

        case 'Mustard':
            r.style.setProperty('--ifm-background-color', "#eddd6e");
            break;

        case 'Peach':
            r.style.setProperty('--ifm-background-color', "#edd1b0");
            break;

        case 'Lemon':
            r.style.setProperty('--ifm-background-color', "#f8fd89");
            break;

        default:
            r.style.setProperty('--ifm-background-color', "white");
            break;
    }
}

// Reverts all select boxes to 'Default'
function resetSelections(event) {
    $("#font, #size, #weight, #letter, #line, #dark-colour, #light-colour, #background").val('Default').change();
}

// Applies accessibility settings to the entire website.
function applySettings(event) {
    sessionStorage.setItem("font", r.style.getPropertyValue('--ifm-font-family-base'));
    sessionStorage.setItem("size", r.style.getPropertyValue('--ifm-font-size-base'));
    sessionStorage.setItem("weight", r.style.getPropertyValue('--ifm-font-weight-base'));
    sessionStorage.setItem("letter", r.style.getPropertyValue('--psdi-letter-spacing-base'));
    sessionStorage.setItem("line", r.style.getPropertyValue('--ifm-line-height-base'));
    sessionStorage.setItem("darkColour", r.style.getPropertyValue('--psdi-dark-text-color-body'));
    sessionStorage.setItem("lightColour", r.style.getPropertyValue('--psdi-light-text-color-body'));
    sessionStorage.setItem("back", r.style.getPropertyValue('--ifm-background-color'));

    sessionStorage.setItem("fontOpt", $("#font").find(":selected").text());
    sessionStorage.setItem("sizeOpt", $("#size").find(":selected").text());
    sessionStorage.setItem("weightOpt", $("#weight").find(":selected").text());
    sessionStorage.setItem("letterOpt", $("#letter").find(":selected").text());
    sessionStorage.setItem("lineOpt", $("#line").find(":selected").text());
    sessionStorage.setItem("darkColourOpt", $("#dark-colour").find(":selected").text());
    sessionStorage.setItem("lightColourOpt", $("#light-colour").find(":selected").text());
    sessionStorage.setItem("backOpt", $("#background").find(":selected").text());
    console.log(sessionStorage.getItem("backOpt"));

    alert("The settings have been applied to the entire website.");
}

