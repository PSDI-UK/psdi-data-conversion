/*
  accessibility.js
  Version 1.0, 7th June 2024

  This is the JavaScript which makes the Accessibility gui work.
*/

const r = document.querySelector(':root');
const s = getComputedStyle(document.documentElement);

const LIGHT_MODE = "light";
const DARK_MODE = "dark";

function toggleMode() {
  let currentMode = document.documentElement.getAttribute("data-theme");
  let new_mode;

  if (currentMode == DARK_MODE) {
    new_mode = LIGHT_MODE;
  } else {
    new_mode = DARK_MODE;
  }

  document.documentElement.setAttribute("data-theme", new_mode);
  sessionStorage.setItem("mode", new_mode);
}

$(document).ready(function() {

    const fontOpt = sessionStorage.getItem("fontOpt"),
          sizeOpt = sessionStorage.getItem("sizeOpt"),
          weightOpt = sessionStorage.getItem("weightOpt"),
          letterOpt = sessionStorage.getItem("letterOpt"),
          lineOpt = sessionStorage.getItem("lineOpt"),
          darkColourOpt = sessionStorage.getItem("darkColourOpt"),
          lightColourOpt = sessionStorage.getItem("lightColourOpt"),
          lightBackOpt = sessionStorage.getItem("lightBackOpt"),
          darkBackOpt = sessionStorage.getItem("darkBackOpt");

    if (fontOpt != null) {

        $("#font").val(fontOpt).change();
        $("#size").val(sizeOpt).change();
        $("#weight").val(weightOpt).change();
        $("#letter").val(letterOpt).change();
        $("#line").val(lineOpt).change();
        $("#dark-colour").val(darkColourOpt).change();
        $("#light-colour").val(lightColourOpt).change();
        $("#light-background").val(lightBackOpt).change();
        $("#dark-background").val(darkBackOpt).change();
    }

    $("#font").change(changeFont);
    $("#size").change(changeFontSize);
    $("#weight").change(changeFontWeight);
    $("#letter").change(changeLetterSpacing);
    $("#line").change(changeLineSpacing);
    $("#dark-colour").change(changeFontColourDark);
    $("#light-colour").change(changeFontColourLight);
    $("#light-background").change(changeLightBackground);
    $("#dark-background").change(changeDarkBackground);
    $("#resetButton").click(resetSelections);
    $("#applyButton").click(applySettings);
});

// Changes the font for accessibility purposes
function changeFont(event) {
    const fontSelection = $("#font").find(":selected");
    const font = fontSelection.text().trim();

    if (font=="Default") {
        r.style.setProperty('--ifm-font-family-base', s.getPropertyValue('--psdi-default-font'));
        r.style.setProperty('--ifm-heading-font-family', s.getPropertyValue('--psdi-default-heading-font'));
    } else {
        // To avoid duplication of font settings, we retrieve the style to apply from what's applied to the font in the
        // selection box
        let fontFamily = fontSelection[0].style['font-family'];
        r.style.setProperty('--ifm-font-family-base', fontFamily);
        r.style.setProperty('--ifm-heading-font-family', fontFamily);
    }
}

// Changes the letter spacing for accessibility purposes.
function changeLetterSpacing(event) {
    const space = $("#letter").find(":selected").text();

    if (space == "Default") {
        r.style.setProperty('--psdi-letter-spacing-base', s.getPropertyValue('--psdi-default-letter-spacing'));
    } else {
        r.style.setProperty('--psdi-letter-spacing-base', space+"px");
    }
}

// Changes the line spacing for accessibility purposes.
function changeLineSpacing(event) {
    const space = $("#line").find(":selected").text();
    
    if (space=="Default") {
        r.style.setProperty('--ifm-line-height-base', s.getPropertyValue('--psdi-default-line-height'));
    } else {
        r.style.setProperty('--ifm-line-height-base', space);
    }
}

// Changes the font size for accessibility purposes.
function changeFontSize(event) {
    const size = $("#size").find(":selected").text();

    if (size=="Default") {
        r.style.setProperty('--ifm-font-size-base', s.getPropertyValue('--psdi-default-font-size'));
    } else {
        r.style.setProperty('--ifm-font-size-base', size+"px");
    }
}

// Changes the font weight for accessibility purposes.
function changeFontWeight(event) {
    const weight = $("#weight").find(":selected").text();

    if (weight=="Default") {
        r.style.setProperty('--ifm-font-weight-base', s.getPropertyValue('--psdi-default-font-weight'));
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

    if (colour==='Default') {
        r.style.setProperty('--psdi-'+lightOrDark+'-text-color-body',
            s.getPropertyValue('--psdi-default-'+lightOrDark+'-text-color-body'));
        r.style.setProperty('--psdi-'+lightOrDark+'-text-color-heading',
            s.getPropertyValue('--psdi-default-'+lightOrDark+'-text-color-heading'));
    } else {
        r.style.setProperty('--psdi-'+lightOrDark+'-text-color-body', colour);
        r.style.setProperty('--psdi-'+lightOrDark+'-text-color-heading', colour);
    }
}

// Changes the background colour for accessibility purposes.
function changeLightBackground(event) {
    const colour = $("#light-background").find(":selected").text();

    if (colour=="Default") {
        r.style.setProperty('--ifm-background-color', s.getPropertyValue('--psdi-default-background-color'));
    } else {
        r.style.setProperty('--ifm-background-color', colour);
    }
}

// Changes the background colour for accessibility purposes.
function changeDarkBackground(event) {
    const colour = $("#dark-background").find(":selected").text();

    if (colour=="Default") {
        r.style.setProperty('--ifm-color-primary', s.getPropertyValue('--psdi-default-color-primary'));
    } else {
        r.style.setProperty('--ifm-color-primary', colour);
    }
}

// Reverts all select boxes to 'Default'
function resetSelections(event) {
    $("#font, #size, #weight, #letter, #line, #dark-colour, #light-colour, #light-background, #dark-background").val('Default').change();
}

// Applies accessibility settings to the entire website.
function applySettings(event) {
    sessionStorage.setItem("font", s.getPropertyValue('--ifm-font-family-base'));
    sessionStorage.setItem("hfont", s.getPropertyValue('--ifm-heading-font-family'));
    sessionStorage.setItem("size", s.getPropertyValue('--ifm-font-size-base'));
    sessionStorage.setItem("weight", s.getPropertyValue('--ifm-font-weight-base'));
    sessionStorage.setItem("letter", s.getPropertyValue('--psdi-letter-spacing-base'));
    sessionStorage.setItem("line", s.getPropertyValue('--ifm-line-height-base'));
    sessionStorage.setItem("darkColour", s.getPropertyValue('--psdi-dark-text-color-body'));
    sessionStorage.setItem("lightColour", s.getPropertyValue('--psdi-light-text-color-body'));
    sessionStorage.setItem("lightBack", s.getPropertyValue('--ifm-background-color'));
    sessionStorage.setItem("darkBack", s.getPropertyValue('--ifm-color-primary'));

    sessionStorage.setItem("fontOpt", $("#font").find(":selected").val());
    sessionStorage.setItem("sizeOpt", $("#size").find(":selected").val());
    sessionStorage.setItem("weightOpt", $("#weight").find(":selected").val());
    sessionStorage.setItem("letterOpt", $("#letter").find(":selected").val());
    sessionStorage.setItem("lineOpt", $("#line").find(":selected").val());
    sessionStorage.setItem("darkColourOpt", $("#dark-colour").find(":selected").val());
    sessionStorage.setItem("lightColourOpt", $("#light-colour").find(":selected").val());
    sessionStorage.setItem("lightBackOpt", $("#light-background").find(":selected").val());
    sessionStorage.setItem("darkBackOpt", $("#dark-background").find(":selected").val());

    alert("The settings have been applied to the entire website.");
}

