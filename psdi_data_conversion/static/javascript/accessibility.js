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

function loadOption(jsName, cssSelector, changeFunc) {
    const opt = sessionStorage.getItem(jsName+"Opt");
    if (opt!=null)
        $(cssSelector).val(opt).change();
    $(cssSelector).change(changeFunc);
}

$(document).ready(function() {

    loadOption("font", "#font", changeFont);
    loadOption("size", "#size", changeFontSize);
    loadOption("weight", "#weight", changeFontWeight);
    loadOption("letter", "#letter", changeLetterSpacing);
    loadOption("line", "#line", changeLineSpacing);
    loadOption("darkColour", "#dark-colour", changeFontColourDark);
    loadOption("lightColour", "#light-colour", changeFontColourLight);
    loadOption("lightBack", "#light-background", changeLightBackground);
    loadOption("darkBack", "#dark-background", changeDarkBackground);

    $("#resetButton").click(resetSelections);
    $("#applyButton").click(applyAllSettings);
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
    ["#font", "#size", "#weight", "#letter", "#line", "#dark-colour", "#light-colour", "#light-background",
        "#dark-background"].forEach(function (selector) {
        // Don't trigger a change event if it's already on Default
        if ($(selector).find(":selected").val() != "Default")
            $(selector).val("Default").change();
    });
}

// Save a setting for one accessibility option to sessionStorage
function applySetting(jsName, cssSelector, cssVar) {

    // Check if set to default and not previously set, in which case don't save anything to storage
    let selectedVal = $(cssSelector).find(":selected").val();
    if (selectedVal=="Default" && sessionStorage.getItem(jsName)==null)
        return;

    sessionStorage.setItem(jsName, s.getPropertyValue(cssVar));
    sessionStorage.setItem(jsName+"Opt", selectedVal);
}

// Applies accessibility settings to the entire website.
function applyAllSettings(event) {
    applySetting("font", "#font", "--ifm-font-family-base");
    applySetting("hfont", "#font", "--ifm-heading-font-family");
    applySetting("size", "#size", "--ifm-font-size-base");
    applySetting("weight", "#weight", "--ifm-font-weight-base");
    applySetting("letter", "#letter", "--psdi-letter-spacing-base");
    applySetting("line", "#line", "--ifm-line-height-base");
    applySetting("darkColour", "#dark-colour", "--psdi-dark-text-color-body");
    applySetting("lightColour", "#light-colour", "--psdi-light-text-color-body");
    applySetting("lightBack", "#light-background", "--ifm-background-color");
    applySetting("darkBack", "#dark-background", "--ifm-color-primary");

    alert("The settings have been applied to the entire website.");
}

