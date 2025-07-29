/*
  accessibility.js
  Version 1.0, 7th June 2024

  This is the JavaScript which makes the Accessibility gui work.
*/

const r = document.querySelector(':root');
const s = getComputedStyle(document.documentElement);

const LIGHT_MODE = "light";
const DARK_MODE = "dark";

function loadOption(jsName, cssSelector, changeFunc) {
    const opt = sessionStorage.getItem(jsName + "Opt");
    if (opt != null)
        $(cssSelector).val(opt).change();
    $(cssSelector).change(changeFunc);
}

$(document).ready(function () {

    // Store default styles in the stylesheet
    r.style.setProperty('--psdi-default-font', sessionStorage.getItem('psdi-default-font'));
    r.style.setProperty('--psdi-default-font-size', sessionStorage.getItem('psdi-default-font-size'));
    r.style.setProperty('--psdi-default-heading-font', sessionStorage.getItem('psdi-default-heading-font'));
    r.style.setProperty('--psdi-default-letter-spacing', sessionStorage.getItem('psdi-default-letter-spacing'));
    r.style.setProperty('--psdi-default-font-weight', sessionStorage.getItem('psdi-default-font-weight'));
    r.style.setProperty('--psdi-default-light-text-color-body',
        sessionStorage.getItem('psdi-default-light-text-color-body'));
    r.style.setProperty('--psdi-default-light-text-color-heading',
        sessionStorage.getItem('psdi-default-light-text-color-heading'));
    r.style.setProperty('--psdi-default-dark-text-color-body',
        sessionStorage.getItem('psdi-default-dark-text-color-body'));
    r.style.setProperty('--psdi-default-dark-text-color-heading',
        sessionStorage.getItem('psdi-default-dark-text-color-heading'));
    r.style.setProperty('--psdi-default-background-color', sessionStorage.getItem('psdi-default-background-color'));
    r.style.setProperty('--psdi-default-color-primary', sessionStorage.getItem('psdi-default-color-primary'));

    // Set up selection boxes
    loadOption("font", "#font", changeFont);
    loadOption("size", "#size", changeFontSize);
    loadOption("weight", "#weight", changeFontWeight);
    loadOption("letter", "#letter", changeLetterSpacing);
    loadOption("line", "#line", changeLineSpacing);
    loadOption("darkColour", "#dark-colour", changeFontColourDark);
    loadOption("lightColour", "#light-colour", changeFontColourLight);
    loadOption("lightBack", "#light-background", changeLightBackground);
    loadOption("darkBack", "#dark-background", changeDarkBackground);

    // Connect buttons
    $("#resetButton").click(resetSelections);
    $("#saveButton").click(saveSettings);
    $("#dmRestoreButton").click(restoreAllSettings);
    $("#lmRestoreButton").click(restoreAllSettings);
});

// Changes the font for accessibility purposes
function changeFont(event) {
    const fontSelection = $("#font").find(":selected");
    const font = fontSelection.text().trim();

    if (font == "Default") {
        r.style.setProperty('--ifm-font-family-base', sessionStorage.getItem('psdi-default-font'));
        r.style.setProperty('--ifm-heading-font-family', sessionStorage.getItem('psdi-default-heading-font'));
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
        r.style.setProperty('--psdi-letter-spacing-base', sessionStorage.getItem('psdi-default-letter-spacing'));
    } else {
        r.style.setProperty('--psdi-letter-spacing-base', space + "px");
    }
}

// Changes the line spacing for accessibility purposes.
function changeLineSpacing(event) {
    const space = $("#line").find(":selected").text();

    if (space == "Default") {
        r.style.setProperty('--ifm-line-height-base', sessionStorage.getItem('psdi-default-line-height'));
    } else {
        r.style.setProperty('--ifm-line-height-base', space);
    }
}

// Changes the font size for accessibility purposes.
function changeFontSize(event) {
    const size = $("#size").find(":selected").text();

    if (size == "Default") {
        r.style.setProperty('--ifm-font-size-base', sessionStorage.getItem('psdi-default-font-size'));
    } else {
        r.style.setProperty('--ifm-font-size-base', size + "px");
    }
}

// Changes the font weight for accessibility purposes.
function changeFontWeight(event) {
    const weight = $("#weight").find(":selected").text();

    if (weight == "Default") {
        r.style.setProperty('--ifm-font-weight-base', sessionStorage.getItem('psdi-default-font-weight'));
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

function changeFontColour(event, lightOrDark = "dark") {

    const colour = $("#" + lightOrDark + "-colour").find(":selected").text();

    if (colour === 'Default') {
        r.style.setProperty('--psdi-' + lightOrDark + '-text-color-body',
            sessionStorage.getItem('psdi-default-' + lightOrDark + '-text-color-body'));
        r.style.setProperty('--psdi-' + lightOrDark + '-text-color-heading',
            sessionStorage.getItem('psdi-default-' + lightOrDark + '-text-color-heading'));
    } else {
        r.style.setProperty('--psdi-' + lightOrDark + '-text-color-body', colour);
        r.style.setProperty('--psdi-' + lightOrDark + '-text-color-heading', colour);
    }
}

// Changes the background colour for accessibility purposes.
function changeLightBackground(event) {
    const colour = $("#light-background").find(":selected").text();

    if (colour == "Default") {
        r.style.setProperty('--ifm-background-color', sessionStorage.getItem('psdi-default-background-color'));
    } else {
        r.style.setProperty('--ifm-background-color', colour);
    }
}

// Changes the background colour for accessibility purposes.
function changeDarkBackground(event) {
    const colour = $("#dark-background").find(":selected").text();

    if (colour == "Default") {
        r.style.setProperty('--ifm-color-primary', sessionStorage.getItem('psdi-default-color-primary'));
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

/**
 * Get the selected value for a selection box and force a change to ensure it's in effect
 * 
 * @param {*} cssSelector The selector for the selection box, e.g. "#font"
 * @returns {string}
 */
function getAndRestoreSetting(cssSelector) {

    let selector = $(cssSelector).find(":selected");
    let selectedVal = selector.val();

    // Just in case something went wrong and the stored options and values aren't aligned, force a change to the
    // selected option to change the value to match
    $(selector).val(selectedVal).change();

    return selectedVal;
}

// Save a setting for one accessibility option to sessionStorage
function applySetting(jsName, cssSelector, cssVar, settingsData) {

    const selectedVal = getAndRestoreSetting(cssSelector);

    let val = s.getPropertyValue(cssVar);

    settingsData[jsName] = val;
    settingsData[jsName + "Opt"] = selectedVal;

    // Check if set to default and not previously set, in which case don't save anything to storage
    if (selectedVal == "Default" && sessionStorage.getItem(jsName) == null)
        return;

    sessionStorage.setItem(jsName, val);
    sessionStorage.setItem(jsName + "Opt", selectedVal);

}

// Applies accessibility settings to the entire website.
function saveSettings(event, show_alert = true) {

    let settingsData = new Object();

    applySetting("font", "#font", "--ifm-font-family-base", settingsData);
    applySetting("hfont", "#font", "--ifm-heading-font-family", settingsData);
    applySetting("size", "#size", "--ifm-font-size-base", settingsData);
    applySetting("weight", "#weight", "--ifm-font-weight-base", settingsData);
    applySetting("letter", "#letter", "--psdi-letter-spacing-base", settingsData);
    applySetting("line", "#line", "--ifm-line-height-base", settingsData);
    applySetting("darkColour", "#dark-colour", "--psdi-dark-text-color-body", settingsData);
    applySetting("lightColour", "#light-colour", "--psdi-light-text-color-body", settingsData);
    applySetting("lightBack", "#light-background", "--ifm-background-color", settingsData);
    applySetting("darkBack", "#dark-background", "--ifm-color-primary", settingsData);

    $.post(`/save_accessibility/`, {
        'data': JSON.stringify(settingsData)
    })
        .done(() => {
            if (show_alert) {
                alert("Your accessibility settings have been saved. If you accidentally save settings which are " +
                    "unreadable and can't find the \"Reset\" button, you can restore the default settings by " +
                    "deleting this site's cookie in your browser's settings.");
            }
        })
        .fail(function (e) {
            alert("ERROR: Could not save accessibility settings. Your settings should still persist for this " +
                "session, but will not be restored for future sessions. Please use the \"Contact\" link in the " +
                "header to report this error to us.");

            // For debugging
            console.error("Error saving accessibility settings", e.status, e.responseText);
        });


}

// Restores the ability to manage settings after using the mode toggle button
import { setMode, loadAccessibility } from "./load_accessibility.js";
function restoreAllSettings(event) {
    setMode("disable");
    loadAccessibility();

    // We save settings here to apply them and save them into session storage
    saveSettings(event, false);
}
