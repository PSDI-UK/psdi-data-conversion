/*
  documentation.js
  Version 1.0, 26th June 2024

  This is the JavaScript which makes the documentation.htm gui work.
*/

const r = document.querySelector(':root');

$(document).ready(function() {

    const font = sessionStorage.getItem("font"),
          size = sessionStorage.getItem("size"),
          weight = sessionStorage.getItem("weight"),
          letter = sessionStorage.getItem("letter"),
          line = sessionStorage.getItem("line"),
          darkColour = sessionStorage.getItem("darkColour"),
          lightColour = sessionStorage.getItem("lightColour"),
          lightBack = sessionStorage.getItem("lightBack"),
          darkBack = sessionStorage.getItem("darkBack");

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

        r.style.setProperty('--ifm-background-color', lightBack);
        r.style.setProperty('--ifm-color-primary', darkBack);
    }
});

