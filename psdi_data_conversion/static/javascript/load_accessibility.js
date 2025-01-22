const r = document.querySelector(':root');
const s = getComputedStyle(document.documentElement);

function setDefault(default_varname, current_varname) {
  if (s.getPropertyValue('--' + default_varname) == "") {
    r.style.setProperty('--' + default_varname, s.getPropertyValue('--' + current_varname))
  }
}

setDefault("psdi-default-font", "ifm-font-family-base");
setDefault("psdi-default-heading-font", "ifm-heading-font-family");

setDefault("psdi-default-font-size", "ifm-font-size-base");

setDefault("psdi-default-font-weight", "ifm-font-weight-base");

setDefault("psdi-default-letter-spacing", "psdi-letter-spacing-base");

setDefault("psdi-default-dark-text-color-body", "psdi-dark-text-color-body");
setDefault("psdi-default-dark-text-color-heading", "psdi-dark-text-color-heading");
setDefault("psdi-default-light-text-color-body", "psdi-light-text-color-body");
setDefault("psdi-default-light-text-color-heading", "psdi-light-text-color-heading");

setDefault("psdi-default-line-height", "ifm-line-height-base");

setDefault("psdi-default-background-color", "ifm-background-color");
setDefault("psdi-default-color-primary", "ifm-color-primary");

// Load values from session storage
const font = sessionStorage.getItem("font"),
  hfont = sessionStorage.getItem("hfont"),
  size = sessionStorage.getItem("size"),
  weight = sessionStorage.getItem("weight"),
  letter = sessionStorage.getItem("letter"),
  line = sessionStorage.getItem("line"),
  darkColour = sessionStorage.getItem("darkColour"),
  lightColour = sessionStorage.getItem("lightColour"),
  lightBack = sessionStorage.getItem("lightBack"),
  darkBack = sessionStorage.getItem("darkBack"),
  mode = sessionStorage.getItem("mode");

function loadProperty(current_varname, value) {
  if (value != null) {
    r.style.setProperty('--' + current_varname, value);
  }
}

loadProperty("ifm-font-family-base", font);
loadProperty("ifm-heading-font-family", hfont);

loadProperty("ifm-font-size-base", size);

loadProperty("ifm-font-weight-base", weight);

loadProperty("psdi-letter-spacing-base", letter);

loadProperty("psdi-dark-text-color-body", darkColour);
loadProperty("psdi-dark-text-color-heading", darkColour);
loadProperty("psdi-light-text-color-body", lightColour);
loadProperty("psdi-light-text-color-heading", lightColour);

loadProperty("ifm-line-height-base", line);

loadProperty("ifm-background-color", lightBack);
loadProperty("ifm-color-primary", darkBack);

if (font != null) {

  r.style.setProperty('--ifm-font-family-base', font);
  r.style.setProperty('--ifm-heading-font-family', hfont);

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

document.documentElement.setAttribute("data-theme", mode);