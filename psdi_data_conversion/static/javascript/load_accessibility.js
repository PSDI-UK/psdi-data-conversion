const r = document.querySelector(':root');
const s = getComputedStyle(document.documentElement);
const initMode = sessionStorage.getItem("mode");

function setDefault(default_varname, current_varname) {
  if (sessionStorage.getItem(default_varname) == null) {
    sessionStorage.setItem(default_varname, s.getPropertyValue('--' + current_varname))
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
let font = sessionStorage.getItem("font"),
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

export function loadProperty(current_varname, value) {
  if (value != null) {
    r.style.setProperty('--' + current_varname, value);
  }
}

function applyStoredAccessibility() {

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

}

if (font != null) {
  applyStoredAccessibility();
}

export function loadAccessibility() {
  $.get(`/load_accessibility/`)
    .done((data) => {

      const oData = JSON.parse(data);

      function getAndSave(key) {
        let value = oData[key];
        if (value != null) {
          sessionStorage.setItem(key, value);
          sessionStorage.setItem(key + "Opt", oData[key + "Opt"]);
          return value;
        } else {
          return sessionStorage.getItem(key);
        }
      }

      font = getAndSave("font");
      hfont = getAndSave("hfont");
      size = getAndSave("size");
      weight = getAndSave("weight");
      letter = getAndSave("letter");
      line = getAndSave("line");
      darkColour = getAndSave("darkColour");
      lightColour = getAndSave("lightColour");
      lightBack = getAndSave("lightBack");
      darkBack = getAndSave("darkBack");
      mode = getAndSave("mode");

      if (initMode == null) {
        // Only load accessibility settings if the light/dark mode hasn't been toggled
        applyStoredAccessibility();
      }
    });

}

loadAccessibility();

document.documentElement.setAttribute("data-theme", mode);

// Connect the color mode toggle button in the header - since we write the header directly in our templates, the
// psdi-common.js code doesn't connect the normal function, and so we use our own custom function here to play nicely
// with our accessibility settings

const LIGHT_MODE = "light";
const DARK_MODE = "dark";

export function setMode(new_mode = null) {

  // If not provide a mode, toggle between modes
  if (new_mode == null) {
    let currentMode = document.documentElement.getAttribute("data-theme");
    if (currentMode == DARK_MODE) {
      new_mode = LIGHT_MODE;
    } else {
      new_mode = DARK_MODE;
    }
  }

  if (new_mode == DARK_MODE) {

    loadProperty("psdi-dark-text-color-body", sessionStorage.getItem("psdi-default-light-text-color-body"));
    loadProperty("psdi-dark-text-color-heading", sessionStorage.getItem("psdi-default-light-text-color-heading"));
    loadProperty("psdi-light-text-color-body", sessionStorage.getItem("psdi-default-dark-text-color-body"));
    loadProperty("psdi-light-text-color-heading", sessionStorage.getItem("psdi-default-dark-text-color-heading"));

    loadProperty("ifm-background-color", s.getPropertyValue("--psdi-dm-bg-color-default"));
    loadProperty("ifm-color-primary", s.getPropertyValue("--psdi-dm-bg-color-primary"));

  } else if (new_mode == LIGHT_MODE) {

    loadProperty("psdi-dark-text-color-body", sessionStorage.getItem("psdi-default-dark-text-color-body"));
    loadProperty("psdi-dark-text-color-heading", sessionStorage.getItem("psdi-default-dark-text-color-heading"));
    loadProperty("psdi-light-text-color-body", sessionStorage.getItem("psdi-default-light-text-color-body"));
    loadProperty("psdi-light-text-color-heading", sessionStorage.getItem("psdi-default-light-text-color-heading"));

    loadProperty("ifm-background-color", sessionStorage.getItem("psdi-default-background-color"));
    loadProperty("ifm-color-primary", sessionStorage.getItem("psdi-default-color-primary"));

  } else {
    loadAccessibility();
  }

  document.documentElement.setAttribute("data-theme", new_mode);
  sessionStorage.setItem("mode", new_mode);
}

function toggleMode() {
  setMode();
}

const lModeToggleButton = document.querySelectorAll(".color-mode-toggle");
lModeToggleButton.forEach(function (modeToggleButton) {
  modeToggleButton.addEventListener("click", toggleMode);
});

// Load the settings for the current mode if it's already been toggled in this session
if (initMode != null) {
  setMode(initMode);
}