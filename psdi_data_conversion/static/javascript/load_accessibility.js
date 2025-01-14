const r = document.querySelector(':root');
const s = getComputedStyle(document.documentElement);

// Check if the default values have been saved yet, and save them if not
if (s.getPropertyValue('--psdi-default-font') == "") {

  r.style.setProperty('--psdi-default-font', s.getPropertyValue('--ifm-font-family-base'));
  r.style.setProperty('--psdi-default-heading-font', s.getPropertyValue('--ifm-heading-font-family'));

  r.style.setProperty('--psdi-default-font-size', s.getPropertyValue('--ifm-font-size-base'));

  r.style.setProperty('--psdi-default-font-weight', s.getPropertyValue('--ifm-font-weight-base'));

  r.style.setProperty('--psdi-default-letter-spacing', s.getPropertyValue('--psdi-letter-spacing-base'));

  r.style.setProperty('--psdi-default-dark-text-color-body',
      s.getPropertyValue('--psdi-dark-text-color-body'));
  r.style.setProperty('--psdi-default-dark-text-color-heading',
      s.getPropertyValue('--psdi-dark-text-color-heading'));
  r.style.setProperty('--psdi-default-light-text-color-body',
      s.getPropertyValue('--psdi-light-text-color-body'));
  r.style.setProperty('--psdi-default-light-text-color-heading',
      s.getPropertyValue('--psdi-light-text-color-heading'));

  r.style.setProperty('--psdi-default-line-height', s.getPropertyValue('--ifm-line-height-base'));

  r.style.setProperty('--psdi-default-background-color', s.getPropertyValue('--ifm-background-color'));
  r.style.setProperty('--psdi-default-color-primary', s.getPropertyValue('--ifm-color-primary'));
}

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