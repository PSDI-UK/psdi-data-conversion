// This file provides common functions related to the PSDI common assets. Most of these are run automatically as
// appropriate, except the first batch of functions exported below. These are exposed to the user so that they can
// customise aspects of the common HTML header

const DEFAULT_TITLE = "";
const DEFAULT_BRAND_LINK_TARGET = "./";
const DEFAULT_HEADER_LINKS_SOURCE = "./header-links.html";

let title = DEFAULT_TITLE;
let brandLinkTarget = DEFAULT_BRAND_LINK_TARGET;
let headerLinksSource = DEFAULT_HEADER_LINKS_SOURCE;

export function setTitle(s) {
  // Public function for the user to set the site title that will appear in the header, to the right of the PSDI logo
  title = s;
}

export function setBrandLinkTarget(s) {
  // Public function for the user to set the target that clicking on the PSDI brand should link to
  brandLinkTarget = s;
}

export function setHeaderLinksSource(s) {
  // Public function to set the name of an HTML file containing the links to appear on the right side of the header
  // for a given page
  headerLinksSource = s;
}

const LIGHT_MODE = "light";
const DARK_MODE = "dark";

// Load color mode from session storage and apply it
let mode = sessionStorage.getItem("mode");
if (!mode) {
  mode = LIGHT_MODE;
}
document.documentElement.setAttribute("data-theme", mode);

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

export function connectModeToggleButton() {
  // Connect the mode toggle function to the button
  const lModeToggleButton = document.querySelectorAll(".color-mode-toggle");
  lModeToggleButton.forEach(function (modeToggleButton) {
    modeToggleButton.addEventListener("click", toggleMode);
  });
}

// Counter for elements that need to be loaded - each we request loading will increment this by 1
let loadSteps = 0;

function finalizeLoad() {
  // Decrement the load steps and check if all steps are finished. If so, remove the cover
  --loadSteps;
  if (loadSteps <= 0) {
    $("#cover").hide();
  }
}

export function addHeaderLinks() {
  // We want to load in the links, but preserve the existing mode toggle button alongside them, so this function
  // handles saving it and re-adding it

  let headerLinksParent = $("#psdi-header .navbar__items--right");
  let modeToggle = $("#psdi-header .color-mode-toggle");

  headerLinksParent.load(headerLinksSource,
    function (_response, status, _xhr) {
      if (status == "error") {
        headerLinksParent[0].textContent = "ERROR: Could not load header links";
      }
      headerLinksParent[0].appendChild(modeToggle[0]);
      connectModeToggleButton();
      finalizeLoad();
    });
}

$(document).ready(function () {

  // Start fading out the cover over one second as a failsafe in case something goes wrong and it never gets removed
  $("#cover").fadeOut(1000);

  // Count the elements we'll need to load first, to avoid prematurely removing the cover
  // We load an element only if it's a pure stub with no children; otherwise we assume it is intended to be used
  // as-is and not overwritten by a load

  const headerStub = $("#psdi-header");
  let loadHeader = false;
  if (headerStub.length > 0 && headerStub[0].childNodes.length == 0) {
    loadHeader = true;
    ++loadSteps;
  }

  const footerStub = $("#psdi-footer");
  let loadFooter = false;
  if (footerStub.length > 0 && footerStub[0].childNodes.length == 0) {
    loadFooter = true;
    ++loadSteps;
  }

  // Load only if the header stub has no children
  if (loadHeader) {
    $("#psdi-header").load("./psdi-common-header.html",
      function (_response, status, _xhr) {
        if (status != "error") {
          $("#psdi-header a.navbar__brand")[0].href = brandLinkTarget;
          $("#psdi-header .navbar__title")[0].textContent = title;
          addHeaderLinks();
        } else {
          $("#psdi-header")[0].textContent = "ERROR: Could not load page header";
          connectModeToggleButton();
          finalizeLoad();
        }
      });
  }

  // Load only if the footer stub has no children
  if (loadFooter) {
    $("#psdi-footer").load("./psdi-common-footer.html",
      function (_response, status, _xhr) {
        if (status == "error") {
          $("#psdi-footer")[0].textContent = "ERROR: Could not load page footer";
        }
        finalizeLoad();
      });
  }

  if (loadSteps==0)
    finalizeLoad();
  
});
