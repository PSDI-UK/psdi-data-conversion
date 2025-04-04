/**
 * @file common.js
 * @date 2025-02-14
 * @author Bryan Gillis
 */

export function initDirtyForms() {
  $("form.gui").dirtyForms();
}

export function cleanDirtyForms() {
  $('form.gui').dirtyForms('setClean');
}

export function dirtyDirtyForms() {
  $('form.gui').dirtyForms('setDirty');
}

export function enableDirtyForms() {
  $('form.gui').removeClass($.DirtyForms.ignoreClass);
}

export function disableDirtyForms() {
  $('form.gui').addClass($.DirtyForms.ignoreClass);
}


/**
 * Gets whether or not the app is operating in "Service mode"
 * 
 * This is the mode used for the public web app.
 *
 * @return {bool} True indicates service mode, False indicates local mode
 */
export function getServiceMode() {
  return sessionStorage.getItem("service_mode");
}

/**
 * Sets the service mode for the CSS document of the current page
 */
export function loadServiceMode() {
  document.documentElement.setAttribute("service-mode", getServiceMode());
}

// Set the service mode variable for each page so that only appropriate elements are shown
loadServiceMode();

/**
 * Gets whether or not the app is operating in "Production mode"
 * 
 * This is the mode used in staging/production deployments, but not dev deployments
 *
 * @return {bool} True indicates production mode, False indicates dev mode
 */
export function getProductionMode() {
  return sessionStorage.getItem("production_mode");
}

/**
 * Sets the production mode for the CSS document of the current page
 */
export function loadProductionMode() {
  document.documentElement.setAttribute("production-mode", getProductionMode());
}

// Set the production mode variable for this page so that only appropriate elements are shown
loadProductionMode();