/**
 * @file common.js
 * @date 2025-02-14
 * @author Bryan Gillis
 */

// Connect the color mode toggle button in the header - since we write the header directly in our templates, the
// psdi-common.js code doesn't run this automatically out of caution
import { connectModeToggleButton } from "./psdi-common.js";
connectModeToggleButton();

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