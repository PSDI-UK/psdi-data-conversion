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