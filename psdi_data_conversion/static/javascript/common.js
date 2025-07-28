/**
 * @file common.js
 * @date 2025-02-14
 * @author Bryan Gillis
 */

const AUTH_REFRESH = 60; // Seconds

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

async function refreshAuthCookie(repeat) {

  const jwt = await import("/static/javascript/jwt.js");

  const keyInfo = await jwt.getKeyInfo("keys");

  if (keyInfo) {

    const auth_token = await jwt.sign({ "type": "refreshAuth" }, keyInfo, AUTH_REFRESH * 2);

    document.cookie = `auth_token=${auth_token}`;

    if (repeat) {
      setTimeout(refreshAuthCookie, AUTH_REFRESH * 1000);
    }
  }
}

window.addEventListener("load", function () {

  const loginLink = document.querySelector("a#loginLink");

  // Add public key cookie when the user clicks login.

  if (loginLink !== null) {

    loginLink.addEventListener("click", async function (event) {

      const hasPublicKey = document.cookie.split(";").some((item) => item.trim().startsWith("public_key="));

      if (hasPublicKey === false) {

        event.preventDefault();

        const jwt = await import("/static/javascript/jwt.js");

        const keyInfo = await jwt.getKeyInfo("keys");

        const publicKey = {
          "kid": keyInfo.thumbprint,
          "crv": keyInfo.exportedPublicKey.crv,
          "kty": keyInfo.exportedPublicKey.kty,
          "x": keyInfo.exportedPublicKey.x,
          "y": keyInfo.exportedPublicKey.y
        };

        document.cookie = `public_key=${encodeURIComponent(JSON.stringify(publicKey))}`;

        await refreshAuthCookie();

        // loginLink.click();

        return;
      }
    });
  }
});

refreshAuthCookie(true);

$(document).ready(function () {

  // Enable all tooltips on the page
  const tooltipTriggerList = document.querySelectorAll('[data-bs-toggle="tooltip"]');
  const tooltipList = [...tooltipTriggerList].map(tooltipTriggerEl => new bootstrap.Tooltip(tooltipTriggerEl));

});