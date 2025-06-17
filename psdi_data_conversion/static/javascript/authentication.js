// authentication.js

const auth_refresh = 60;

async function refreshAuthCookie(repeat) {

    const jwt = await import("/static/javascript/jwt.js");

    const keyInfo = await jwt.getKeyInfo("keys");

    if (keyInfo) {

        const auth_token = await jwt.sign({ "type": "refreshAuth" }, keyInfo, auth_refresh * 2);

        document.cookie = `auth_token=${auth_token}`;

        if (repeat) {
            setTimeout(refreshAuthCookie, auth_refresh * 1000);
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
