"""auth.py

This module contains the OpenID Connect and JSON Web Token handling.
"""

import os
import jwt
import json
import requests

from datetime import datetime
from cachetools.func import ttl_cache
from flask import Flask, abort, request, redirect
from urllib.parse import urlencode, unquote


# Authentication settings
KEYCLOAK_URL = os.getenv("KEYCLOAK_URL", None)
KEYCLOAK_REALM = os.getenv("KEYCLOAK_REALM", None)
KEYCLOAK_CLIENT_ID = os.getenv("KEYCLOAK_CLIENT_ID", None)
KEYCLOAK_SECRET = os.getenv("KEYCLOAK_SECRET", None)
KEYCLOAK_REDIRECT_URL = os.getenv("KEYCLOAK_REDIRECT_URL", None)
SESSION_TIMEOUT_SECONDS = int(os.getenv("SESSION_TIMEOUT_SECONDS", None))

user_keys = {}


@ttl_cache(ttl=60)
def get_keycloak_public_key():
    # Get JSON Web Key Set from Keycloak so we can verify tokens.
    # This needs to run periodically.

    jwks_url = f"{KEYCLOAK_URL}/realms/{KEYCLOAK_REALM}/protocol/openid-connect/certs"
    jwks = requests.get(jwks_url).json()
    public_keys = {}

    for key in jwks['keys']:
        public_keys[key['kid']] = jwt.algorithms.RSAAlgorithm.from_jwk(json.dumps(key))

    return public_keys


def get_login_url():

    query = {
        'client_id': KEYCLOAK_CLIENT_ID,
        'redirect_uri': KEYCLOAK_REDIRECT_URL,
        'response_type': 'code',
        'scope': 'openid',
    }

    return f"{KEYCLOAK_URL}/realms/{KEYCLOAK_REALM}/protocol/openid-connect/auth?{urlencode(query)}"


def get_logout_url():

    return "/logout"


def oidc_callback():
    # Get the _code_ parameter to use for Keycloak communication

    code = request.args.get('code')

    if not code:
        abort(400)

    # Make request to Keycloak for a id / access token
    keycloak_data = {
        "grant_type": "authorization_code",
        "code": code,
        "client_id": KEYCLOAK_CLIENT_ID,
        "client_secret": KEYCLOAK_SECRET,
        "redirect_uri": KEYCLOAK_REDIRECT_URL,
        "scope": "openid",
    }

    keycloak_url = f"{KEYCLOAK_URL}/realms/{KEYCLOAK_REALM}/protocol/openid-connect/token"

    token_data = requests.post(keycloak_url, data=keycloak_data).json()
    access_token = token_data.get("access_token")

    try:
        # Verify and decode the access token
        headers = jwt.get_unverified_header(access_token)
        public_key = get_keycloak_public_key().get(headers['kid'])

        decoded_access_token = jwt.decode(
            access_token,
            key=public_key,
            audience="account",
            algorithms=['RS256']
        )

        user_public_key_string = unquote(request.cookies.get('public_key'))

        user_public_key = json.loads(user_public_key_string)

        kid = user_public_key["kid"]

        user_keys[kid] = {
            "last_used": datetime.utcnow(),
            "access_token": decoded_access_token,
            "public_key": jwt.PyJWK.from_json(user_public_key_string)
        }

        return redirect("/")

    except jwt.InvalidTokenError as e:
        print(f"Failed to verify access token: {e}")
        abort(400)


def logout():

    authenticated_user = get_authenticated_user()

    if authenticated_user != None:

        sid = authenticated_user["sid"]

        for kid in user_keys.copy():

            if user_keys[kid]["access_token"]["sid"] == sid:

                del user_keys[kid]

    return redirect("/")


def get_authenticated_user():

    authenticated_user = None

    auth_token_string = request.cookies.get('auth_token')

    if auth_token_string != None:

        auth_token = unquote(auth_token_string)

        unverified_headers = jwt.get_unverified_header(auth_token)

        kid = unverified_headers['kid']

        if kid in user_keys:

            user_key = user_keys[kid]

            now = datetime.utcnow()
            elapsed = now - user_key["last_used"]

            if elapsed.seconds > SESSION_TIMEOUT_SECONDS:

                del user_keys[kid]

            else:

                try:

                    decoded_user_key = jwt.decode(
                        auth_token,
                        key=user_key["public_key"],
                        algorithms=['ES256']
                    )

                    authenticated_user = user_key["access_token"]
                    user_key["last_used"] = datetime.utcnow()

                except jwt.InvalidTokenError as e:

                    print(f"Failed to verify session token: {e}")

    return authenticated_user


def init_authentication(app: Flask):
    """Connect the provided Flask app to each of the post methods
    """

    app.route('/oidc_callback', methods=['GET'])(oidc_callback)
    app.route('/logout', methods=['GET'])(logout)
