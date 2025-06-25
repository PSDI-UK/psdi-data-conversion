"""auth.py

This module contains the OpenID Connect and JSON Web Token handling.
"""

import json
import sys
from datetime import datetime, timezone
from typing import Any
from urllib.parse import unquote, urlencode

import jwt
import requests
from cachetools.func import ttl_cache
from flask import Flask, abort, redirect, request

from psdi_data_conversion.gui.env import get_env

d_all_user_keys: dict[str, dict | str] = {}


@ttl_cache(ttl=60)
def get_keycloak_public_key():
    # Get JSON Web Key Set from Keycloak so we can verify tokens.
    # This needs to run periodically.

    env = get_env()

    jwks_url = f"{env.keycloak_url}/realms/{env.keycloak_realm}/protocol/openid-connect/certs"
    d_jwks: dict[str, dict | str] = requests.get(jwks_url).json()
    public_keys: dict[str, str] = {}

    for d_key in d_jwks['keys']:
        public_keys[d_key['kid']] = jwt.algorithms.RSAAlgorithm.from_jwk(json.dumps(d_key))

    return public_keys


def get_login_url():

    env = get_env()

    query = {
        'client_id': env.keycloak_client_id,
        'redirect_uri': env.keycloak_redirect_url,
        'response_type': 'code',
        'scope': 'openid',
    }

    return f"{env.keycloak_url}/realms/{env.keycloak_realm}/protocol/openid-connect/auth?{urlencode(query)}"


def get_logout_url():

    return "/logout"


def oidc_callback():
    # Get the _code_ parameter to use for Keycloak communication

    code = request.args.get('code')

    if not code:
        abort(400)

    env = get_env()

    # Make request to Keycloak for a id / access token
    keycloak_data = {
        "grant_type": "authorization_code",
        "code": code,
        "client_id": env.keycloak_client_id,
        "client_secret": env.get_keycloak_secret(),
        "redirect_uri": env.keycloak_redirect_url,
        "scope": "openid",
    }

    keycloak_url = f"{env.keycloak_url}/realms/{env.keycloak_realm}/protocol/openid-connect/token"

    d_token_data: dict[str, dict | str] = requests.post(keycloak_url, data=keycloak_data).json()
    access_token: str = d_token_data.get("access_token")

    if not access_token:
        print(f"ERROR: Access token not granted: {repr(d_token_data)}", file=sys.stderr)
        abort(400)

    try:
        # Verify and decode the access token
        d_header: dict[str, Any] = jwt.get_unverified_header(access_token)
        public_key = get_keycloak_public_key()[d_header['kid']]

        decoded_access_token: str = jwt.decode(
            access_token,
            key=public_key,
            audience="account",
            algorithms=['RS256']
        )

        user_public_key_string = unquote(request.cookies['public_key'])

        d_user_public_key: dict[str, dict | str] = json.loads(user_public_key_string)

        kid: str = d_user_public_key["kid"]

        d_all_user_keys[kid] = {
            "last_used": datetime.now(timezone.utc),
            "access_token": decoded_access_token,
            "public_key": jwt.PyJWK.from_json(user_public_key_string)
        }

        return redirect("/")

    except jwt.InvalidTokenError as e:
        print(f"ERROR: Failed to verify access token: {e}", file=sys.stderr)
        abort(400)


def logout():

    d_authenticated_user: dict[str, dict | str] | None = get_authenticated_user()

    if d_authenticated_user is None:

        return redirect("/")

    sid = d_authenticated_user["sid"]

    # Iterate over keys in a shallow copy so the iteration isn't disrupted when we delete a key
    for kid in d_all_user_keys.copy():

        if d_all_user_keys[kid]["access_token"]["sid"] == sid:

            del d_all_user_keys[kid]

    return redirect("/")


def get_authenticated_user():

    auth_token_string = request.cookies.get('auth_token')

    if auth_token_string is None:
        return None

    auth_token = unquote(auth_token_string)

    d_unverified_header: dict[str, dict | str] = jwt.get_unverified_header(auth_token)

    kid: str = d_unverified_header['kid']

    d_authenticated_user: dict[str, dict | str] | None = None

    if kid in d_all_user_keys:

        d_user_key: dict[str, dict | str] = d_all_user_keys[kid]

        now = datetime.now(timezone.utc)
        elapsed = now - d_user_key["last_used"]
        timeout = get_env().session_timeout_seconds

        if timeout > 0 and elapsed.seconds > timeout:

            del d_all_user_keys[kid]

        else:

            try:

                d_authenticated_user = d_user_key["access_token"]
                d_user_key["last_used"] = datetime.now(timezone.utc)

            except jwt.InvalidTokenError as e:

                print(f"Failed to verify session token: {e}", file=sys.stderr)

    return d_authenticated_user


def init_authentication(app: Flask):
    """Connect the provided Flask app to each of the post methods
    """

    app.route('/oidc_callback', methods=['GET'])(oidc_callback)
    app.route('/logout', methods=['GET'])(logout)
