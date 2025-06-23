"""auth.py

This module contains the OpenID Connect and JSON Web Token handling.
"""

import json
import sys
from datetime import datetime
from urllib.parse import unquote, urlencode

import jwt
import requests
from cachetools.func import ttl_cache
from flask import Flask, abort, redirect, request

from psdi_data_conversion.gui.env import get_env

user_keys = {}


@ttl_cache(ttl=60)
def get_keycloak_public_key():
    # Get JSON Web Key Set from Keycloak so we can verify tokens.
    # This needs to run periodically.

    env = get_env()

    jwks_url = f"{env.keycloak_url}/realms/{env.keycloak_realm}/protocol/openid-connect/certs"
    jwks = requests.get(jwks_url).json()
    public_keys = {}

    for key in jwks['keys']:
        public_keys[key['kid']] = jwt.algorithms.RSAAlgorithm.from_jwk(json.dumps(key))

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
        print(f"Failed to verify access token: {e}", file=sys.stderr)
        abort(400)


def logout():

    authenticated_user = get_authenticated_user()

    if authenticated_user is not None:

        sid = authenticated_user["sid"]

        for kid in user_keys.copy():

            if user_keys[kid]["access_token"]["sid"] == sid:

                del user_keys[kid]

    return redirect("/")


def get_authenticated_user():

    authenticated_user = None

    auth_token_string = request.cookies.get('auth_token')

    if auth_token_string is not None:

        auth_token = unquote(auth_token_string)

        unverified_headers = jwt.get_unverified_header(auth_token)

        kid = unverified_headers['kid']

        if kid in user_keys:

            user_key = user_keys[kid]

            now = datetime.utcnow()
            elapsed = now - user_key["last_used"]
            timeout = get_env().session_timeout_seconds

            if timeout > 0 and elapsed.seconds > timeout:

                del user_keys[kid]

            else:

                try:

                    authenticated_user = user_key["access_token"]
                    user_key["last_used"] = datetime.utcnow()

                except jwt.InvalidTokenError as e:

                    print(f"Failed to verify session token: {e}", file=sys.stderr)

    return authenticated_user


def init_authentication(app: Flask):
    """Connect the provided Flask app to each of the post methods
    """

    app.route('/oidc_callback', methods=['GET'])(oidc_callback)
    app.route('/logout', methods=['GET'])(logout)
