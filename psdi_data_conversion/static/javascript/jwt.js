// IndexedDB routines for key storage.

const authDatabaseName = "authKeys";

const issuedAtAdjustedSeconds = -5;

function getDatabase() {

    return new Promise((resolve, reject) => {

        const request = indexedDB.open(authDatabaseName, 1);

        request.onsuccess = event => resolve(event.target.result);
        request.onupgradeneeded = event => event.target.result.createObjectStore("keys", { keyPath: "id" });
        request.onerror = reject;
    });
}

function putItemInDatabase(database, store, item) {

    return new Promise((resolve, reject) => {

        const transaction = database.transaction([store], "readwrite");
        const objectStore = transaction.objectStore(store);

        const req = objectStore.put(item);

        req.onsuccess = resolve;
        req.onerror = reject;
    });
}

function getItemFromDatabase(database, store, key) {

    return new Promise((resolve, reject) => {

        const transaction = database.transaction([store], "readonly");
        const objectStore = transaction.objectStore(store);

        const req = objectStore.get(key);

        req.onsuccess = event => resolve(event.target.result);
        req.onerror = reject;
    });
}

// ECDSA-signed JWT implemention.

const subtle = self.crypto.subtle;
const textEncoder = new TextEncoder();

// These are the ECDSA curve/hash combinations described in RFC 7518.

const algorithmDetails = {
    'ES256': { crv: 'P-256', hash: 'SHA-256' },
    'ES384': { crv: 'P-384', hash: 'SHA-384' },
    'ES512': { crv: 'P-521', hash: 'SHA-512' },
};

function convertUint8ArrayToBase64Url(array) {

    const outChars = [
        'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H',
        'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P',
        'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X',
        'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f',
        'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n',
        'o', 'p', 'q', 'r', 's', 't', 'u', 'v',
        'w', 'x', 'y', 'z', '0', '1', '2', '3',
        '4', '5', '6', '7', '8', '9', '-', '_'
    ];

    const numSextets = Math.ceil(array.length * 8 / 6);

    let result = "";

    for (let arrayIndex = 0, sextetIndex = 0; sextetIndex < numSextets; sextetIndex++) {

        let charIndex;

        switch (sextetIndex % 4) {

            case 0:

                charIndex = (array[arrayIndex] & 0xfc) >> 2;
                break;

            case 1:

                charIndex = ((array[arrayIndex] & 0x03) << 4) | ((array[arrayIndex + 1] & 0xf0) >> 4);
                arrayIndex++;
                break;

            case 2:

                charIndex = ((array[arrayIndex] & 0x0f) << 2) | ((array[arrayIndex + 1] & 0xc0) >> 6);
                arrayIndex++;
                break;

            case 3:

                charIndex = array[arrayIndex] & 0x3f;
                arrayIndex++;
                break;
        }

        result += outChars[charIndex];
    }

    return result;
}

function getAlgorithmDetails(jsonWebKey) {

    for (const algorithm in algorithmDetails) {

        const details = algorithmDetails[algorithm];

        if (jsonWebKey.crv === details.crv) {
            return Object.assign({ name: algorithm }, details);
        }
    }

    throw "Unsupported JWK";
}

async function generateKeyPair(algorithm) {

    const details = algorithmDetails[algorithm];

    if (details === undefined) {
        throw `Unsupported algorithm: ${algorithm}`;
    }

    const keyType = {
        name: 'ECDSA',
        namedCurve: details.crv
    };

    return subtle.generateKey(keyType, false, ['sign', 'verify']);
}

async function getThumbprint(publicKey) {

    const details = getAlgorithmDetails(publicKey);

    const digestDetails = textEncoder.encode(JSON.stringify({
        crv: publicKey.crv,
        kty: 'EC',
        x: publicKey.x,
        y: publicKey.y
    }));

    const hash = await subtle.digest({ name: details.hash }, digestDetails);

    return convertUint8ArrayToBase64Url(new Uint8Array(hash));
}

export async function hasKeyInfo(keyId) {

    const database = await getDatabase();

    let keyInfo = await getItemFromDatabase(database, "keys", keyId);

    return keyInfo !== undefined;
}

export async function getKeyInfo(keyId) {

    const database = await getDatabase();

    let keyInfo = await getItemFromDatabase(database, "keys", keyId);

    if (keyInfo === undefined) {

        const keyPair = await generateKeyPair('ES256');
        const exportedPublicKey = await subtle.exportKey('jwk', keyPair.publicKey);
        const thumbprint = await getThumbprint(exportedPublicKey);

        keyInfo = {
            id: keyId,
            publicKey: keyPair.publicKey,
            privateKey: keyPair.privateKey,
            exportedPublicKey: exportedPublicKey,
            thumbprint: thumbprint
        };

        putItemInDatabase(database, "keys", keyInfo);
    }

    return keyInfo;
}

export async function sign(claims, keyInfo, expireSeconds) {

    const privateKey = keyInfo.privateKey;
    const exportedPublicKey = keyInfo.exportedPublicKey;
    const thumbprint = keyInfo.thumbprint;

    const details = getAlgorithmDetails(exportedPublicKey);

    const headers = {
        typ: 'JWT',
        alg: details.name,
        kid: thumbprint
    }

    const keyType = {
        name: 'ECDSA',
        namedCurve: exportedPublicKey.crv,
    };

    const signatureType = {
        name: 'ECDSA',
        hash: { name: details.hash }
    };

    claims.iat = Math.floor(Date.now() / 1000) + issuedAtAdjustedSeconds;

    if (expireSeconds !== undefined) {
        claims.exp = claims.iat + expireSeconds;
    }

    const encodedHeaders = convertUint8ArrayToBase64Url(textEncoder.encode(JSON.stringify(headers)));
    const encodedPayload = convertUint8ArrayToBase64Url(textEncoder.encode(JSON.stringify(claims)));

    const dataToSign = textEncoder.encode(`${encodedHeaders}.${encodedPayload}`);

    const signature = await subtle.sign(signatureType, privateKey, dataToSign);

    const encodedSignature = convertUint8ArrayToBase64Url(new Uint8Array(signature));

    const jws = `${encodedHeaders}.${encodedPayload}.${encodedSignature}`;

    return jws;
}
