// data.js

const response = await fetch("/static/data/data.json");

let data = null;
try {
    data = JSON.parse(await response.text());
} catch (err) {
    console.log("ERROR: Database could not be loaded, error: " + err)
}

/**
 * Indicates whether or not the database was successfully loaded
 * @returns {boolean}
 */
export function databaseLoaded() {
    return data !== null;
}

const collator = new Intl.Collator();

// This compare function takes multiple pairs of strings so that it can be
// used to sort with multiple levels.

function compare(...pairs) {

    for (const pair of pairs) {

        const pairResult = collator.compare(pair[0], pair[1]);

        if (pairResult !== 0) {
            return pairResult;
        }
    }

    return 0;
}

async function getFormatByExtensionAndNote(extension, note) {

    const formats = data.formats.filter(format =>
        (format.extension === extension) && (format.note === note)
    );

    return formats[0];
}

export async function getInputFormats() {

    const inFormatIds = new Set(data.converts_to.map(record => record.in_id));

    const inFormats = data.formats.filter(format => inFormatIds.has(format.id));

    return inFormats.sort((a, b) => compare([a.extension, b.extension], [a.note, b.note]))
}

export async function getOutputFormats() {

    const outFormatIds = new Set(data.converts_to.map(record => record.out_id));

    const outFormats = data.formats.filter(format => outFormatIds.has(format.id));

    return outFormats.sort((a, b) => compare([a.extension, b.extension], [a.note, b.note]))
}

/**
 * Gets the ID for a format, given its extension and note
 * 
 * @param {string} extension - The extension of the format, e.g. 'pdb'
 * @param {string} note - The note of the format, e.g. 'Protein Databank'
 * @returns {(int|null)} - The ID of the format if found, or else null
 */
export async function getFormatId(extension, note) {

    var format = (data.formats.filter(format => (format.extension === extension) && (format.note === note)));

    if (format === undefined) {
        return null;
    }

    return format[0].id;
}

export async function getOutputFormatsForInputFormat(inExtension, inNote) {

    const inputFormat = (data.formats.filter(format => (format.extension === inExtension) && (format.note === inNote)))[0];

    if (inputFormat === undefined) {
        return [];
    }

    const outFormatIds = new Set(data.converts_to.filter(record => record.in_id === inputFormat.id).map(record => record.out_id));

    const outFormats = data.formats.filter(format => outFormatIds.has(format.id));

    return outFormats.sort((a, b) => compare([a.extension, b.extension], [a.note, b.note]))
}

export async function getInputFormatsForOutputFormat(outExtension, outNote) {

    const outputFormat = (data.formats.filter(format => (format.extension === outExtension) && (format.note === outNote)))[0];

    if (outputFormat === undefined) {
        return [];
    }

    const inFormatIds = new Set(data.converts_to.filter(record => record.out_id === outputFormat.id).map(record => record.in_id));

    const inFormats = data.formats.filter(format => inFormatIds.has(format.id));

    return inFormats.sort((a, b) => compare([a.extension, b.extension], [a.note, b.note]))
}

export async function getConverters(inExtension, inNote, outExtension, outNote) {

    const inputFormat = (data.formats.filter(format => (format.extension === inExtension) && (format.note === inNote)))[0];
    const outputFormat = (data.formats.filter(format => (format.extension === outExtension) && (format.note === outNote)))[0];

    if ((inputFormat === undefined) || (outputFormat === undefined)) {
        return [];
    }

    const convertersById = new Map(data.converters.map(converter => [converter.id, converter]));

    const conversions = data.converts_to.filter(record => (record.in_id === inputFormat.id) && record.out_id === outputFormat.id);

    const convertersWithDegreeOfSuccess = conversions.map(conversion => ({
        name: convertersById.get(conversion.converters_id).name
    }));

    return convertersWithDegreeOfSuccess.sort((a, b) => compare(a.name, b.name));
}

export async function getConverterByName(name) {
    return (data.converters.filter(converter => (converter.name === name)))[0];
}

export async function getAllFormats() {
    return data.formats.sort((a, b) => compare([a.extension, b.extension], [a.note, b.note]));
}

export async function getInputFlags(extension, note) {

    const format = (data.formats.filter(format => (format.extension === extension) && (format.note === note)))[0];

    if (format !== undefined) {

        const obFlagsIn = new Set(data.obformat_to_flags_in.filter(entry => entry.formats_id === format.id).map(entry => entry.obflags_in_id));

        return data.obflags_in.filter(entry => obFlagsIn.has(entry.id)).sort((a, b) => compare([a.flag, b.flag]));
    }
}

export async function getOutputFlags(extension, note) {

    const format = (data.formats.filter(format => (format.extension === extension) && (format.note === note)))[0];

    if (format !== undefined) {

        const obFlagsOut = new Set(data.obformat_to_flags_out.filter(entry => entry.formats_id === format.id).map(entry => entry.obflags_out_id));

        return data.obflags_out.filter(entry => obFlagsOut.has(entry.id)).sort((a, b) => compare([a.flag, b.flag]));
    }
}

export async function getInputArgFlags(extension, note) {

    const format = await getFormatByExtensionAndNote(extension, note);

    if (format === undefined) {
        throw `Can't find format for ${extension} and ${note}`;
    }

    const argFlagEntries = data.obformat_to_argflags_in
        .filter(entry => entry.formats_id === format.id)
        .map(entry => entry.obargflags_in_id);

    const argFlags = data.obargflags_in
        .filter(entry => argFlagEntries.indexOf(entry.id) !== -1);

    return argFlags.sort((a, b) => compare([a.flag, b.flag]));
}

export async function getOutputArgFlags(extension, note) {

    const format = await getFormatByExtensionAndNote(extension, note);

    if (format === undefined) {
        throw `Can't find format for ${extension} and ${note}`;
    }

    const argFlagEntries = data.obformat_to_argflags_out
        .filter(entry => entry.formats_id === format.id)
        .map(entry => entry.obargflags_out_id);

    const argFlags = data.obargflags_out
        .filter(entry => argFlagEntries.indexOf(entry.id) !== -1);

    return argFlags.sort((a, b) => compare([a.flag, b.flag]));
}

export async function getLevelChemInfo(inExtension, inNote, outExtension, outNote) {

    const inputFormat = (data.formats.filter(format => (format.extension === inExtension) && (format.note === inNote)))[0],
        outputFormat = (data.formats.filter(format => (format.extension === outExtension) && (format.note === outNote)))[0];

    return [inputFormat, outputFormat];
}
