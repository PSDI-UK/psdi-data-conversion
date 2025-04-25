/*
  convert_common.js
  Version 1.0, 17th December 2024
*/

import { disableDirtyForms, enableDirtyForms, initDirtyForms } from "./common.js";

const SECOND = 1000; // Milliseconds
const CONVERT_TIMEOUT = 60 * SECOND;
const MEGABYTE = 1024 * 1024;

const ZIP_EXT = "zip"
const TAR_EXT = "tar"
const GZ_EXT = "gz"
const BZ_EXT = "bz"
const XZ_EXT = "gz"
const TARGZ_EXT = "tar.gz"
const TARBZ_EXT = "tar.bz"
const TARXZ_EXT = "tar.xz"

// Short list of safe allowed characters:
// \w: All letters and digits
// \s: All whitespace characters
// .: Period
// \-: Hyphen
// :: Colon
// +: Plus symbol
// *: Asterisk
// =: Equals sign
// $: Dollar sign
// /: Forward-slash
// \\: Backslash
export const SAFE_CHAR_REGEX = "[\\w\\s.\\-:+*=$\\/\\\\]*"

// Whether or not file extensions will be checked
let extCheck = true;

// Whether or not the user wants to also get the log file
let requestLog = false;

var token = "",
    max_file_size = 0,
    in_ext = "",
    out_ext = "",
    in_str = "",
    out_str = "";

export function commonConvertReady(converter) {
    token = sessionStorage.getItem("token");

    // Open Babel uniquely has its own maximum file size
    if (converter == "Open Babel") {
        max_file_size = sessionStorage.getItem("max_file_size_ob");
    } else {
        max_file_size = sessionStorage.getItem("max_file_size");
    }

    // Set the text for displaying the maximum size
    if (max_file_size > 0) {
        $(".max-file-size").text(" (max size " + (max_file_size / MEGABYTE).toFixed(2) + " MB)");
    }

    in_str = sessionStorage.getItem("in_str");
    out_str = sessionStorage.getItem("out_str");

    // $$$$$ FUNCTION FOR THIS? $$$$$
    const in_str_array = in_str.split(": ");
    in_ext = in_str_array[0];                // e.g. "ins"
    const in_note = in_str_array[1];         // e.g. "ShelX"

    const out_str_array = out_str.split(": ");
    out_ext = out_str_array[0];
    const out_note = out_str_array[1];

    $("#heading").html("Convert from \'" + in_ext + "\' (" + in_note + ") to \'" + out_ext + "\' (" + out_note +
        ") using " + converter);

    // Connect the buttons to events
    $("#extCheck").click(setExtCheck);
    $("#requestLog").click(setRequestLog);
    $("#clearUpload").click(clearUploadedFile);

    // Connect the file upload to event and limit the types it can accept
    $("#fileToUpload").change(checkFile);
    limitFileType();

    initDirtyForms();

    return [token, in_str, in_ext, out_str, out_ext];
}

// Converts user-supplied file to another format and downloads the resulting file
export function convertFile(form_data, download_fname, fname) {

    showSpinner();
    disableConvertButton();

    let convertTimedOut = false;

    var jqXHR = $.ajax({
        url: `/convert/`,
        type: "POST",
        data: form_data,
        processData: false,
        contentType: false,
        timeout: CONVERT_TIMEOUT,
        success: async function () {

            hideSpinner();
            enableConvertButton();
            clearUploadedFile();
            disableDirtyForms();

            if (!convertTimedOut) {
                await downloadFile(`../downloads/${download_fname}`, download_fname)

                let msg = "To the best of our knowledge, this conversion has worked. A download prompt for your " +
                    "converted file should now be open for you.";
                if (requestLog) {
                    msg += " A download prompt for the log will appear when you close this box.\n\n" +
                        "You may need to tell your browser to allow this site to download multiple files and try the " +
                        "conversion again if your browser initially disallows the download of the log file.";
                }
                alert(msg);

                if (requestLog) {
                    await downloadFile(`../downloads/${fname}.log.txt`, fname + '.log.txt')
                }
            }

            var fdata = new FormData();

            fdata.append("filename", download_fname);
            fdata.append("logname", fname + '.log.txt');

            $.ajax({
                url: `/delete/`,
                type: "POST",
                data: fdata,
                processData: false,
                contentType: false
            })
                .fail(function (e) {
                    // For debugging
                    console.log("Error deleting remote files after download");
                    console.log(e.status);
                    console.log(e.responseText);
                })
        },
        error: function (xmlhttprequest, textstatus, message) {
            hideSpinner();
            enableConvertButton();
            if (textstatus === "timeout") {
                convertTimedOut = true;
                alert("ERROR: Conversion attempt timed out. This may be because the conversion is too complicated, " +
                    "or because the server is currently busy.");
                console.log("ERROR: Conversion timed out")
            }
        }
    })
        .fail(function (e, textstatus, message) {
            let errLog = `/static/downloads/${fname}.log.txt`;

            fetch(errLog, { cache: "no-store" })
                .then(function (response) {
                    if (response.status == 404) {
                        alert("An unknown error occurred, which produced no error log. If you are using the web " +
                            "app, please provide feedback on the conversion that you were attempting by clicking on " +
                            "'Contact' in the navigation bar.\n" +
                            "If you were trying to run this locally, check your terminal for error messages which " +
                            "may help explain what went wrong.");
                        return "";
                    }
                    else if (!convertTimedOut) {
                        return response.text();
                    }
                })
                .then(function (text) {
                    if (text != "" && text != null)
                        alert("ERROR: Request to the backend converter returned status '" + textstatus +
                            "' and message: " + message + "\n" + text);
                })

            // For debugging
            console.log("Error converting file");
            console.log(e.status);
            console.log(e.responseText);
        })
}

function setExtCheck(event) {
    extCheck = this.checked;

    // Toggle whether or not the file upload limits uploaded type based on whether or not this box is ticked
    if (extCheck) {
        limitFileType();
    } else {
        unlimitFileType();
    }
}

export function getExtCheck() {
    return extCheck;
}

function setRequestLog(event) {
    requestLog = this.checked;
}

export function splitArchiveExt(filename) {
    const filename_segments = filename.split(".");
    // Check for extreme cases
    if (filename_segments.length == 0) {
        return ["", ""];
    } else if (filename_segments.length == 1) {
        return [filename, ""];
    }

    let base;
    let ext = filename_segments.at(-1);

    if ([GZ_EXT, BZ_EXT, XZ_EXT].includes(ext)) {
        // In the case that the extension is one of the second parts of tarball extensions, check if the prior
        // extension is "tar"
        let prior_ext = filename_segments.at(-2);
        if (prior_ext == TAR_EXT && filename_segments.length > 2) {
            base = filename_segments.slice(0, -2).join(".");
            ext = prior_ext + "." + ext;
        } else {
            base = filename_segments.slice(0, -1).join(".");
        }
    } else {
        base = filename_segments.slice(0, -1).join(".");
    }

    return [base, ext];
}

export function isArchiveExt(ext) {
    return [ZIP_EXT, TAR_EXT, TARGZ_EXT, TARBZ_EXT, TARXZ_EXT].includes(ext);
}

// Check that the file meets requirements for upload
function checkFile(event) {

    // Enable dirty form checking whenever a file is uploaded
    enableDirtyForms();

    let allGood = true;
    let file = this.files[0];
    let message = "";

    // Check file has the proper extension if checking is enabled
    if (extCheck) {

        const file_name = file.name;

        const [_, ext] = splitArchiveExt(file_name);

        if (![in_ext, ZIP_EXT, TAR_EXT, TARGZ_EXT, TARBZ_EXT, TARXZ_EXT].includes(ext)) {
            message += "The file extension is not " + in_ext + " or a zip or tar archive extension" +
                ": If you're confident this file is the correct type, untick the box above to disable file extension " +
                "enforcement (note that zip/tar archives MUST have the correct extension, regardless of this " +
                "tickbox). Otherwise, please select another file or change the 'from' format on the 'Home' page.";
            allGood = false;
        }
    }

    // Check file does not exceed maximum size
    if (max_file_size > 0 && file.size > max_file_size) {
        if (message !== "")
            message += "\n\n";
        message += "The file exceeds the maximum size limit of " + (max_file_size / MEGABYTE).toFixed(2) +
            " MB; its size is " + (file.size / MEGABYTE).toFixed(2) + " MB. Please either log in for an increased " +
            "file size limit, or see the Downloads page to run the app locally with no limit.";
        allGood = false;
    }

    if (allGood) {
        enableConvertButton();
    } else {
        disableConvertButton();
        alert(message);
    }
}

/**
 * Allow the file upload to only accept the expected type of file, plus archives
 */
function limitFileType() {
    $("#fileToUpload")[0].accept = "." + in_ext + ", .zip, .tar, .tar.gz, .tar.xz, .tar.bz";
}

/**
 * Allow the file upload to accept any type of file
 */
function unlimitFileType() {
    $("#fileToUpload")[0].accept = "*";
}

/**
 * Clear any uploaded file
 */
function clearUploadedFile() {
    $("#fileToUpload").val('');
    disableConvertButton();
}

/**
 * Enable the "Convert" button
 */
function enableConvertButton() {
    $("#uploadButton").css({ "background-color": "var(--ifm-color-primary)", "color": "var(--ifm-hero-text-color)" });
    $("#uploadButton").prop({ disabled: false });
}


/**
 * Disable the "Convert" button
 */
function disableConvertButton() {
    $("#uploadButton").css({ "background-color": "var(--psdi-bg-color-secondary)", "color": "gray" });
    $("#uploadButton").prop({ disabled: true });
}

/**
 * Start a download of a file
 * 
 * The file is first fetched as a blob, then a link to download it is created, clicked, and removed
 * 
 * @param {str} path The path to the file to be downloaded
 * @param {str} filename The desired filename of the downloaded file
 * @returns {Promise<Response>}
 */
async function downloadFile(path, filename) {
    return fetch(path)
        .then(res => res.blob())
        .then(data => {
            var a = document.createElement("a");
            a.href = window.URL.createObjectURL(data);
            a.download = filename;
            a.click();
            a.remove();
        });
}

/**
 * Show the loading spinner
 */
function showSpinner() {
    $(".loading-spinner").css({ display: "inherit" });
}

/**
 * Hide the loading spinner
 */
function hideSpinner() {
    $(".loading-spinner").css({ display: "none" });
}