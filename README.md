# CHEMISTRY FILE FORMAT CONVERSION DATABASE

This is the repository for the Pathfinder 2 Chemistry File Format Conversion source code.


## Directory structure of the Python Flask app website:

app.py
requirements.txt
  (azure-functions, requests, openbabel-wheel:
   installed by yaml workflow file on GitHub)
startup.sh*
static
    downloads
        dummy
          (a blank file)
    images
        favicon.ico
        PSDI_Logo_CMYK_282c.svg
    javascript
        byte.js
          (byte array representation of SQLite database)
        convert.py
          (converts format.db to byte.js - not automatic)
        format.db
          (native representation of SQLite database)
        format.js
        node_modules
          (containing SQL.js)
    styles
        format.css
    uploads
        dummy
          (a blank file)
templates
    index.htm

*Must be copied to another location on the Azure site. Once there, it persists. Installs libraries missing from the Docker container (required by Open Babel). Go to     Development Tools > SSH > Go    to arrive at a terminal for the site. At    /home    place file startup.sh. Go to    Settings > Configurations > General settings    and enter    /home/startup.sh    in Startup Command. Note that startup.sh can be placed anywhere below    /home    (provided the paths match).


## Running the Python Flask app hosted on the Microsoft Azure site

Enter    https://psdidev2.azurewebsites.net    in a browser.


## Running the Python/Flask app locally

Python and Open Babel must be installed.

From the command line (in the appropriate directory), enter:
export FLASK_APP=app.py

To enable debug mode, if required, enter:
export FLASK_ENV=development

Run server app.py by entering    Python3.9 -m flask run    (the Python version may differ) at the command line (in the appropriate directory).

Run the website by entering    127.0.0.1:5000    in a browser.

In case of problems when using Chrome, try opening Chrome from the command line:
open -a "Google Chrome.app" --args --allow-file-access-from-files


## Using the website

Typing in the text boxes at the top of 'Convert from' and Convert to' filters the lists in the lower text areas. Selecting from both of these text areas populates the ‘Conversion success’ dropdown list appropriately. Selecting from the latter causes converter details to be displayed, including a link to its website. Currently, because the database is incomplete, the 'Conversion success' box may remain empty.

If there are no converters listed, or if a file format is not found, the user is able to provide feedback. The server writes the feedback to a server-side file. Note that this does not make complete sense at the moment because of the incomplete database.

If a conversion is supported by Open Babel, such a conversion is offered. If the user clicks on the 'Yes' radio button, read (input) and/or write (output) option flags are listed (or not) depending on the selected formats. Options may be selected if required, and the user uploads a file with the appropriate extension. A file in the new format downloads automatically.

The user can send an email to psdi@soton.ac.uk (top of screen), and there is a button at the bottom for showing/hiding notes.

## Database

If changed, the SQLite database format.db must be converted to a byte array so that SQL.js can read it and thus register the changes. The conversion is effected by entering     python3.9 convert.py    at the command line (with format.db and convert.py in the same directory - the Python version may differ). This produces the file byte.js, which is used by module format.js. It is intended that PostgreSQL will replace SQLite in due course.

The website user is not able to update the database, but instead provides feedback.
