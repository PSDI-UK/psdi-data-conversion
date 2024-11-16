# CHEMISTRY FILE FORMAT CONVERSION DATABASE

This is the repository for the Pathfinder 2 Chemistry File Format Conversion source code.


## Directory structure of the Python Flask app website:

app.py
requirements.txt
  (azure-functions, requests, openbabel-wheel:
   installed by yaml workflow file on GitHub)
startup.sh*
static
    content
        accessibility.htm
        convert.htm
        convertato.htm
        documentation.htm
        feedback.htm
        header-links.htm
        index-header-links.htm
        report.htm
    downloads (created by app.py if not extant)
    images
        PSDI_Logo_CMYK_282c.svg
    javascript
        convert.js
        convert.py
        convertato.js
        format.js
        load_accessibility.js
        report.js
    styles
        format.css
    uploads (created by app.py if not extant)
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

The database can only be accessed from the University of Southampton or when using Global Connect. The current IP address must be added to the database's firewall rules on the Azure site.

In case of problems when using Chrome, try opening Chrome from the command line:
open -a "Google Chrome.app" --args --allow-file-access-from-files


## Using the website

Guidance on usage is given on each page of the website.


## Database

A PostgreSQL database is hosted on the same Azure site. The website user is not able to update the database, but instead provides feedback.
