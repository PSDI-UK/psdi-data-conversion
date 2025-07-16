"""app.py

Entry-point module for when Flask is called directly to start the server

TEST CHANGE
"""

from psdi_data_conversion.gui.setup import get_app

app = get_app()
