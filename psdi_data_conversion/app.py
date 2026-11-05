"""app.py

Entry-point module for when Flask is called directly to start the server
"""

from psdi_data_conversion.gui.setup import get_app_wrapper

app = get_app_wrapper()
