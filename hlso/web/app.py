"""Setup of Haplotype-Lso Dash application."""

import os

import dash
import flask

from . import callbacks, settings
from .ui import build_layout
from hlso import __version__

#: Path to assets.
ASSETS_FOLDER = os.path.join(os.path.dirname(__file__), "assets")

#: The Flask application to use.
app_flask = flask.Flask(__name__)

# Setup URL prefix for Flask.
app_flask.config["APPLICATION_ROOT"] = "%s/" % settings.PUBLIC_URL_PREFIX

#: The Dash application to run.
app = dash.Dash(
    __name__,
    # Use our specific Flask app
    server=app_flask,
    # All files from "assets/" will be served as "/assets/*"
    assets_folder=ASSETS_FOLDER,
    # The visualization will be served below "/dash"
    routes_pathname_prefix="/dash/",
    requests_pathname_prefix="%s/dash/" % settings.PUBLIC_URL_PREFIX,
)

# Set app title
app.title = "hlso v%s" % __version__

# Initialize UI
app.layout = build_layout()

# Serve assets locally
app.css.config.serve_locally = True
app.scripts.config.serve_locally = True

# TODO: register callbacks
callbacks.register_upload(app)
callbacks.register_computation_complete(app)
callbacks.register_row_clicks(app)

# Add redirection for root.
@app_flask.route("/")
def redirect_root():
    return flask.redirect("%s/dash/" % settings.PUBLIC_URL_PREFIX)
