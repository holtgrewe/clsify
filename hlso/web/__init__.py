"""Code for the web interface to Haplotype-Lso."""

import os

from logzero import logger

from . import settings


def run(parser, args):
    """Run the ``hlso`` command line interface."""
    logger.info("Configuring settings from arguments %s", args)
    settings.HOST = args.host
    settings.PORT = args.port
    settings.PUBLIC_URL_PREFIX = args.public_url_prefix

    logger.info("Running server...")
    from .app import app  # noqa

    app.run_server(
        host=settings.HOST, port=settings.PORT, debug=args.debug, dev_tools_hot_reload=args.debug
    )
    logger.info("Web server stopped. Have a nice day!")


def add_parser(subparser):
    """Configure the ``argparse`` sub parser."""
    parser = subparser.add_parser("web")
    parser.set_defaults(func=run)

    parser.add_argument(
        "--host", help="Server host", default=os.environ.get("HLSO_HOST", "0.0.0.0")
    )
    parser.add_argument(
        "--port", help="Server port", type=int, default=int(os.environ.get("HLSO_PORT", "8050"))
    )
    parser.add_argument(
        "--public-url-prefix",
        default=os.environ.get("HLSO_URL_PREFIX", ""),
        help="The prefix that this app will be served under (e.g., if behind a reverse proxy.)",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        default=(os.environ.get("DEBUG", "0").lower() in ("1", "true", "yes")),
        help="Whether or not to enable debugging.",
    )
