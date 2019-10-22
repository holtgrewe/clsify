"""Classification of C. liberibacter solanacearum following IPPC standards.

This module is the main entry point.
"""

import argparse
from subprocess import check_output, CalledProcessError
import sys

from . import cli, web


def main(argv=None):
    """Main entrypoint (before parsing command line arguments)."""
    parser = argparse.ArgumentParser(description="Classify Lso Sanger reads.")
    subparsers = parser.add_subparsers()
    cli.add_parser(subparsers)
    web.add_parser(subparsers)

    for prog in ("blastn", "bcftools"):
        found_all = True
        try:
            check_output(["which", prog])
        except CalledProcessError:
            found_all = False
            print("ERROR: Required program %s not found!" % repr(prog), file=sys.stderr)
    if not found_all:
        parser.exit(1)

    args = parser.parse_args(argv)
    if not hasattr(args, "func"):
        parser.print_help()
        parser.exit(1)
    else:
        return args.func(parser, args)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
