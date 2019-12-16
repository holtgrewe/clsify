"""Download of seed references."""

import os
import shutil
import tempfile

from logzero import logger

from .cli import _proc_args
from .conversion import convert_seqs


def run(parser, args):
    """Perform conversion."""
    args = _proc_args(parser, args)
    logger.info("Converting sequences...")
    logger.info("Args = %s", args)
    with tempfile.TemporaryDirectory() as tmpdir:
        seq_files = convert_seqs(args.seq_files, tmpdir)
        for seq_file in seq_files:
            shutil.copy(seq_file, os.path.join(args.out_dir, os.path.basename(seq_file)))
    logger.info("All done, have a nice day!")


def add_parser(subparser):
    """Configure the ``argparse`` sub parser."""
    parser = subparser.add_parser("convert")
    parser.set_defaults(func=run)

    parser.add_argument(
        "--file-name-as-seq-name",
        action="store_true",
        default=False,
        help="Set file name to sample name",
    )
    parser.add_argument("out_dir", help="Path to output directory.")
    parser.add_argument("seq_files", nargs="+", default=[], action="append")
