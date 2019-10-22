"""Code for the command line interface to ``clsify``."""

import tempfile

from logzero import logger

from .conversion import convert_seqs
from .export import write_excel
from .workflow import blast_and_haplotype_many, result_pairs_to_data_frames

# from .analysis import convert_seqs, infer_from_file, write_ref
# from .blast import run_blastn
# from .haplo import HAPLOTYPE_NAMES, HaplotypeMatchScore, run_haplotyping, PlusMinus
from .web.settings import SAMPLE_REGEX


def _proc_args(parser, args):
    seq_files = []
    for lst in args.seq_files:
        seq_files += [x for y in lst for x in y.split(" ") if x]
    args.seq_files = seq_files
    return args


def run(parser, args):
    """Run the ``clsify`` command line interface."""
    args = _proc_args(parser, args)
    logger.info("Starting Lso classification.")
    logger.info("Arguments are %s", vars(args))
    with tempfile.TemporaryDirectory() as tmpdir:
        logger.info("Converting sequences (if necessary)...")
        seq_files = convert_seqs(args.seq_files, tmpdir)
        logger.info("Running BLAST and haplotyping...")
        result_pairs = blast_and_haplotype_many(seq_files)
        logger.info("Converting results into data frames...")
        df_summary, df_blast, df_haplotyping = result_pairs_to_data_frames(
            result_pairs, args.sample_regex
        )
        logger.info("Writing XLSX file to %s", args.output)
        write_excel(df_summary, df_blast, df_haplotyping, args.output)
    logger.info("All done. Have a nice day!")


def add_parser(subparser):
    """Configure the ``argparse`` sub parser."""
    parser = subparser.add_parser("cli")
    parser.set_defaults(func=run)

    parser.add_argument(
        "--sample-regex",
        default=SAMPLE_REGEX,
        help="Regular expression to match file name to sample name.",
    )
    parser.add_argument("-o", "--output", default="clsified.xlsx", help="Path to output file")
    parser.add_argument("seq_files", nargs="+", default=[], action="append")
