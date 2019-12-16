"""Code for the command line interface to ``hlso``."""

import tempfile
import typing

import attr
from logzero import logger

from .conversion import convert_seqs
from .export import write_excel
from .workflow import blast_and_haplotype_many, results_to_data_frames

# from .analysis import convert_seqs, infer_from_file, write_ref
# from .blast import run_blastn
# from .haplo import HAPLOTYPE_NAMES, HaplotypeMatchScore, run_haplotyping, PlusMinus
from .web.settings import SAMPLE_REGEX


def _proc_args(parser, args):
    seq_files = []
    for lst in args.seq_files:
        seq_files += lst
    args.seq_files = seq_files
    return args


@attr.s(auto_attribs=True, frozen=True)
class Config:
    """Configuration of the command line interface."""

    #: The paths to the input files.
    input_paths: typing.Tuple[str]
    #: The path to the output file.
    output_path: str
    #: Whether or not the sample names are to be inferred from the command line.
    #: If not, then the sequence names are used.
    sample_name_from_file: bool = False
    #: The regular expression to parse information from the sample.
    sample_regex: str = SAMPLE_REGEX


def run(parser, args):
    """Run the ``hlso`` command line interface."""
    args = _proc_args(parser, args)
    config = Config(
        input_paths=tuple(args.seq_files),
        output_path=args.output,
        sample_name_from_file=args.sample_name_from_file,
        sample_regex=args.sample_regex,
    )
    logger.info("Starting Lso classification.")
    logger.info("Arguments are %s", config)
    with tempfile.TemporaryDirectory() as tmpdir:
        logger.info("Converting sequences (if necessary)...")
        seq_files = convert_seqs(args.seq_files, tmpdir, config.sample_name_from_file)
        config = Config(**{**attr.asdict(config), "input_paths": tuple(sorted(seq_files))})
        logger.info("Running BLAST and haplotyping...")
        results = blast_and_haplotype_many(seq_files)
        logger.info("Converting results into data frames...")
        df_summary, df_blast, df_haplotyping = results_to_data_frames(results, args.sample_regex)
        logger.info("Summary:\n%s", df_summary)
        logger.info("Writing XLSX file to %s", args.output)
        write_excel(df_summary, df_blast, df_haplotyping, args.output)
    logger.info("All done. Have a nice day!")


def add_parser(subparser):
    """Configure the ``argparse`` sub parser."""
    parser = subparser.add_parser("cli")
    parser.set_defaults(func=run)

    parser.add_argument(
        "--sample-name-from-file",
        action="store_true",
        default=False,
        help="Use sample name instead of file name",
    )
    parser.add_argument(
        "--sample-regex",
        default=SAMPLE_REGEX,
        help="Regular expression to match file name to sample name.",
    )
    parser.add_argument("-o", "--output", default="clsified.xlsx", help="Path to output file")
    parser.add_argument("seq_files", nargs="+", default=[], action="append")
