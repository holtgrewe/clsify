"""Prepare modified copies of reference sequence by inserting BLAST matches into them."""

import os
import tempfile
import textwrap

from logzero import logger

from .common import load_fasta, revcomp
from .conversion import convert_seqs
from .workflow import blast_and_haplotype_many, REF_FILE
from .cli import _proc_args


def do_paste(match, ref_seqs=None):
    ref_seqs = ref_seqs or load_fasta(REF_FILE)
    seq = ref_seqs[match.database]
    qseq = "".join([c for c in match.alignment.qseq if c in "cgatnCGATN"])
    if match.query_strand == "-":
        qseq = revcomp(seq)
    seq = seq[: match.database_start] + qseq + seq[match.database_end :]
    return "".join([c for c in seq if c != "N"]).upper()


def write_pasted(results, output_prefix):
    logger.info("Loading reference sequences...")
    ref_seqs = load_fasta(REF_FILE)
    logger.info("Writing pasted sequences...")
    for result in results:
        for match in result["matches"]:
            if match.database is None:
                logger.info("  => no match for %s", match.query)
                continue
            query = match.query.split(" ")[0]
            out_path = "{}{}/{}-{}.fasta".format(
                output_prefix, match.database, query, match.database
            )
            if not os.path.exists(os.path.dirname(out_path)):
                os.makedirs(os.path.dirname(out_path))
            logger.info("  - %s", out_path)
            with open(out_path, "wt") as outputf:
                seq = do_paste(match, ref_seqs)
                # print(">{}-{}".format(query, match.database), file=outputf)
                print(">{} ({} bp)".format(query, len(seq)), file=outputf)
                print("\n".join(textwrap.wrap(seq)), file=outputf)


def run(parser, args):
    """Run the ``hlso`` command line interface."""
    args = _proc_args(parser, args)
    logger.info("Starting Lso classification.")
    logger.info("Arguments are %s", vars(args))
    with tempfile.TemporaryDirectory() as tmpdir:
        logger.info("Converting sequences (if necessary)...")
        seq_files = convert_seqs(args.seq_files, tmpdir)
        logger.info("Running BLAST (and haplotyping)...")
        results = blast_and_haplotype_many(seq_files)
        logger.info("Writing out pasted sequence")
        write_pasted(results, args.output_prefix)
    logger.info("All done. Have a nice day!")


def add_parser(subparser):
    """Configure the ``argparse`` sub parser."""
    parser = subparser.add_parser("paste")
    parser.set_defaults(func=run)

    parser.add_argument(
        "--keep-masked",
        default=False,
        help="Keep masked reference sequence (mostly useful for debugging purposes).",
    )
    parser.add_argument(
        "-o", "--output-prefix", default="hlso_paste_out.d/", help="Prefix for output files"
    )
    parser.add_argument("seq_files", nargs="+", default=[], action="append")
