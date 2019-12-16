"""Running seeds through WWWBLAST."""

from multiprocessing import Pool
import os
import textwrap

from Bio.Blast import NCBIWWW
from logzero import logger

from .ref_download import load_tsv
from .cli import _proc_args


def do_blast(x):
    record, args = x

    out_dir = os.path.dirname(args.in_tsv)
    out_file = os.path.basename(record["path"])[: -len(".fasta")] + ".blast.xml"
    out_path = os.path.join(out_dir, out_file)

    logger.info("Loading %s", record["path"])
    with open(os.path.join(out_dir, record["path"]), "rt") as inputf:
        fasta = inputf.read()

    logger.info("Running NCBI WWW BLAST for %s", fasta.split()[0][1:])
    blast_res = NCBIWWW.qblast(program="blastn", database="nt", sequence=fasta, megablast=True)
    with open(out_path, "wt") as outputf:
        logger.info("Writing NCBI WWW BLAST result for %s", fasta.split()[0][1:])
        outputf.write(blast_res.read())


def run(parser, args):
    """Perform NCBIWWW of the seeds."""
    logger.info("Starting NCBI WWW BLAST of seed sequences and reference.")
    logger.info("Arguments are %s", vars(args))

    out_dir = os.path.dirname(args.in_tsv)

    jobs = []

    logger.info("Running NCBI WWW BLAST for each sequence...")
    header, records = load_tsv(args.in_tsv)
    for record in records:
        out_file = os.path.basename(record["path"])[: -len(".fasta")] + ".blast.xml"
        out_path = os.path.join(out_dir, out_file)
        if os.path.exists(out_path):
            logger.warn("Output path already exists: %s", out_path)
            continue
        else:
            jobs.append((record, args))

    with Pool(args.num_threads) as p:
        p.map(do_blast, jobs, chunksize=1)

    logger.info("All done. Have a nice day!")


def add_parser(subparser):
    """Configure the ``argparse`` sub parser."""
    parser = subparser.add_parser("ref_blast")
    parser.set_defaults(func=run)

    parser.add_argument(
        "--num-threads", default=8, type=int, help="Number of parallel searches to run"
    )
    parser.add_argument("in_tsv", help="Path to input TSV file.")
