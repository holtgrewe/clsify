"""Download of seed references."""

import os
import textwrap

from Bio import Entrez
from logzero import logger

from .cli import _proc_args
from .common import load_tsv

#: Known FASTA extensions.
FASTA_EXTS = (".fasta", ".fa", ".fna")

#: Reference sequences.
REF_SEQS = {
    "16S": {"accession": "EU812559.1", "start": 0, "end": 1225},
    "16S-23S": {"accession": "EU812559.1", "start": 1225, "end": 2515},
    "50S": {"accession": "EU834131.1", "start": 0, "end": 1714},
}


def download_references(_parser, args):
    """Download references."""
    out_dir = os.path.dirname(args.out_tsv)

    logger.info("Downloading references...")
    for name, vals in REF_SEQS.items():
        out_path_fasta = os.path.join(out_dir, "ref_" + name + ".fasta")
        if os.path.exists(out_path_fasta):
            logger.warn("- already exists: %s", vals["accession"])
        else:
            logger.info("- downloading: %s", vals["accession"])
            res = Entrez.efetch(
                db="nucleotide", id=vals["accession"], rettype="fasta", retmode="text"
            )
            header, seq = res.read().rstrip().split("\n", 1)
            seq = "".join([c for c in seq if c in "cgatnCGATN"])
            prefix = "N" * vals["start"]
            suffix = "N" * (len(seq) - vals["end"])
            seq = prefix + seq[vals["start"] : vals["end"]] + suffix
            with open(out_path_fasta, "wt") as outputf:
                print(header, file=outputf)
                print("\n".join(textwrap.wrap(seq)), file=outputf)


def download_seeds(parser, args):
    """Download seed sequences"""
    out_dir = os.path.dirname(args.out_tsv)
    out_path_tsv = os.path.join(out_dir, "seeds_paths.tsv")
    if os.path.exists(out_path_tsv):
        logger.warn("Output file already exists: %s", out_path_tsv)
        logger.warn("Skipping reference download.")
        return

    logger.info("Loading input TSV...")
    header, records = load_tsv(args.in_tsv)

    logger.info("Downloading sequences...")
    for record in records:
        if os.path.exists(os.path.join(out_dir, record["accession"])):
            logger.info("- already exists: %s", os.path.join(out_dir, record["accession"]))
            record["path"] = record["accession"]
            out_path = os.path.join(out_dir, record["accession"])
            if not os.path.exists(out_path):
                parser.exit(1, "File not found: %s" % out_path)
        else:
            logger.info("- downloading: %s", record["accession"])
            record["path"] = record["accession"] + ".fasta"
            out_path_fasta = os.path.join(out_dir, record["path"])
            res = Entrez.efetch(
                db="nucleotide", id=record["accession"], rettype="fasta", retmode="text"
            )
            record = res.read()
            with open(out_path_fasta, "wt") as outputf:
                outputf.write(record.rstrip() + "\n")

    logger.info("Writing output TSV: %s", out_path_tsv)
    header[-1] = "path"
    with open(out_path_tsv, "wt") as outputf:
        print("\t".join(header), file=outputf)
        for record in records:
            print("\t".join([record[k] for k in header]), file=outputf)


def run(parser, args):
    """Perform the download of seed sequences and references."""
    logger.info("Starting download of seed sequences and reference.")
    logger.info("Arguments are %s", vars(args))

    download_references(parser, args)
    download_seeds(parser, args)

    logger.info("All done. Have a nice day!")


def add_parser(subparser):
    """Configure the ``argparse`` sub parser."""
    parser = subparser.add_parser("ref_download")
    parser.set_defaults(func=run)

    parser.add_argument("in_tsv", help="Path to output TSV file.")
    parser.add_argument("out_tsv", help="Path to output TSV file.")
