"""Helpers for converting read files."""

import os

from bioconvert.scf2fasta import SCF2FASTA
from bioconvert.abi2fasta import ABI2FASTA
from logzero import logger


def convert_seqs(seq_files, tmpdir):
    """Convert SRF and AB1 files to FASTQ."""
    logger.info("Running file conversion...")
    result = []
    for seq_file in seq_files:
        path_fasta = os.path.join(tmpdir, os.path.basename(seq_file)[:-4]) + ".fasta"
        if seq_file.endswith(".scf"):
            logger.info("Converting SCF file %s...", seq_file)
            SCF2FASTA(seq_file, path_fasta)()
            result.append(path_fasta)
        elif seq_file.endswith(".ab1"):
            logger.info("Converting ABI file %s...", seq_file)
            ABI2FASTA(seq_file, path_fasta)()
            result.append(path_fasta)
        else:
            result.append(seq_file)
    logger.info("Done converting files.")
    return result
