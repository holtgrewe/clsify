"""Helpers for converting read files."""

import os

from bioconvert.scf2fastq import SCF2FASTQ
from bioconvert.abi2fastq import ABI2FASTQ
from logzero import logger


def convert_seqs(seq_files, tmpdir):
    """Convert SRF and AB1 files to FASTQ."""
    logger.info("Running file conversion...")
    result = []
    for seq_file in seq_files:
        path_fastq = os.path.join(tmpdir, os.path.basename(seq_file)[:-4])
        if seq_file.endswith(".scf"):
            logger.info("Converting SCF file %s...", seq_file)
            SCF2FASTQ(seq_file, path_fastq)()
            result.append(path_fastq)
        elif seq_file.endswith(".ab1"):
            logger.info("Converting ABI file %s...", seq_file)
            ABI2FASTQ(seq_file, path_fastq)()
            result.append(path_fastq)
        else:
            result.append(seq_file)
    logger.info("Done converting files.")
    return result
