"""Helpers for converting read files."""

import os
import tempfile
import typing

from bioconvert.scf2fasta import SCF2FASTA
from bioconvert.abi2fasta import ABI2FASTA
from bioconvert.fastq2fasta import FASTQ2FASTA
from logzero import logger

from .common import load_fasta


def convert_seqs(
    seq_files: typing.Iterable[str], tmpdir: str, sample_name_from_file_name: bool = False
) -> typing.List[str]:
    """Convert SRF and AB1 files to FASTQ."""
    logger.info("Running file conversion...")
    result = []

    for seq_file in seq_files:
        file_basename = os.path.basename(seq_file)[: -len(".fasta")]
        path_fasta = os.path.join(tmpdir, file_basename) + ".fasta"
        result.append(path_fasta)

        path_fasta_tmp = tempfile.NamedTemporaryFile(dir=tmpdir, delete=False)
        if seq_file.endswith(".scf"):
            logger.info("Converting SCF file %s...", seq_file)
            SCF2FASTA(seq_file, path_fasta_tmp)()
            result.append(path_fasta_tmp)
        elif seq_file.endswith(".ab1"):
            logger.info("Converting ABI file %s...", seq_file)
            ABI2FASTA(seq_file, path_fasta_tmp)()
            result.append(path_fasta_tmp)
        elif seq_file.endswith(".fastq"):
            logger.info("Converting FASTAQ file %s...", seq_file)
            FASTQ2FASTA(seq_file, path_fasta_tmp)()
            result.append(path_fasta_tmp)
        else:
            path_fasta_tmp = seq_file

        fasta_content = load_fasta(path_fasta_tmp)

        if sample_name_from_file_name and len(fasta_content) != 1:
            prefix_no = 1
        else:
            prefix_no = 0

        with open(path_fasta, "wt") as outputf:
            for name, seq in fasta_content.items():
                if sample_name_from_file_name:
                    name = file_basename
                if prefix_no:
                    prefix = "_%d" % prefix_no
                    prefix_no += 1
                else:
                    prefix = ""
                print(">%s%s\n%s" % (prefix, name, seq), file=outputf)

    logger.info("Done converting files.")
    return result
