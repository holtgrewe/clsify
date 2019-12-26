"""Code for running blastn and handling matches.

We store matches in ``BlastMatch`` objects that can be easily converted into dicts and
subsequently Pandas dataframe rows.
"""

import os
import subprocess
import typing
import shlex
import json

import attr
from Bio.Blast import NCBIXML
import io
from logzero import logger

from .common import revcomp, rev


@attr.s(auto_attribs=True, frozen=True)
class Alignment:
    """Representation of an alignment."""

    #: (database) hit sequence
    hseq: str
    #: alignment mid line
    midline: str
    #: query sequence
    qseq: str

    def revcomp(self) -> typing.TypeVar("Alignment"):
        return Alignment(revcomp(self.hseq), "".join(reversed(self.midline)), revcomp(self.qseq))

    def wrapped(self, qry_start: int, db_start: int, line_length: int = 60) -> str:
        result = []
        for offset in range(0, len(self.hseq), line_length):
            end = min(len(self.hseq), offset + line_length)
            result += [
                ("Sbjct %%4d %%-%ds %%s" % line_length)
                % (db_start + offset + 1, self.hseq[offset:end], db_start + offset + line_length),
                ("      %%4s %%-%ds" % line_length) % ("", "".join(self.midline[offset:end])),
                ("Query %%4d %%-%ds %%s" % line_length)
                % (qry_start + offset + 1, self.qseq[offset:end], qry_start + offset + line_length),
                "",
            ]
        return "\n".join(result)

    @staticmethod
    def build_empty() -> typing.TypeVar("Alignment"):
        return Alignment(hseq="", midline="", qseq="")


@attr.s(auto_attribs=True, frozen=True)
class BlastMatch:
    """Representation of a match."""

    #: query file name
    path: typing.Optional[str]
    #: query sequence name
    query: str
    #: database sequence name
    database: str
    #: identity fraction
    identity: float
    #: bit score
    bits: float
    #: strand of query ("+" or "-")
    query_strand: str
    #: 0-based start position
    query_start: int
    #: 0-based end position
    query_end: int
    #: strand of database ("+" or "-")
    database_strand: str
    #: 0-based start position
    database_start: int
    #: 0-based end position
    database_end: int
    #: CIGAR string of match
    match_cigar: str
    #: matching sequence
    match_seq: str
    #: Alignment
    alignment: typing.Optional[Alignment]

    @property
    def database_length(self):
        return self.database_end - self.database_start

    @property
    def query_length(self):
        return self.query_end - self.query_start

    @property
    def is_match(self) -> bool:
        return bool(self.database)

    @staticmethod
    def build_nomatch(query, path: str = None) -> typing.TypeVar("BlastMatch"):
        return BlastMatch(
            path=path,
            query=query,
            database=None,
            identity=0.0,
            bits=0.0,
            query_strand=".",
            query_start=0,
            query_end=0,
            database_strand=".",
            database_start=0,
            database_end=0,
            match_cigar="",
            match_seq="",
            alignment=Alignment.build_empty(),
        )


def is_nucl(c: str) -> bool:
    """Helper, returns wether ``c`` is a nucleotide character."""
    return c in "acgtnACGTN"


def match_cigar(
    qseq: str, hseq: str, query_start: int, query_end: int, query_len: int
) -> typing.Tuple[typing.Tuple[int, str]]:
    """Return CIGAR string for the match with the given paraemters."""
    match_cigar = []
    if query_start > 0:
        match_cigar.append([query_start, "H"])
    for q, h in zip(qseq, hseq):
        if q == "-":
            op = "D"
        elif h == "-":
            op = "I"
        else:
            op = "M"
        if match_cigar and match_cigar[-1][1] == op:
            match_cigar[-1][0] += 1
        else:
            match_cigar.append([1, op])
    if query_end != query_len:
        match_cigar.append([query_len - query_end, "H"])
    return tuple(tuple(t) for t in match_cigar)


def only_dna(seq: str) -> str:
    """Return ``str`` with all ACGTN characters from ``seq``."""
    return "".join(filter(lambda x: x in "ACGTNacgtn", seq))


def parse_blastn_xml(
    blastn_xml: str, path_query: typing.Optional[str] = None
) -> typing.Tuple[BlastMatch]:
    """Parse BLASTN output in JSON format and return list of ``BlastMatch``."""
    result = []
    for blast_record in NCBIXML.parse(io.StringIO(blastn_xml)):
        for blast_alignment in blast_record.alignments:
            hsp = blast_alignment.hsps[0]
            qseq = only_dna(hsp.query)
            hseq = only_dna(hsp.sbjct)
            db_strand = "+" if hsp.strand[1] == "Plus" else "-"
            db_start = min(hsp.sbjct_start, hsp.sbjct_end) - 1
            db_end = max(hsp.sbjct_start, hsp.sbjct_end)
            query_strand = "+" if hsp.strand[0] == "Plus" else "-"
            query_start = min(hsp.query_start, hsp.query_end) - 1
            query_end = max(hsp.query_start, hsp.query_end)
            if db_strand == "-":
                db_strand = "+"
                query_strand = "-" if query_strand == "+" else "+"
                qseq = revcomp(qseq)
                hseq = revcomp(hseq)
                alignment = Alignment(
                    hseq=revcomp(hsp.sbjct), midline=rev(hsp.match), qseq=revcomp(hsp.query)
                )
            else:
                alignment = Alignment(hseq=hsp.sbjct, midline=hsp.match, qseq=hsp.query)
            cigar = match_cigar(qseq, hseq, query_start, query_end, blast_record.query_letters)
            result.append(
                BlastMatch(
                    path=path_query,
                    query=blast_record.query,
                    database=blast_alignment.title.split()[1],
                    bits=hsp.bits,
                    identity=hsp.identities / hsp.align_length,
                    query_strand=query_strand,
                    query_start=query_start,
                    query_end=query_end,
                    database_strand=db_strand,
                    database_start=db_start,
                    database_end=db_end,
                    match_cigar="".join(["".join(map(str, x)) for x in cigar]),
                    match_seq="".join([x for x in filter(is_nucl, qseq)]),
                    alignment=alignment,
                )
            )
    return tuple(result)


def run_blast(database: str, query: str) -> typing.Tuple[BlastMatch]:
    """Run blastn on FASTA query ``query`` to database sequence at ``database``."""
    cmd = ("blastn", "-db", database, "-query", query, "-outfmt", "16")
    logger.info("Executing %s", repr(" ".join(cmd)))
    return parse_blastn_xml(subprocess.check_output(cmd).decode("utf-8"), path_query=query)


def run_makeblastdb(path: str, dbtype: str = "nucl") -> str:
    """Run blastn on FASTA query ``query`` to database sequence at ``database``."""
    cmd = ("makeblastdb", "-in", path, "-dbtype", "nucl")
    logger.info("Executing %s", repr(" ".join(cmd)))
    return subprocess.check_output(cmd).decode("utf-8")
