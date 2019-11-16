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
from logzero import logger


def revcomp(seq):
    m = {"a": "t", "A": "T", "t": "a", "T": "A", "c": "g", "C": "G", "g": "c", "G": "C"}
    return "".join(reversed(list(map(lambda x: m.get(x, x), seq))))


@attr.s(auto_attribs=True, frozen=True)
class Alignment:
    """Representation of an alignment."""

    #: (database) hit sequence
    hseq: str
    #: alignment mid line
    midline: str
    #: query sequence
    qseq: str

    def revcomp(self):
        return Alignment(revcomp(self.hseq), "".join(reversed(self.midline)), revcomp(self.qseq))

    def wrapped(self, qry_start, db_start, line_length=60):
        result = []
        for offset in range(0, len(self.hseq), line_length):
            end = min(len(self.hseq), offset + line_length)
            result += [
                ("Sbjct %%4d %%-%ds %%s" % line_length)
                % (db_start + offset + 1, self.hseq[offset:end], db_start + offset + line_length),
                ("      %%4s %%-%ds" % line_length) % ("", self.midline[offset:end]),
                ("Query %%4d %%-%ds %%s" % line_length)
                % (qry_start + offset + 1, self.qseq[offset:end], qry_start + offset + line_length),
                "",
            ]
        return "\n".join(result)

    @staticmethod
    def build_empty():
        return Alignment(hseq="", midline="", qseq="")


@attr.s(auto_attribs=True, frozen=True)
class BlastMatch:
    """Representation of a match."""

    #: query sequence name
    query: str
    #: database sequence name
    database: str
    #: identity fraction
    identity: float
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
    #: string with original JSON of match
    match_json: str

    @property
    def is_match(self):
        return bool(self.database)

    @staticmethod
    def build_nomatch(query, match_json):
        return BlastMatch(
            query=query,
            database=None,
            identity=0.0,
            query_strand=".",
            query_start=0,
            query_end=0,
            database_strand=".",
            database_start=0,
            database_end=0,
            match_cigar="",
            match_seq="",
            alignment=Alignment.build_empty(),
            match_json=match_json,
        )


def is_nucl(c):
    """Helper, returns wether ``c`` is a nucleotide character."""
    return c in "acgtnACGTN"


def match_cigar(qseq, hseq, query_start, query_end, query_len):
    """Return CIGAR string for the match with the given paraemters."""
    match_cigar = []
    if query_start > 0:
        match_cigar.append([query_start, "H"])
    for i, (q, h) in enumerate(zip(qseq, hseq)):
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
    return match_cigar


def parse_blastn_json(blastn_json):
    """Parse BLASTN output in JSON format and return ``BlastMatch``."""
    match_dict = json.loads(blastn_json)
    search = match_dict["BlastOutput2"][0]["report"]["results"]["search"]
    hits = search["hits"]
    result = []
    for hit in hits:
        hsp = hit["hsps"][0]
        qseq = "".join(filter(lambda x: x in "ACGTNacgtn", hsp["qseq"]))
        hseq = "".join(filter(lambda x: x in "ACGTNacgtn", hsp["hseq"]))
        db_strand = "+" if hsp["hit_strand"] == "Plus" else "-"
        db_start = min(hsp["hit_from"], hsp["hit_to"]) - 1
        query_start = min(hsp["query_from"], hsp["query_to"]) - 1
        query_end = max(hsp["query_from"], hsp["query_to"])
        query_strand = "+" if hsp["query_strand"] == "Plus" else "-"
        alignment = Alignment(hseq=hsp["hseq"], midline=hsp["midline"], qseq=hsp["qseq"])
        if db_strand != "+":  # ensure database is always forward
            logger.info("Reverse-complementing match!")
            db_strand = "+"
            query_strand = {"+": "-", "-": "+"}[query_strand]
            qseq = revcomp(qseq)
            hseq = revcomp(hseq)
            alignment = alignment.revcomp()
        cigar = match_cigar(qseq, hseq, query_start, query_end, search["query_len"])
        result.append(
            BlastMatch(
                query=search["query_title"],
                database=hit["description"][0]["title"].split()[0],
                identity=hsp["identity"] / hsp["align_len"],
                query_strand=query_strand,
                query_start=query_start,
                query_end=query_end,
                database_strand=db_strand,
                database_start=db_start,
                database_end=max(hsp["hit_from"], hsp["hit_to"]),
                match_cigar="".join(["".join(map(str, x)) for x in cigar]),
                match_seq="".join([x for x in filter(is_nucl, qseq)]),
                alignment=alignment,
                match_json=blastn_json,
            )
        )
    return result


def run_blast(database, query):
    """Run blastn on FASTA query ``query`` to database sequence at ``database``."""
    cmd = ("blastn", "-db", shlex.quote(database), "-query", shlex.quote(query), "-outfmt", "15")
    logger.info("Executing %s", repr(" ".join(cmd)))
    blastn_json = subprocess.check_output(cmd).decode("utf-8")
    return parse_blastn_json(blastn_json)
