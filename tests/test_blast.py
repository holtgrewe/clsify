import os

from clsify import blast, workflow

import pytest

# Path to directory with files pointed to from Nelson et al., (2012).
PATH_NELSON = os.path.join(os.path.dirname(__file__), "data", "Nelson2012")


def test_revcomp():
    assert blast.revcomp("CGAT") == "ATCG"
    assert blast.revcomp("C-GAT") == "ATC-G"


def test_alignment_revcomp():
    ali = blast.Alignment("CGAT", "| ||", "C-AT")
    expected = blast.Alignment("ATCG", "|| |", "AT-G")
    assert ali.revcomp() == expected


def test_alignment_wrapped():
    ali = blast.Alignment("CGATCGAT", "||||||||", "CGATCGAT")
    expected = [
        "Sbjct    2 CGAT 5",
        "           ||||",
        "Query    1 CGAT 4",
        "",
        "Sbjct    6 CGAT 9",
        "           ||||",
        "Query    5 CGAT 8",
        "",
    ]
    assert ali.wrapped(0, 1, 4) == "\n".join(expected)


def test_alignment_build_empty():
    assert str(blast.Alignment.build_empty()) == "Alignment(hseq='', midline='', qseq='')"


def test_blast_match_build_nomatch():
    assert not blast.BlastMatch.build_nomatch("query", "{}").is_match


def test_match_cigar():
    assert blast.match_cigar("CGATCG-T", "CG-TCGAT", 2, 10, 12) == [
        [2, "H"],
        [2, "M"],
        [1, "I"],
        [3, "M"],
        [1, "D"],
        [1, "M"],
        [2, "H"],
    ]


def test_parse_blastn_json_haplotype_a_16s():
    """Parse ``blastn`` results from Haplotype A specimen for 16S/23S sequence."""
    query = "EU812559.1"
    subject = "EU812559.1"
    path_json = os.path.join(PATH_NELSON, "%s-%s.json" % (query, subject))
    with open(path_json, "rt") as inputf:
        match_json = inputf.read()
    match = blast.parse_blastn_json("%s.fasta" % query, subject, match_json)

    assert match
    assert match.query == query
    assert match.database == subject
    assert match.match_cigar == "2515M"


def test_parse_blastn_json_haplotype_a_16s_nomatch():
    """Parse ``blastn`` results with no match (database was 50S sequence)."""
    query = "EU812559.1"
    subject = "EU834131.1"
    path_json = os.path.join(PATH_NELSON, "%s-%s.json" % (query, subject))
    with open(path_json, "rt") as inputf:
        match_json = inputf.read()
    match = blast.parse_blastn_json("%s.fasta" % query, subject, match_json)

    assert match


def test_parse_blastn_json_haplotype_a_50s():
    """Parse ``blastn`` results from Haplotype A specimen for 50S sequence."""
    query = "EU834131.1"
    subject = "EU834131.1"
    path_json = os.path.join(PATH_NELSON, "%s-%s.json" % (query, subject))
    with open(path_json, "rt") as inputf:
        match_json = inputf.read()
    match = blast.parse_blastn_json("%s.fasta" % query, subject, match_json)

    assert match
    assert match.query == query
    assert match.database == subject
    assert match.match_cigar == "1714M"


def test_run_blast():
    """Parse ``blastn`` results from Haplotype A specimen for 50S sequence."""
    ref_files = list(workflow.REF_FILES.values())
    result = blast.run_blast(ref_files[0], os.path.join(PATH_NELSON, "EU812559.1.fa"))
    assert result
