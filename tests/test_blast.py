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
    path_json = os.path.join(PATH_NELSON, "EU812559.1.json")
    with open(path_json, "rt") as inputf:
        match_json = inputf.read()

    matches = blast.parse_blastn_json(match_json)

    assert matches[0].query == "EU812559.1"
    assert matches[0].database.endswith("/16S-23S")
    assert matches[0].match_cigar == "1225H1290M"

    assert matches[1].query == "EU812559.1"
    assert matches[1].database.endswith("/16S")
    assert matches[1].match_cigar == "1225M1290H"


def test_run_blast():
    """Parse ``blastn`` results from Haplotype A specimen for 50S sequence."""
    result = blast.run_blast(workflow.REF_FILE, os.path.join(PATH_NELSON, "EU812559.1.fa"))
    assert len(result) == 2
    assert result[0].database.endswith("/16S-23S")
    assert result[1].database.endswith("/16S")
