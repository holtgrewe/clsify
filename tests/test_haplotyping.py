import os

from clsify import blast, haplotyping, workflow

import pytest

# Path to the data directory.
PATH_DATA = os.path.join(os.path.dirname(__file__), "data")
# Path to the samples.tsv file.
PATH_SAMPLES = os.path.join(PATH_DATA, "samples.tsv")
# Path to directory with files pointed to from Nelson et al., (2012).
PATH_NELSON = os.path.join(PATH_DATA, "Nelson2012")


def load_samples(path):
    """Load ```samples.tsv``` file."""
    header = None
    result = []
    with open(path, "rt") as inputf:
        for line in inputf:
            if line.startswith("#"):
                continue
            arr = line.strip().split("\t")
            if header is None:
                header = arr
            else:
                result.append(dict(zip(header, arr)))
    return result


# A list of samples, loaded from ``samples.tsv``.***
SAMPLES = load_samples(PATH_SAMPLES)


def test_haplotyping():
    """Test haplotyping for all samples."""
    samples_seen = set()
    samples_haplotyped = set()
    bad = []
    for ref in workflow.REF_FILES:
        for record in SAMPLES:
            path_json = os.path.join(
                PATH_DATA, record["folder"], "%s-%s.json" % (record["accession"], ref)
            )
            with open(path_json, "rt") as inputf:
                str_json = inputf.read()
            match = blast.parse_blastn_json(ref, "%s.fa" % record["accession"], str_json)
            haplo_result = haplotyping.run_haplotyping(match)
            res = haplo_result.asdict()
            samples_seen.add(record["accession"])
            if res["best_haplotypes"] != "-":
                samples_haplotyped.add(record["accession"])
                assert record["haplotype"] in res["best_haplotypes"], record["accession"]
    assert samples_seen == samples_haplotyped


def test_haplotyping_result():
    results = []
    for ref in workflow.REF_FILES:
        for accession in ["EU812559.1", "EU834131.1"]:
            path_json = os.path.join(
                PATH_DATA, "Nelson2012", "%s-%s.json" % (accession, ref)
            )
            with open(path_json, "rt") as inputf:
                str_json = inputf.read()
            match = blast.parse_blastn_json(ref, "%s.fa" % accession, str_json)
            haplo_result = haplotyping.run_haplotyping(match)
            assert len(haplo_result.asdict()) > 2
            assert len(haplo_result.asdict(only_summary=True)) == 2
            results.append(haplo_result)
    merged = results[0].merge(results[3])
    assert results[0].asdict(True) == {'best_haplotypes': 'A,C', 'best_score': 14}
    assert results[1].asdict(True) == {'best_haplotypes': '-', 'best_score': 0}
    assert results[2].asdict(True) == {'best_haplotypes': '-', 'best_score': 0}
    assert results[3].asdict(True) == {'best_haplotypes': 'A', 'best_score': 14}
    assert merged.asdict(True) == {'best_haplotypes': 'A', 'best_score': 28}
