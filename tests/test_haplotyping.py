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
    for record in SAMPLES:
        path_json = os.path.join(PATH_DATA, record["folder"], "%s.json" % record["accession"])
        with open(path_json, "rt") as inputf:
            str_json = inputf.read()
        matches = blast.parse_blastn_json(str_json)
        haplo_results = haplotyping.run_haplotyping(matches)
        assert len(haplo_results) == 1
        haplo_result = list(haplo_results.values())[0]
        res = haplo_result.asdict()
        samples_seen.add(record["accession"])
        if res["best_haplotypes"] != "-":
            samples_haplotyped.add(record["accession"])
            assert record["haplotype"] in res["best_haplotypes"], record["accession"]
    assert samples_seen == samples_haplotyped
