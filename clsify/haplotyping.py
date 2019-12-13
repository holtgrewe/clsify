"""Code for reading the haplotyping table information.

This contains the informative positions for haplotyping of calls.
"""

import os
import shlex
import subprocess
import tempfile
import typing

import attr
import json
from logzero import logger

from .blast import BlastMatch
from .common import call_variants, normalize_var


@attr.s(auto_attribs=True, frozen=True)
class HaplotypingPos:
    """A haplotyping position."""

    #: database reference sequence name
    reference: str
    #: 0-based position on database reference
    position: int
    #: mapping from haplotype name to value
    haplo_values: typing.Dict[str, str]


def load_haplotyping_table(path: str) -> typing.Dict[typing.Tuple[str, int], HaplotypingPos]:
    """Load haplotyping table from the given ``path``."""
    # logger.debug("Loading haplotyping table from %s", path)
    result = {}

    if not os.path.exists(path):
        logger.warning("Could not load haplotyping table from %s", path)
        return result

    header = None
    with open(path, "rt") as inputf:
        for line in inputf:
            if line.startswith("#"):
                continue
            arr = line.rstrip().split("\t")
            if not header:
                header = arr
            else:
                record = dict(zip(header, arr))
                key = (record["reference"], int(record["pos"]) - 1)
                result[key] = HaplotypingPos(
                    reference=key[0], position=key[1], haplo_values=dict(zip(header[2:], arr[2:]))
                )
    # logger.debug("Done loading %d records", len(result))
    return result


#: The haplotype table.
HAPLOTYPE_TABLE = load_haplotyping_table(
    os.path.join(os.path.dirname(__file__), "data", "haplotype_table.txt")
)
# logger.debug("haplotype table = %s", HAPLOTYPE_TABLE)

#: The haplotype names
HAPLOTYPE_NAMES = "ABCDE"


@attr.s(auto_attribs=True, frozen=True)
class HaplotypingResult:
    """A haplotyping result."""

    #: The file name used for haplotyping
    filename: str
    #: The query name
    query: str
    #: mapping from ``(reference, zero_based_pos)`` to allele value
    informative_values: typing.Dict[typing.Tuple[str, int], str]

    def merge(
        self, other: typing.TypeVar("HaplotypingResult")
    ) -> typing.TypeVar("HaplotypingResult"):
        keys = list(
            sorted(set(self.informative_values.keys()) | set(other.informative_values.keys()))
        )
        merged = {}
        for key in keys:
            if key in self.informative_values and key in other.informative_values:
                here = self.informative_values[key]
                there = other.informative_values[key]
                if here == there:
                    merged[key] = here
            merged[key] = self.informative_values.get(key, other.informative_values.get(key))
        return HaplotypingResult(filename="-", query="-", informative_values=merged)  # post-merging

    def asdict(self, only_summary=False) -> typing.Dict:
        informative = {}
        scores = {}
        for name in HAPLOTYPE_NAMES:
            plus, minus = self.compare(name)
            informative["%s_pos" % name] = plus
            informative["%s_neg" % name] = minus
            scores[name] = plus - minus
        best_score = max(scores.values())
        if best_score > 0:
            best_haplotypes = ",".join(
                [key for key, value in scores.items() if value == best_score]
            )
        else:
            best_haplotypes = "-"
        if only_summary:
            return {"best_haplotypes": best_haplotypes, "best_score": best_score}
        else:
            return {
                "filename": self.filename,
                "best_haplotypes": best_haplotypes,
                "best_score": best_score,
                **informative,
                **{
                    "%s:%d" % (key[0], key[1] + 1): self.informative_values.get(key)
                    for key in HAPLOTYPE_TABLE
                },
            }

    def compare(self, haplotype: str) -> typing.Tuple[int, int]:
        """Return ``(match_count, mismatch_count)`` for the given ``haplotype``."""
        positive = 0
        negative = 0
        for key, value in self.informative_values.items():
            if HAPLOTYPE_TABLE[key].haplo_values[haplotype] == value:
                positive += 1
            else:
                negative += 1
        return (positive, negative)

    @classmethod
    def fromdict(self, dict_: typing.Dict) -> typing.TypeVar("HaplotypingResult"):
        informative_values = {}
        for key, value in dict_.items():
            if ":" in key and value is not None:
                arr = key.split(":", 1)
                informative_values[(arr[0], int(arr[1]) - 1)] = value
        return HaplotypingResult(
            filename=dict_["filename"], query=dict_["query"], informative_values=informative_values
        )


@attr.s(auto_attribs=True, frozen=True)
class HaplotypingResultWithMatches:
    """Result from haplotyping."""

    #: The actual ``HaplotypingResult``.
    result: typing.Optional[HaplotypingResult]
    #: The ``BlastMatch``es that the haplotyping is based on.
    matches: typing.Optional[typing.Tuple[BlastMatch]]

    @staticmethod
    def build_empty() -> typing.TypeVar("HaplotypingResultWithMatches"):
        return HaplotypingResultWithMatches(None, None)


def run_haplotyping(
    matches: typing.Iterable[BlastMatch]
) -> typing.Dict[str, HaplotypingResultWithMatches]:
    """Perform the haplotyping based on the match."""
    results_matches = {}
    results_haplo = {}

    # TODO: properly handle overlapping changes; directly look into alignment
    for match in matches:
        ref = match.database
        if "_" in ref:
            ref = ref.split("_")[0]

        calls = call_variants(match.alignment.hseq, match.alignment.qseq, match.database_start)

        informative_values = {}
        for (h_ref, h_pos), variant in HAPLOTYPE_TABLE.items():
            if ref == h_ref and h_pos >= match.database_start and h_pos < match.database_end:
                if h_pos + 1 in calls:
                    informative_values[(ref, h_pos)] = variant.haplo_values["alt"]
                else:
                    informative_values[(ref, h_pos)] = variant.haplo_values["ref"]

        result = HaplotypingResult(
            filename=match.path, query=match.query, informative_values=informative_values
        )
        if result.filename in results_haplo:
            results_matches[result.filename].append(match)
            results_haplo[result.filename] = results_haplo[result.filename].merge(result)
        else:
            results_matches[result.filename] = [match]
            results_haplo[result.filename] = result

    return {
        filename: HaplotypingResultWithMatches(
            result=results_haplo[filename], matches=results_matches[filename]
        )
        for filename in results_matches
    }
