"""Code for reading the haplotyping table information.

This contains the informative positions for haplotyping of calls.
"""

import os
import shlex
import subprocess
import tempfile
import typing

import attr
from logzero import logger
import vcfpy


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
    header = None
    with open(path, "rt") as inputf:
        for line in inputf:
            arr = line.rstrip().split("\t")
            if not header:
                header = arr
            else:
                record = dict(zip(header, arr))
                key = (record["reference"], int(record["position"]) - 1)
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
HAPLOTYPE_NAMES = ("A", "B", "C", "D", "E")


@attr.s(auto_attribs=True, frozen=True)
class HaplotypingResult:
    """A haplotyping result."""

    #: The file name used for haplotyping
    filename: str
    #: mapping from ``(reference, zero_based_pos)`` to allele value
    informative_values: typing.Dict[typing.Tuple[str, int], str]

    def merge(self, other):
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
        return HaplotypingResult(filename="-", informative_values=merged)  # post-merging

    def asdict(self, only_summary=False):
        informative = {}
        scores = {}
        for name in HAPLOTYPE_NAMES:
            plus, minus = self.compare(name)
            informative["%s_pos" % name] = plus
            informative["%s_neg" % name] = minus
            scores[name] = plus - minus
        best_score = max(scores.values())
        best_haplotypes = ",".join([key for key, value in scores.items() if value == best_score])
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

    def compare(self, haplotype):
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
    def fromdict(self, dict_):
        informative_values = {}
        for key, value in dict_.items():
            if ":" in key and value is not None:
                arr = key.split(":", 1)
                informative_values[(arr[0], int(arr[1]) - 1)] = value
        return HaplotypingResult(filename=dict_["filename"], informative_values=informative_values)


def run_bcftools(ref_path, path_sam, tmpdir):
    """Perform variant calling."""
    logger.info("running variant calling on SAM file %s", path_sam)
    cmd_mpileup = (
        "bcftools",
        "mpileup",
        "--fasta-ref",
        shlex.quote(ref_path),
        shlex.quote(path_sam),
    )
    cmd_call = ("bcftools", "call", "-c", "--ploidy", "1")
    with subprocess.Popen(cmd_mpileup, stdout=subprocess.PIPE) as mpileup:
        path_vcf = os.path.join(tmpdir, path_sam + ".vcf")
        with open(path_vcf, "wb") as vcf_file:
            with subprocess.Popen(cmd_call, stdin=mpileup.stdout, stdout=vcf_file):
                pass
        logger.info("=> VCF file %s", path_vcf)
    return path_vcf


def run_haplotyping(path_ref, blast_matches):
    """Perform haplotying based on the ``BlastMatch`` objects in ``blast_matches``."""
    logger.info("Reading reference FAI...")
    refs = {}
    with open("%s.fai" % path_ref, "rt") as fai:
        for line in fai:
            arr = line.rstrip().split("\t")
            refs[arr[0]] = {"name": arr[0], "length": int(arr[1])}
    logger.info("=> %d reference sequences", len(refs))

    logger.info("Performing haplotyping for %d BLAST matches...", len(blast_matches))
    result = {}
    with tempfile.TemporaryDirectory() as tmpdir:
        for match in blast_matches:
            logger.info("Haplotyping for %s", match.query)
            # Write out BLAST match as SAM file.
            path_sam = os.path.join(tmpdir, "tmp.sam")
            with open(path_sam, "wt") as sam:
                print("@HD\tVN:1.6\tSO:coordinate", file=sam)
                for ref in refs.values():
                    print("@SQ\tSN:%s\tLN:%d" % (ref["name"], ref["length"]), file=sam)
                if match.is_match:
                    print("\t".join(match.to_sam_tsv()), file=sam)

            # Perform variant calling.
            path_vcf = run_bcftools(path_ref, path_sam, tmpdir)

            # Extract sample bases.
            sample_bases = {}
            with vcfpy.Reader.from_path(path_vcf) as vcf_file:
                for record in vcf_file:
                    key = (record.CHROM, record.POS - 1)
                    if key in HAPLOTYPE_TABLE:
                        sample_bases[key] = record.calls[0].gt_bases[0]
            # Store results.
            result[match.query] = HaplotypingResult(match.query, sample_bases)

    logger.info("=> done haplotyping")
    return result


if __name__ == "__main__":
    from .blast import run_blast

    ref = "clsify/data/EU834131.1.fasta"
    match = run_blast(ref, "examples/CA_SanJoaquinValley_Tm4_16S_LsoF.fasta")
    # print(match)
    res = run_haplotyping(ref, [match])
    # print([r.asdict() for r in res.values()])
    # print([HaplotypingResult.fromdict(r.asdict()) for r in res.values()])
