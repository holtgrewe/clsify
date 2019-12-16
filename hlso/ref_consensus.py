"""Perform consensus building for creating the reference step.

The following actions are performed for this step:

1. Build consensus sequence for each seed:  for each seed, consider the NCBI WWWBLAST result
   and build replace the seed segment by the database hit and perform an end-to-end MSA for
   this for each seed sequence.Generate consensus seed sequence for this.

2. Perform a BLAST search of the consensus seed sequences and generate type reference sequences
   for each haplotype.

3. Perform MSA of the type reference sequences and generate a haplotyping table from this.
"""

# TODO: replace dict representing variant by attrs-based class.

import itertools
import os
import shlex
import shutil
import subprocess
import tempfile
import textwrap

from Bio.Blast import NCBIXML
from Bio.Align import Applications
from logzero import logger

from .ref_download import REF_SEQS
from .ref_blast import load_tsv
from .cli import _proc_args
from .common import call_variants, describe, load_fasta, normalize_var, only_bases
from .paste import REF_FILE, do_paste
from .workflow import only_blast


def paste_query_seq(query_seq, blast_alignment):
    """Paste database sequence from ``blast_alignment`` into ``query_seq``."""
    hsp = blast_alignment.hsps[0]

    if hsp.strand[0] == "Plus":  # query strand is forward
        prefix = query_seq[: hsp.query_start]
        sbjct_seq = "".join([c for c in hsp.sbjct if c in "cgatnCGATN"])
        suffix = query_seq[hsp.query_end :]
        seq = prefix + sbjct_seq + suffix
    else:
        raise Exception("Cannot handle matches on query reverse strand")

    return "".join([c for c in seq if c != "N"]).upper()


def consensus(cs):
    counter = {}
    for c in cs:
        counter.setdefault(c, 0)
        counter[c] += 1
    max_val = max(counter.values())
    for key, value in counter.items():
        if value == max_val:
            return key
    return "N"


def write_consensus(path, query_name, query_seq, alignments, args):
    if os.path.exists(path):
        logger.warn("Consensus sequence exists, skipping: %s", path)
        return

    logger.info("Creating consensus for %s", query_name)
    single_hsps = [a for a in alignments if len(a.hsps) == 1]
    logger.info("- removed %d alignments with HSP count != 1", len(alignments) - len(single_hsps))
    sorted_alis = sorted(
        single_hsps, key=lambda a: a.hsps[0].identities / a.hsps[0].align_length, reverse=True
    )

    good_alis = [
        a for a in sorted_alis if a.hsps[0].identities >= a.hsps[0].align_length - args.max_errors
    ]
    logger.info("- picked %d alignments with up to %d errors", len(good_alis), args.max_errors)
    logger.info("  %s", [a.title.split("|")[3] for a in good_alis])
    if not good_alis:
        logger.error("Could not pick any alignments for %s. Skipping", query_name)
        return

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_fasta = os.path.join(tmpdir, "seqs.fasta")
        logger.info("Writing to %s...", tmp_fasta)
        with open(tmp_fasta, "wt") as tmpf:
            for ali_count, ali in enumerate(good_alis):
                pasted_seq = "\n".join(textwrap.wrap(paste_query_seq(query_seq, ali)))
                print(">%s\n%s" % (ali.title.split("|")[3], pasted_seq), file=tmpf)

        path_alignment = path[: -len(".consensus.fasta")] + ".clustalw.fasta"
        if os.path.exists(path_alignment):
            logger.warn("Already exists: %s", path_alignment)
            logger.warn(" => not running ClustalW")
        elif ali_count == 0:
            logger.warn("Only one sequence, no need to align, just copy")
            shutil.copy(tmp_fasta, path_alignment)
        else:
            logger.info("Running ClustalW")
            cmd = [
                "clustalw",
                "-INFILE=%s" % tmp_fasta,
                "-ALIGN",
                "-OUTFILE=%s" % path_alignment,
                "-OUTPUT=FASTA",
                "-OUTORDER=INPUT",
                "-SEED=42",
            ]
            logger.info("- %s", " ".join(cmd))
            res = subprocess.check_call(cmd)
            if res != 0:
                logger.error("Runing ClustalW failed with retval %d", res)
            else:
                logger.info("Done running ClustalW!")

    logger.info("Computing consensus sequence to %s", path)
    fasta = load_fasta(path_alignment)
    consensus_seq = "".join(filter(lambda x: x in "ACGTN", map(consensus, zip(*fasta.values()))))
    with open(path, "wt") as outputf:
        print(
            ">consensus-%s\n%s" % (query_name, "\n".join(textwrap.wrap(consensus_seq))),
            file=outputf,
        )


def build_seed_consensus(records, args):
    """Build consensus sequence for each seed."""
    logger.info("Building seed consensus sequences...")

    base_dir = os.path.dirname(args.in_tsv)

    for record in records:
        logger.info("Handling %s", record["path"])
        basename = os.path.basename(record["path"])[: -len(".fasta")]
        fasta = load_fasta(os.path.join(base_dir, basename + ".fasta"))

        path_blast = os.path.join(base_dir, basename + ".blast.xml")
        with open(path_blast, "rt") as inputf:
            for record in NCBIXML.parse(inputf, debug=1 if args.verbose else 0):
                # Compute consensus for each query sequence if it does not exist already.
                path_consensus = os.path.join(base_dir, basename + ".consensus.fasta")
                if os.path.exists(path_consensus):
                    logger.warn("Skipping, already exists: %s", path_consensus)
                    continue
                query_name = record.query.split()[0]
                write_consensus(
                    path=path_consensus,
                    query_name=query_name,
                    query_seq=fasta[query_name],
                    alignments=record.alignments,
                    args=args,
                )


def build_haplotype_sequences(records, args):
    """Build consensus sequence for each seed."""
    logger.info("Building building haplotype sequences...")

    ref_seqs = load_fasta(REF_FILE)
    base_dir = os.path.dirname(args.in_tsv)

    results = {}  # actually consensus seeds
    for record in records:
        out_path = os.path.join(base_dir, record["region"] + ".fasta")
        if os.path.exists(out_path):
            logger.warn("Output path exists, skipping: %s", out_path)
            continue

        key = (record["haplotype"], record["region"])
        logger.info("Handling %s", record["path"])
        basename = os.path.basename(record["path"])[: -len(".fasta")]
        matches = only_blast(os.path.join(base_dir, basename + ".consensus.fasta"))
        for match in matches:
            if match.database is None:
                logger.info("  => no match for %s", match.query)
            else:
                logger.info("Pasting region=%s haplotype=%s", record["region"], record["haplotype"])
                seq = do_paste(match, ref_seqs)
                results.setdefault(record["region"], {})[record["haplotype"]] = seq

    for region, region_results in results.items():
        fasta = os.path.join(base_dir, region + ".fasta")
        logger.info("Writing padded consensus seed to to %s", fasta)
        with open(fasta, "wt") as outputf:
            for haplotype, seq in region_results.items():
                wrapped = textwrap.wrap(seq)
                print(">%s_%s\n%s" % (region, haplotype, "\n".join(wrapped)), file=outputf)

    return results


def strip_n(s):
    return "".join([c for c in s if c not in "nN"])


def build_haplotyping_table(records, args):
    """Build consensus sequence for each seed."""
    logger.info("Building haplotyping table...")

    ref_seqs = load_fasta(REF_FILE)
    base_dir = os.path.dirname(args.in_tsv)

    ali_input = {
        key.split("_")[1]: {key.split("_")[1]: strip_n(value)} for key, value in ref_seqs.items()
    }

    for region, seqs in ali_input.items():
        fasta_path = os.path.join(base_dir, "%s.fasta" % region)
        if not os.path.exists(fasta_path):
            logger.warn("Consensus sequences FASTA not found, skipping: %s", fasta_path)
        else:
            logger.info("Loading consensus sequences FASTA: %s", fasta_path)
            seqs.update(load_fasta(fasta_path))

    result = {}

    for region, seqs in ali_input.items():
        for seq_name, ref_seq in ref_seqs.items():
            ref_name, ref_region = seq_name.split("_", 1)
            if ref_region == region:
                ref_offset = 0
                for ref_offset, c in enumerate(ref_seq):
                    if c == "N":
                        break
                break  # => ref_name, ref_seq, ref_offset are set, else Exception is raised
        else:
            raise Exception("Could not find reference for %s" % ref_name)

        path_alignment = os.path.join(base_dir, "%s.clustalw.fasta" % region)
        if len(seqs) == 1:
            logger.warn("Only reference and no consensus, skipping: %s", path_alignment)
            continue

        logger.info("Running ClustalW on consensus alignments for %s", region)
        if os.path.exists(path_alignment):
            logger.warn("Ref+consensus alignment already exists, skipping: %s", path_alignment)
        else:
            with tempfile.TemporaryDirectory() as tmpdir:
                tmp_fasta = os.path.join(tmpdir, "seqs.fasta")
                logger.info("Writing to %s...", tmp_fasta)
                with open(tmp_fasta, "wt") as tmpf:
                    for name, seq in seqs.items():
                        print(">%s\n%s" % (name, seq), file=tmpf)

                logger.info("Running ClustalW")
                cmd = [
                    "clustalw",
                    "-INFILE=%s" % tmp_fasta,
                    "-ALIGN",
                    "-OUTFILE=%s" % path_alignment,
                    "-OUTPUT=FASTA",
                    "-OUTORDER=INPUT",
                    "-SEED=42",
                ]
                logger.info("- %s", " ".join(cmd))
                res = subprocess.check_call(cmd)
                if res != 0:
                    logger.error("Runing ClustalW failed with retval %d", res)
                else:
                    logger.info("Done running ClustalW!")

        logger.info("Loading MSA from %s", path_alignment)
        alignment = load_fasta(path_alignment)

        logger.info("Generating variant table for each haplotype")
        variants = {
            key: call_variants(alignment[region], alignment[key], REF_SEQS[region]["start"])
            for key in sorted(alignment.keys())
            if key != region
        }
        xs = set()
        for calls in variants.values():
            for vals in calls.values():
                pos = vals["pos"]
                ref = only_bases(vals["ref_bases"])
                alt = only_bases(vals["alt_bases"])
                desc = vals["description"]
                xs.add((pos, ref, alt, desc))

        table = {}
        haplotypes = [h.split("_")[1] for h in variants.keys()]

        header = ["row", "pos", "ref_bases", "alt_bases", "description"] + haplotypes
        lines = []
        for row, (pos, ref, alt, desc) in enumerate(sorted(xs), 1):
            ref_seq, region = seq_name.split("_", 1)
            record = {"reference": ref_seq, "region": region, "position": pos}
            line = [row, pos, ref, alt, desc]

            # Find reference allele bases.
            ref = "N"
            for name in variants.keys():
                if pos in variants[name] and variants[name][pos]["description"] == desc:
                    ref = only_bases(variants[name][pos]["ref"])
                    break

            for name in variants.keys():
                if pos in variants[name] and variants[name][pos]["description"] == desc:
                    record[name.split("_", 1)[1]] = only_bases(variants[name][pos]["alt"])
                    line.append("+")
                else:
                    record[name.split("_", 1)[1]] = only_bases(ref)
                    line.append("-")
            table[(pos, ref, alt, desc)] = record
            lines.append(line)

        result[(ref_name, ref_region)] = dict(sorted(table.items()))

        print("Haplotyping Table")
        print("-----------------")
        print()
        print("\t".join(map(str, header)).replace("_bases", "").upper() + "\n")
        for line in lines:
            print("\t".join(map(str, line)))

    if args.output_table:
        logger.info("Writing haplotyping table to %s", args.output_table)
        separator = (
            "# --------\t-------\t------\t------\t------\t---------------"
            + "\t------" * len(haplotypes)
        )
        with open(args.output_table, "wt") as haplof:
            print(separator, file=haplof)
            print(
                "\t".join(["reference", "region", "pos", "ref", "alt", "description"] + haplotypes),
                file=haplof,
            )
            print(separator, file=haplof)
            for (ref_name, ref_region), records in result.items():
                for (pos, ref, alt, desc), record in records.items():
                    line = [ref_name, ref_region, pos, ref, alt, desc] + [
                        record[h] for h in haplotypes
                    ]
                    print("\t".join(map(str, line)), file=haplof)
                print(separator, file=haplof)
    else:
        logger.info("Writing no haplotype table (use --output-table to specify the path")


def run(parser, args):
    """Perform the consensus computation."""
    logger.info("Finishing reference computation.")
    logger.info("Arguments are %s", vars(args))

    _header, records = load_tsv(args.in_tsv)

    build_seed_consensus(records, args)
    build_haplotype_sequences(records, args)
    build_haplotyping_table(records, args)

    logger.info("All done. Have a nice day!")


def add_parser(subparser):
    """Configure the ``argparse`` sub parser."""
    parser = subparser.add_parser("ref_consensus")
    parser.set_defaults(func=run)

    parser.add_argument("--verbose", action="store_true", default=False, help="Enable verbose mode")
    parser.add_argument(
        "--max-errors",
        type=int,
        default=0,  # TODO: make sequence-specific!
        help="Maximal number of mismatches to accept in seed consensus computation",
    )
    parser.add_argument("--output-table", default=None, help="Path to output haplotype table.")
    parser.add_argument("in_tsv", help="Path to output TSV file.")
