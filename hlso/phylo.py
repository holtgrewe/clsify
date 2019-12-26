"""Phylogenetics analysis"""

import os
import tempfile
import typing

from logzero import logger
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.spatial import distance
from scipy.cluster import hierarchy

from .blast import run_blast, run_makeblastdb
from .common import write_fasta, load_fasta


def phylo_analysis(
    df: pd.DataFrame, *, path_out: typing.Optional[str] = None
) -> typing.Dict[str, typing.Dict[str, object]]:
    logger.info("Performing phylogenetics analysis on\n%s", df)
    result = {}

    for region, group in df.groupby("region"):
        seqs = dict(zip(group["query"], group["orig_sequence"]))
        path_ref = os.path.join(os.path.dirname(__file__), "reference/%s.fasta" % region)
        if os.path.exists(path_ref):
            logger.info("Loading reference %s", path_ref)
            seqs.update(load_fasta(path_ref))
        keys = tuple(sorted(seqs.keys()))

        if len(seqs) == 1:  # skip if only one sequence given
            result[region] = {"message": "Cannot compute dendrogram for one sequence."}
            continue

        with tempfile.TemporaryDirectory() as tmp_dir:
            logger.info("Performing phylogenetics analysis for %s region", region)
            path_seqs = os.path.join(tmp_dir, "seqs.fasta")
            with open(path_seqs, "wt") as tmpf:
                write_fasta(seqs, file=tmpf)
                tmpf.flush()

            logger.info("Running all-to-all BLAST")
            run_makeblastdb(path_seqs)
            results = run_blast(path_seqs, path_seqs)

            logger.info("Performing clustering and dendrogram...")
            identities = sorted((m.database, m.query, m.identity) for m in results)
            difference = list(map(lambda x: 100.0 * (1.0 - x[2]), identities))
            dist_sq = np.asarray(difference, dtype=np.float64).reshape(len(keys), len(keys))
            dist_cd = distance.squareform(dist_sq)
            clustering = hierarchy.average(dist_cd)
            result[region] = {
                "labels": keys,
                "dist": dist_cd.tolist(),
                "linkage": clustering.tolist(),
            }

            if path_out:
                plot_phylo(clustering, keys, region, path_out % region)

    return result


def plot_phylo(linkage, labels, region, fname, *, format=None, dpi=None):
    plt.figure()
    hierarchy.dendrogram(
        linkage, labels=labels, orientation="left", link_color_func=lambda x: "black"
    )
    plt.suptitle("UPGMA for %s region" % region)
    plt.xlabel("difference [%]")
    ax = plt.gca()
    for k in ("left", "right", "top"):
        ax.spines[k].set_color("none")
    plt.tight_layout()
    plt.savefig(fname, format=format, dpi=dpi)
