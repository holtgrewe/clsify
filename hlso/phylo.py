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
from .common import write_fasta


def phylo_analysis(
    df: pd.DataFrame, *, path_out: typing.Optional[str] = None
) -> typing.Dict[str, typing.Dict[str, object]]:
    logger.info("Performing phylogenetics analysis on\n%s", df)
    result = {}

    for region, group in df.groupby("region"):
        seqs = dict(zip(group["query"], group["orig_sequence"]))
        keys = tuple(sorted(group["query"]))
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
            result[region] = {"labels": keys, "dist": dist_cd, "linkage": clustering}

            if path_out:
                region_path_out = path_out % region
                logger.info("Writing Dendrogram to %s", region_path_out)
                plt.figure()
                hierarchy.dendrogram(
                    clustering, labels=keys, orientation="left", link_color_func=lambda x: "black"
                )
                # plt.subplots_adjust(right=0.7)
                plt.suptitle("UPGMA for %s region" % region)
                plt.xlabel("distance [%]")
                ax = plt.gca()
                for k in ("left", "right", "top"):
                    ax.spines[k].set_color("none")
                plt.tight_layout()
                plt.savefig(region_path_out)

    return result
