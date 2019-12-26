"""Implementation of the workflow for (single and batches of) reads.

Each read is aligned to both reference sequences.  The best match is taken and haplotyping is
performed by considering variants at the informative positions if any match could be found.
The result is a pair of ``blast.BlastMatch`` and ``haplotyping.HaplotypingResult`` (or ``None``)
that can be joined by the read file name.  Further, the sample information can be derived
from this with a regexp.
"""

import os
import re
import typing

import attr
from logzero import logger
import pandas as pd
import tempfile

from .blast import run_blast, BlastMatch
from .common import load_fasta
from .haplotyping import run_haplotyping, HaplotypingResultWithMatches

#: Default minimal quality to consider a match as true.
DEFAULT_MIN_IDENTITY = 0.5

#: Default regular expression to use for inferring the sample information.
# TODO: change a bit...?
DEFAULT_PARSE_RE = r"^(?P<sample>[^_]+_[^_]+_[^_]+)_(?P<primer>.*?)\.fasta"

#: The reference files to use.
REF_FILE = os.path.join(os.path.dirname(__file__), "data", "ref_seqs.fasta")


@attr.s(auto_attribs=True, frozen=True)
class NamedSequence:
    """A named sequence."""

    #: the sequence name
    name: str
    #: the sequence
    sequence: str


def only_blast(path_query: str) -> typing.Tuple[BlastMatch]:
    """Run BLAST and haplotyping for the one file at ``path_query``."""
    logger.info("Running BLAST on all references for %s...", path_query)
    return run_blast(REF_FILE, path_query)


def blast_and_haplotype(path_query: str) -> typing.Dict[str, HaplotypingResultWithMatches]:
    return run_haplotyping(only_blast(path_query))


def blast_and_haplotype_many(
    paths_query: typing.Iterable[str]
) -> typing.Dict[str, HaplotypingResultWithMatches]:
    """Run BLAST and haplotyping for all files at ``paths_query``.

    Return list of dicts with keys "best_match" and "haplo_result".
    """
    logger.info("Running BLAST and haplotyping for all queries...")
    result = {}
    for path_query in paths_query:
        path_result = blast_and_haplotype(path_query)
        if path_result:
            result.update(path_result)
        else:
            result[path_query] = HaplotypingResultWithMatches.build_empty()
    return result


def strip_ext(s: str) -> str:
    return s.rsplit(".", 1)[0]


def results_to_data_frames(
    results: typing.Dict[str, HaplotypingResultWithMatches], regex: str, column: str = "query"
) -> typing.Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Convert list of dicts with best_match/haplo_result to triple of Pandas DataFrame.

    The three DataFrame will contain the following information:

    1. A summary data frame showing best BLAST match target and identity plus haplotype.
    2. A data frame showing BLAST result details.
    3. A data frame showing haplotyping result details.
    """
    r_summary = []
    r_blast = []
    r_haplo = []
    for path, result in results.items():
        haplo_result = result.result
        if not haplo_result:
            for query, query_seq in load_fasta(path).items():
                r_summary.append(
                    {"query": query, "database": ".", "identity": 0, "orig_sequence": query_seq}
                )
                r_blast.append({"query": query})
                r_haplo.append({"query": query})
        else:
            query_seq = load_fasta(path)[haplo_result.query]
            haplo_matches = result.matches
            best_match = list(sorted(haplo_matches, key=lambda m: m.identity, reverse=True))[0]

            r_summary.append(
                {
                    "query": best_match.query,
                    "database": best_match.database,
                    "identity": 100.0 * best_match.identity,
                    **{
                        key: value
                        for key, value in haplo_result.asdict().items()
                        if key in ("best_haplotypes", "best_score")
                    },
                    "orig_sequence": query_seq,
                }
            )
            r_blast.append(
                {
                    "query": best_match.query,
                    "database": best_match.database,
                    "identity": 100.0 * best_match.identity,
                    "q_start": best_match.query_start,
                    "q_end": best_match.query_end,
                    "q_str": best_match.query_strand,
                    "db_start": best_match.database_start,
                    "db_end": best_match.database_end,
                    "db_str": best_match.database_strand,
                    "alignment": best_match.alignment.wrapped(
                        best_match.query_start, best_match.database_start
                    ),
                    "orig_sequence": query_seq,
                }
            )
            r_haplo.append(
                {
                    "query": best_match.query,
                    **{
                        key: value
                        for key, value in haplo_result.asdict().items()
                        if "_pos" in key
                        or "_neg" in key
                        or key in ("best_haplotypes", "best_score")
                    },
                }
            )

    dfs = pd.DataFrame(r_summary), pd.DataFrame(r_blast), pd.DataFrame(r_haplo)
    dfs = list(map(lambda df: match_sample_in_data_frame(df, regex, column), dfs))
    dfs[0] = augment_summary(
        dfs[0], results, regex, column, "sample" if "sample" in dfs[0].columns else "query"
    )
    for df in dfs:
        df.index = range(df.shape[0])
        df.insert(0, "id", df.index)
    return tuple(dfs)


SUMMARY_SUFFIX = "_ZZZ *** SUMMARY ***"


def augment_summary(
    df: pd.DataFrame,
    results: typing.Dict[str, HaplotypingResultWithMatches],
    regex: str,
    column: str,
    group_by: str,
):
    grouped = {}
    for record in results.values():
        if not record.result:
            continue
        m = re.match(regex, record.result.query)
        if m and m.group(group_by):
            if m.group(group_by) not in grouped:
                grouped[m.group(group_by)] = record.result
            else:
                grouped[m.group(group_by)] = grouped[m.group(group_by)].merge(record.result)
    rows = []
    for key, value in grouped.items():
        rows.append(
            {
                "query": "%s%s" % (key, SUMMARY_SUFFIX),
                group_by: key,
                **value.asdict(only_summary=True),
            }
        )
    for key in set(df[group_by].values) - grouped.keys():  # fill for those without matches
        rows.append({"query": "%s%s" % (key, SUMMARY_SUFFIX), group_by: key})
    orig_columns = list(df.columns.values)
    if "samples" in df.columns:
        keys = ["sample", "query"]
    else:
        keys = ["query"]
    df = df.append(pd.DataFrame(rows))[orig_columns].sort_values(keys).fillna("-")
    df["query"] = df["query"].str.replace(re.escape(SUMMARY_SUFFIX), "")
    return df


def match_sample_in_data_frame(df, regex, column):
    """Use ``regex`` to parse out information from the ``"query"`` column of the ``df`` DataFrame.

    Return an augmented DataFrame.
    """
    if not df.shape[0]:
        return df  # short-circuit empty
    # Get shortcut to column with query/filenames.
    col_query = df.loc[:, column]

    # Obtain list of new column names, sorted by occurence, must use dict for that.
    names = {}
    for query in col_query:
        m = re.match(regex, query)
        if m:
            for key in m.groupdict():
                names[key] = True
    names = list(names.keys())

    # Build new column values.
    columns = {n: [] for n in names}
    for query in col_query:
        m = re.match(regex, query)
        if m:
            for name in names:
                columns[name].append(m.groupdict().get(name))
        else:
            for name in names:
                columns[name].append(None)

    # Insert new columns into df.
    idx = df.columns.get_loc(column)
    for i, (key, column) in enumerate(columns.items()):
        df.insert(idx + i + 1, key, column)

    return df
