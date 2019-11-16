"""Implementation of the workflow for (single and batches of) reads.

Each read is aligned to both reference sequences.  The best match is taken and haplotyping is
performed by considering variants at the informative positions if any match could be found.
The result is a pair of ``blast.BlastMatch`` and ``haplotyping.HaplotypingResult`` (or ``None``)
that can be joined by the read file name.  Further, the sample information can be derived
from this with a regexp.
"""

import os
import re

import attr
from logzero import logger
import pandas as pd

from .blast import run_blast
from .haplotyping import run_haplotyping

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


def blast_and_haplotype(path_query):
    """Run BLAST and haplotyping for the one file at ``path_query``."""
    logger.info("Running BLAST on all references for %s...", path_query)
    with open(path_query, "rt") as inputf:
        line_name = None
        lines_seq = []
        for line in inputf:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if line_name is not None:
                    break
                else:
                    line_name = line
            elif line_name is not None:
                lines_seq.append(line)
    matches = run_blast(REF_FILE, path_query)
    haplo_result = run_haplotyping(matches)
    return matches, haplo_result


def blast_and_haplotype_many(paths_query):
    """Run BLAST and haplotyping for all files at ``paths_query``.

    Return list of dicts with keys "best_match" and "haplo_result".
    """
    logger.info("Running BLAST and haplotyping for all queries...")
    result = []
    for path_query in sorted(paths_query):
        matches, haplo_result = blast_and_haplotype(path_query)
        result.append({"matches": matches, "haplo_result": haplo_result})
    return result


def results_to_data_frames(results, regex, column="query"):
    """Convert list of dicts with best_match/haplo_result to triple of Pandas DataFrame.

    The three DataFrame will contain the following information:

    1. A summary data frame showing best BLAST match target and identity plus haplotype.
    2. A data frame showing BLAST result details.
    3. A data frame showing haplotyping result details.
    """
    r_summary = []
    r_blast = []
    r_haplo = []
    for result in results:
        query = result["query"]
        best_match = result["best_match"]
        haplo_result = result["haplo_result"]

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
                "orig_sequence": query.sequence,
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
                "orig_sequence": query.sequence,
            }
        )
        r_haplo.append(
            {
                "query": haplo_result.filename,
                **{
                    key: value
                    for key, value in haplo_result.asdict().items()
                    if "_pos" in key or "_neg" in key or key in ("best_haplotypes", "best_score")
                },
            }
        )

    dfs = pd.DataFrame(r_summary), pd.DataFrame(r_blast), pd.DataFrame(r_haplo)
    dfs = list(map(lambda df: match_sample_in_data_frame(df, regex, column), dfs))
    dfs[0] = augment_summary(dfs[0], results, regex, column, "sample")
    for df in dfs:
        df.index = range(df.shape[0])
        df.insert(0, "id", df.index)
    return dfs


def augment_summary(df, results, regex, column, group_by):
    grouped = {}
    for pair in results:
        haplo_result = pair["haplo_result"]
        m = re.match(regex, haplo_result.filename)
        if m and m.group(group_by):
            if m.group(group_by) not in grouped:
                grouped[m.group(group_by)] = haplo_result
            else:
                grouped[m.group(group_by)] = grouped[m.group(group_by)].merge(haplo_result)
    rows = []
    for key, value in grouped.items():
        rows.append({group_by: key, **value.asdict(only_summary=True)})
    orig_columns = list(df.columns.values)
    return df.append(pd.DataFrame(rows))[orig_columns].sort_values(["sample", "query"]).fillna("-")


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


if __name__ == "__main__":
    results = blast_and_haplotype_many(
        [
            "examples/CA_SanJoaquinValley_Tm3_16S_LsoF.fasta",
            "examples/CA_SanJoaquinValley_Tm3_50S_CL514F.fasta",
            "examples/CA_SanJoaquinValley_Tm4_16S_LsoF.fasta",
            "examples/CA_SanJoaquinValley_Tm4_50S_CL514F.fasta",
        ]
    )
    df_summary, df_blast, df_haplotyping = results_to_data_frames(results, DEFAULT_PARSE_RE)
    print(df_summary)
    print(df_blast)
    print(df_haplotyping)
    from .export import write_excel

    write_excel(df_summary, df_blast, df_haplotyping, "foo.xlsx")
