"""Code for exporting the data frames generated from the ``workflow`` module."""

import pandas as pd

from .web.settings import (
    MIN_IDENTITY_GREEN,
    MIN_IDENTITY_YELLOW,
    FONT_COLOR_GREEN,
    FONT_COLOR_YELLOW,
    FONT_COLOR_RED,
    BG_COLOR_GREEN,
    BG_COLOR_YELLOW,
    BG_COLOR_RED,
)

#: Sheet name for summary sheet.
SHEET_SUMMARY = "Summary"
#: Sheet name for BLAST results.
SHEET_BLAST = "BLAST"
#: Sheet name for haplotyping results.
SHEET_HAPLOTYPING = "Haplotyping"


def write_excel(df_summary, df_blast, df_haplotyping, path):
    """Write the sheets to the given ``path``."""
    writer = pd.ExcelWriter(path, engine="xlsxwriter")

    # Write out data frames.
    df_summary.to_excel(writer, sheet_name=SHEET_SUMMARY, index=False)
    df_blast.to_excel(writer, sheet_name=SHEET_BLAST, index=False)
    df_haplotyping.to_excel(writer, sheet_name=SHEET_HAPLOTYPING, index=False)

    # Get workbook and sheets.
    workbook = writer.book
    sheet_summary = writer.sheets[SHEET_SUMMARY]
    sheet_blast = writer.sheets[SHEET_BLAST]

    # Setup formats.
    format_green = workbook.add_format({"bg_color": BG_COLOR_GREEN, "font_color": FONT_COLOR_GREEN})
    format_yellow = workbook.add_format(
        {"bg_color": BG_COLOR_YELLOW, "font_color": FONT_COLOR_YELLOW}
    )
    format_red = workbook.add_format({"bg_color": BG_COLOR_RED, "font_color": FONT_COLOR_RED})
    for f in format_green, format_yellow, format_red:
        f.set_num_format("0.0")

    # Setup conditional formats.
    cond_green = {
        "type": "cell",
        "criteria": "between",
        "maximum": 100,
        "minimum": MIN_IDENTITY_GREEN,
        "format": format_green,
    }
    cond_yellow = {
        "type": "cell",
        "criteria": "between",
        "maximum": MIN_IDENTITY_GREEN,
        "minimum": MIN_IDENTITY_YELLOW,
        "format": format_yellow,
    }
    cond_red = {
        "type": "cell",
        "criteria": "between",
        "maximum": MIN_IDENTITY_YELLOW,
        "minimum": 0,
        "format": format_red,
    }

    # Apply conditional formatting.
    c_summary = chr(ord("A") + df_summary.columns.get_loc("identity"))
    c_blast = chr(ord("A") + df_blast.columns.get_loc("identity"))
    for cond in cond_green, cond_yellow, cond_red:
        sheet_summary.conditional_format(
            "%s2:%s%d" % (c_summary, c_summary, df_summary.shape[0] + 1), cond
        )
        sheet_blast.conditional_format("%s2:%s%d" % (c_blast, c_blast, df_blast.shape[0] + 1), cond)

    writer.save()
