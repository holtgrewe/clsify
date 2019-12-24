"""Haplotype-Lso web front-end settings."""

# Not Configurable ===========================================================

#: String to display on home screen.
HOME_BRAND = "Haplotype-LSO"

# TODO: should this come from the UI?
#: The sample regex
SAMPLE_REGEX = r"^(?P<sample>[^\.]+)\.(?P<region>[^\.]+)(?:\.(?P<primer>.*))?"

#: Whether to write file name as sample name.
FILE_NAME_TO_SAMPLE_NAME = True

#: Known regions.
KNOWN_REGIONS = ("16S", "16S-23S", "50S")
#: Known primers.
KNOWN_PRIMERS = ("LsoF", "Ol2cR", "LsoTX1623F", "LsoTX1623R", "Cl514F", "Cl514R")

#: Font color for red conditional formatting.
FONT_COLOR_RED = "#CC0000"
#: Background color for red conditional formatting.
BG_COLOR_RED = "#FFCCCC"
#: Font color for yellow conditional formatting.
FONT_COLOR_YELLOW = "#996600"
#: Background color for yellow conditional formatting.
BG_COLOR_YELLOW = "#FFFFCC"
#: Font color for green conditional formatting.
FONT_COLOR_GREEN = "#006600"
#: Background color for green conditional formatting.
BG_COLOR_GREEN = "#CCFFCC"

#: Minimal identity (in percent) for green color.
MIN_IDENTITY_GREEN = 98
#: Minimal identity (in percent) for yellow color.
MIN_IDENTITY_YELLOW = 96


# Configurable ===============================================================

#: The host to listen on.
HOST = "0.0.0.0"
#: The port to listen on.
PORT = 8050

#: The public URL prefix to use.
PUBLIC_URL_PREFIX = ""
