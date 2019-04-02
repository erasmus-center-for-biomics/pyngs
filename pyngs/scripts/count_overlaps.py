import sys
import gzip
import pyngs.bed.parsers
import pyngs.bed as bed


def select_parser(name):
    """Select the parser based on a name."""
    if name == "BED3":
        return pyngs.bed.parsers.parse_bed3
    elif name == "BED4":
        return pyngs.bed.parsers.parse_bed4
    elif name == "BED5":
        return pyngs.bed.parsers.parse_bed5
    elif name == "GTF":
        return pyngs.bed.parsers.parse_gff
    return None


def count_overlaps(args):
    """Count the overlaps between a BED and GTF file."""
    # open the input and output files
    instream = sys.stdin
    if args.input != "stdin":
        if args.input.endswith(".gz"):
            instream = gzip.open(args.input, "rt")
        else:
            instream = open(args.input, "rt")
    outstream = sys.stdout
    if args.output != "stdout":
        if args.output.endswith(".gz"):
            outstream = gzip.open(args.output, "wt")
        else:
            outstream = open(args.output, "wt")

    #
