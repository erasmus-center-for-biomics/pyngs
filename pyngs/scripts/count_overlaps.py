import sys
import gzip


def select_parser(name):
    """Select the parser based on a name."""
    if name == "BED3":
        pass
    elif name == "BED4":
        pass
    elif name == "BED5":
        pass
    elif name == "GTF":
        pass


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
