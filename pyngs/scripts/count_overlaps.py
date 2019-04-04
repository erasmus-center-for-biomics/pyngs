import sys
import argparse
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

if __name__ == "__main__":

    # parse the commandline arguments
    sparser = argparse.ArgumentParser(
        prog="count_overlaps",
        description=""".""")
    sparser.add_argument(
        "-i", "--input", dest="input",
        type=str, nargs="?", default="stdin",
        help="The file with overlaps.")
    sparser.add_argument(
        "-o", "--output", dest="out",
        type=str, nargs="?", default="stdout",
        help="The output tab-delimited text file.")
    sparser.add_argument(
        "--type-a", dest="type_a",
        type=str, choices=["GTF", "BED3", "BED4", "BED5"], default="BED5",
        help="The input format of file a.")
    sparser.add_argument(
        "--type-b", dest="type_b",
        type=str, choices=["GTF", "BED3", "BED4", "BED5"], default="GTF",
        help="The input format of file b.")
    sparser.add_argument(
        "-s", "--strand", dest="strand",
        choices=["equals", "opposite", "eiter"], default="either",
        type=str, help="How to use the strand information of the overlaps.")
    sparser.add_argument(
        "-a", "--aggregate-per", dest="aggregate", default=["gene_id"],
        type=str, help="Aggregate the counts per object in the id columns.")
    sparser.add_argument(
        "-w", "--what", dest="what",
        default=["READNAME"],
        type=str, help="What to count.")
    sparser.set_defaults(func=pyngs.scripts.count_overlaps)

    # parse the arguments
    args = sparser.parse_args()
    args.func(args)
