import sys
import argparse
import gzip
import pyngs.bed.regions as regions


def strandfactory(cmp="any"):
    """Get the strand comparison."""
    def any(a, b):
        return True

    def same(a, b):
        return a == b

    def opposite(a, b):
        return a != b

    if cmp == "any":
        return any
    elif cmp == "same":
        return same
    elif cmp == "opposite":
        return opposite
    else:
        return None


def parse_gtf_comment(comment):
    """Parse the GTF comment."""
    fields = {}
    for fld in comment.split(";"):
        tokens = fld.split(" ", 1)
        if len(tokens) == 2:
            fields[tokens[0]] = tokens[1].replace('"', "")
    return fields


def parse_comment(comment):
    """Parse the comment."""
    fields = {}
    for fld in comment.split(";"):
        tokens = fld.split("=", 1)
        if len(tokens) == 2:
            fields[tokens[0]] = tokens[1]
    return fields


def next_overlap(instream, sep="\t"):
    """Get the next overlap from a stream."""
    pa_bed = regions.get_parser("BED6")
    pa_gff = regions.get_parser("GFF")

    # for each line
    for line in instream:
        line = line.rstrip("\n")
        tokens = line.split(sep)

        # if there are insufficient columns
        if len(tokens) < 15:
            continue

        # parse the bed and gff portions of the line
        bed, tokens = pa_bed(tokens)
        gff, tokens = pa_gff(tokens)

        # yield the bed and gff regions
        yield bed, gff


class ExtractData:

    def __init__(self, strand, which, who, what, is_gtf):
        self.strand = lambda a, b: True
        self.which = which
        self.who = who
        self.what = what
        self.is_gtf = is_gtf

    def __call__(self, instream, outstream, sep="\t"):
        """."""
        # header = ["/".join(self.which)]
        # header.extend(self.who)
        # header.append(self.what)

        # outstream.write("#" + "\t".join(header) + "\n")

        for bed, gff in next_overlap(instream, sep=sep):

            # skip entries with the wrong orientation
            if not self.strand(
                    bed.content["strand"],
                    gff.content["strand"]):
                continue

            # skip other types
            which = gff.content["type"]
            if which not in self.which:
                continue

            # parse the GFF comment and get the whos
            if self.is_gtf:
                gfields = parse_gtf_comment(gff.content["comment"])
            else:
                gfields = parse_comment(gff.content["comment"])
            whos = [gfields[w] for w in self.who if w in gfields.keys()]
            if len(whos) != len(self.who):
                continue

            # parse the BED comment and get the whats
            bfields = parse_comment(bed.content["comment"])
            if self.what not in bfields.keys():
                continue
            what = bfields[self.what]

            entry = [which]
            entry.extend(whos)
            entry.append(what)
            outstream.write("\t".join(entry) + "\n")


def extract_overlaps(args):
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
    extractdata = ExtractData(
        strandfactory(args.strand),
        args.which, args.who, args.what, args.gtf)
    extractdata(instream, outstream)

    # close the streams
    if instream is not sys.stdin:
        instream.close()
    if outstream is not sys.stdout:
        outstream.close()


if __name__ == "__main__":

    # parse the commandline arguments
    sparser = argparse.ArgumentParser(
        prog="extract_overlaps",
        description="""Parses  bedtools intersect -wao -sorted
            -a [BED5 file] -b [GTF file] output.""")
    sparser.add_argument(
        "-i", "--input", dest="input",
        type=str, nargs="?", default="stdin",
        help="The file with overlaps.")
    sparser.add_argument(
        "-o", "--output", dest="output",
        type=str, nargs="?", default="stdout",
        help="The output tab-delimited text file.")
    sparser.add_argument(
        "-s", "--strand", dest="strand",
        choices=["same", "opposite", "any"], default="any",
        type=str, help="How to use the strand information of the overlaps.")
    sparser.add_argument(
        "--which", dest="which", default=["exon"],
        type=str, help="Which types in the GTF file to count.")
    sparser.add_argument(
        "--who", dest="who", default=["Parent"],
        type=str, help="""The field with the ids in the GFF
        file to aggregate over.""")
    sparser.add_argument(
        "--what", dest="what", default="READNAME",
        type=str, help="What to count.")
    sparser.add_argument(
        "--gtf", dest="gtf", default=False,
        type=bool, help="the intersect was performed with a GTF file")
    sparser.set_defaults(func=extract_overlaps)

    # parse the arguments
    args = sparser.parse_args()
    args.func(args)
