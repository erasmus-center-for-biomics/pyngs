import sys
import argparse
import gzip
import pyngs.bed.regions as regions


def parse_gtf_comment(comment):
    """Parse the GTF comment."""
    fields = {}
    for fld in comment.split(";"):
        tokens = fld.split(" ", 1)
        if len(tokens) == 2:
            fields[tokens[0]] = tokens[1].replace('"', "")
    return fields


def parse_gff_comment(comment):
    """Parse the comment."""
    fields = {}
    for fld in comment.split(";"):
        tokens = fld.split("=", 1)
        if len(tokens) == 2:
            fields[tokens[0]] = tokens[1]
    return fields


class ExtractData:

    def __init__(self, columns, encodings, fields):
        self.columns = columns
        self.encodings = encodings
        self.fields = fields

        # add the parsers
        self.parsers = []
        for encoding in self.encodings:
            if encoding == "none":
                self.parsers.append(None)
            elif encoding == "gff":
                self.parsers.append(parse_gff_comment)
            elif encoding == "gtf":
                self.parsers.append(parse_gtf_comment)

    def __call__(self, instream, outstream, sep="\t"):
        """Process the data in the input stream."""
        last_column = max(self.columns)

        # for each line in the input stream
        for line in instream:

            # remove the last endline
            line = line.rstrip("\n")
            tokens = line.split(sep)

            # skip lines with insufficient columns
            if len(tokens) < last_column:
                continue

            # parse the fields
            fields = []
            for idx, column in enumerate(self.columns):
                if self.parsers[idx] is None:
                    fields.append(tokens[column])
                else:
                    data = self.parsers[idx](tokens[column])
                    for key in self.fields[idx]:
                        if key in data.keys():
                            fields.append(data[key])
                        else:
                            fields.append("__NA__")

            # add the data to the output
            outstream.write("{0}\n".format("\t".join(fields)))


def extract_from_columns(args):
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

    # extract data from columns
    fields = [fld.split(",") for fld in args.fields]
    extractdata = ExtractData(
        columns=args.columns,
        fields=fields,
        encodings=args.encodings)
    extractdata(instream, outstream)

    # close the streams
    if instream is not sys.stdin:
        instream.close()
    if outstream is not sys.stdout:
        outstream.close()


if __name__ == "__main__":

    # parse the commandline arguments
    sparser = argparse.ArgumentParser(
        prog="extract_from_columns",
        description="""Parses  bedtools intersect -wao -sorted
            -a [file] -b [file] output.""")
    sparser.add_argument(
        "-i", "--input", dest="input",
        type=str, nargs="?", default="stdin",
        help="The file with overlaps.")
    sparser.add_argument(
        "-o", "--output", dest="output",
        type=str, nargs="?", default="stdout",
        help="The output tab-delimited text file.")
    sparser.add_argument(
        "--columns", dest="columns",
        nargs="+", type=int,
        help="The columns to parse from the input stream.")
    sparser.add_argument(
        "--fields", dest="fields",
        nargs="+", type=str,
        help="""The field to extract per column. If multiple fields
        should be extracted per column, commas can be used to
        separate the values.

        If the column has no encoding any
        string can be specified as long as it is not empty.

        If a field is not present in a column, __NA__ will be printed.""")
    sparser.add_argument(
        "--encodings", dest="encodings", choices=["none", "gff", "gtf"],
        nargs="+", type=str, help="""The encoding in the columns.""")
    sparser.set_defaults(func=extract_from_columns)

    # parse the arguments
    args = sparser.parse_args()
    args.func(args)
