import sys
import argparse
import gzip


def count(instream):
    """Count entries."""
    pkey = None
    pvalue = None
    score = 0
    for line in instream:

        # parse the list
        tokens = line.rstrip().split("\t")
        key = tokens[0:(len(tokens) - 1)]
        value = tokens[-1]

        if pkey and key != pkey:
            yield pkey, score
            score = 0

        if value is not pvalue:
            score += 1
        pkey, pvalue = key, value

    if pkey:
        yield pkey, score


def count_last_column(args):
    """Count the number of entries in the last column."""
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

    for key, score in count(instream):
        outstream.write("{0}\t{1}\n".format(
            "\t".join(key),
            score))

    # close the streams
    if instream is not sys.stdin:
        instream.close()
    if outstream is not sys.stdout:
        outstream.close()


if __name__ == "__main__":

    # parse the commandline arguments
    sparser = argparse.ArgumentParser(
        prog="count_last_column",
        description="""Count the number of
        different entries in the last column from a sorted
        tab-delimited text file.""")
    sparser.add_argument(
        "-i", "--input", dest="input",
        type=str, nargs="?", default="stdin",
        help="A tab-delimited text file.")
    sparser.add_argument(
        "-o", "--output", dest="output",
        type=str, nargs="?", default="stdout",
        help="The output tab-delimited text file.")

    sparser.set_defaults(func=count_last_column)

    # parse the arguments
    args = sparser.parse_args()
    args.func(args)
