import sys
import gzip
import argparse
from typing import Generator, Tuple
from pyngs.fastq import fastq


def subreads(provider: Generator[Tuple[str, str, str], None, None], start: int=0, end: int=-1):
    """Shorten the reads from the provider."""
    for name, sequence, quality in provider:
        if len(sequence) < start:
            yield name, "", ""
        elif end < 0 or len(sequence) < end:
            yield name, sequence[start:], quality[start:]
        else:
            yield name, sequence[start:end], quality[start:end]


def extract_subreads(args):
    """Shorten the reads in a fastq file to length."""
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

    # check the start of the read
    if args.start < 0:
        raise ValueError("The sub-read start cannot be smaller than 0")

    # shorten the reads
    provider = fastq(instream)
    for name, sequence, quality in subreads(provider, args.start, args.end):
        outstream.write("@{0}\n{1}\n+\n{2}\n".format(
            name, sequence, quality))
    
    # close the in and output files
    if instream is not sys.stdin:
        instream.close()
    if outstream is not sys.stdout:
        outstream.close()


if __name__ == "__main__":

    sparser = argparse.ArgumentParser(
        prog="extract_subread_from_fastq",
        description="""Extract a subread from a FastQ file.""")
    sparser.add_argument(
        "-i", "--input", dest="input",
        type=str, default="stdin",
        help="""The input FastQ file.""")
    sparser.add_argument(
        "-o", "--output", dest="output",
        type=str, default="stdout",
        help="The output FastQ file.")
    sparser.add_argument(
        "-s", "--start", dest="start",
        type=int, default=0,
        help="""The start of the sub-read.""")
    sparser.add_argument(
        "-e", "--end", dest="end",
        type=int, default=-1,
        help="""The end of the sub-read.""")

    sparser.set_defaults(func=extract_subreads)
    args = sparser.parse_args()
    args.func(args)
