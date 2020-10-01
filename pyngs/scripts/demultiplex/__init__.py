import argparse
import sys
import multiprocessing
from collections import namedtuple


def read_baits(stream):
    """Read bait sequences from a stream."""
    retval = []
    Bait = namedtuple("Bait", ["sample", "sequence"])

    for line in stream:
        line = line.rstrip()
        if not line or line.startswith("#"):
            continue
        tokens = line.split("\t")
        retval.append(Bait(*tokens))
    return retval


def demultiplex(args):
    """."""
    pass


if __name__ == "__main__":

    # parse the commandline arguments
    sparser = argparse.ArgumentParser(
        prog="demultiplex",
        description="""Demultiplex reads using the Levenshstein algorithm.""")
    sparser.add_argument(
        "-b", "--baits", dest="baits",
        type=str,
        help="A set of bait sequences.")
    sparser.add_argument(
        "-f", "--fastq", dest="fastq",
        type=str, nargs="*", default="stdin",
        help="An input FastQ file with single reads.")
    sparser.add_argument(
        "-i", "--index", dest="index_read",
        type=str, nargs="?", default=None,
        help="An (optional) input FastQ file with the index reads.")
    sparser.add_argument(
        "-o", "--output", dest="output",
        type=str, nargs="?", default="stdout",
        help="An output table with the read assignments.")
    sparser.add_argument(
        "-m", "--maximum-score", dest="maximum_score",
        type=int, default=-1,
        help="The maximum score.")
    sparser.add_argument(
        "-f", "--follow-up-score", dest="followup",
        type=int, default=-1,
        help="The minimum score difference between the first and second hits.")
    sparser.add_argument(
        "-w", "--workers", dest="workers",
        type=int, default=4,
        help="The number of concurrent worker processes.")
