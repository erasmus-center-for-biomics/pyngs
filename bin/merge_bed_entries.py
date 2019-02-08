"""A script to merge BED entries that represent the same DNA fragment.

This script should be used together with
cigar_to_bed.py.
"""

import argparse
import gzip
import sys
from operator import attrgetter
from itertools import groupby
from collections import namedtuple
from pyngs import bed


# a BED entry
BED = namedtuple("BED", [
    "chromosome", "start", "end",
    "comment", "score", "strand"])


def decode_comment(comment, psep=";", ssep="="):
    """Decode the comment."""
    retval = {}
    for token in comment.split(psep):
        parts = token.split(ssep)
        if len(parts) < 2:
            retval[parts[0]] = ""
        else:
            retval[parts[0]] = "=".join(parts[1:])
    return retval


def partition_bed(instream, target_tag="READNAME"):
    """Partition the BED file based on a comment field."""
    batch = []
    curtag = None

    # for each bed entry
    for entry in bed.reader(instream, ftype="bed5", sep="\t"):

        # decode the comment and
        tags = decode_comment(entry[3])

        # skip entries without the tag
        if target_tag not in tags:
            continue

        # check if tag matches
        if tags[target_tag] != curtag:

            # yield the batch if present
            if len(batch):
                yield curtag, batch

            # set for the next tag
            curtag = tags[target_tag]
            batch = []
        # add the current entry to the batch
        batch.append(BED(*entry))

    # yield the last data in the batch
    if len(batch):
        yield curtag, batch


def merge_entries(batch, comment):
    """Merge the entries for a batch."""
    def split_batch_on_chromosome(batch):
        """Split read batches on chromosome."""
        batch.sort(key=attrgetter("chromosome"))
        for chrom, values in groupby(attrgetter("chromosome")):
            yield chrom, values

    # split entries per chromosome
    for chrom, entries in split_batch_on_chromosome(batch):
        mrg = BED(
            chrom,
            min([en.start for en in entries]),
            max([en.end for en in entries]),
            comment, 0, ".")
        yield mrg


def merge_paired_bed(instream, outstream, tag="READNAME"):
    """Merge paired BED entries."""
    outstring = "{chromosome}\t{start}\t{end}\t{comment}\t{score}\t{strand}\n"
    for value, group in partition_bed(instream, target_tag=tag):
        for merged in merge_entries(group, value):
            outstream.write(
                outstring.format(
                    merged.chromosome,
                    merged.start,
                    merged.end,
                    merged.comment,
                    merged.score,
                    merged.starnd))


def main():
    """Parse the commandline parameters and run to_bed."""
    parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        description="""
        A script to merge BED entries that
        represent the same DNA fragment.

        This script should be used together with
        cigar_to_bed.py.
        """
    )
    parser.add_argument(
        "-i", "--input", dest="input",
        type=str, nargs="?", default="stdin",
        help="""The BED file with the entries to process.""")
    parser.add_argument(
        "-o", "--ouput", dest="ouput",
        type=str, nargs="?", default="stdout",
        help="The path to the BED with the merged entries.")
    parser.add_argument(
        "-t", "--tag", dest="tag",
        type=str, nargs="?", default="READNAME",
        help="The tag on which to merge the entries.")
    args = parser.parse_args()

    # open the in and output files
    instream = sys.stdin
    if args.sam != "stdin":
        if args.sam.endswith(".gz"):
            instream = gzip.open(args.sam, "rt")
        else:
            instream = open(args.sam, "rt")

    outstream = sys.stdout
    if args.bed != "stdout":
        if args.bed.endswith(".gz"):
            outstream = gzip.open(args.bed, "wt")
        else:
            outstream = open(args.bed, "wt")

    # write the BED entries
    merge_paired_bed(instream, outstream, args.tag)

    # close the streams
    if instream != sys.stdin:
        instream.close()
    if outstream != sys.stdout:
        outstream.close()


if __name__ == "__main__":

    # run the main program entry point
    main()
