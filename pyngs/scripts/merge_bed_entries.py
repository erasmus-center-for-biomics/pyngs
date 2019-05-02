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


def unique(seq):
    """Get and count the unique entries in seq."""
    prev = None
    count = 0

    for item in seq:
        if count and item != prev:
            yield prev, count
            count = 0
            prev = item
        count += 1

    if count:
        yield prev, count


def merge_entries(batch, comment):
    """Merge the entries for a batch."""
    # split entries per chromosome
    batch.sort(key=attrgetter("chromosome"))
    for chrom, values in groupby(batch, key=attrgetter("chromosome")):

        #
        starts, ends, strands = [], [], []
        for entry in values:
            starts.append(entry.start)
            ends.append(entry.end)
            strands.append(entry.strand)

        # determine the strand
        strands.sort()
        ustrand = [e for e in unique(strands)]
        strand = ustrand[0][0] if len(ustrand) == 1 else "."

        # print the merged entry
        mrg = BED(chrom, min(starts), max(ends), comment, 0, strand)
        yield mrg


def merge_paired_bed(instream, outstream, tag="READNAME"):
    """Merge paired BED entries."""
    outstring = "{chromosome}\t{start}\t{end}\t{comment}\t{score}\t{strand}\n"
    for value, group in partition_bed(instream, target_tag=tag):
        for merged in merge_entries(group, value):
            outstream.write(
                outstring.format(
                    chromosome=merged.chromosome,
                    start=merged.start,
                    end=merged.end,
                    comment=merged.comment,
                    score=merged.score,
                    strand=merged.strand))


def merge_bed_entries(args):
    """
    Merge BED entries that represent the same DNA fragment.

    This function is meant to be used together with cigar-to-bed.
    """
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

    # write the BED entries
    merge_paired_bed(instream, outstream, args.tag)

    # close the streams
    if instream != sys.stdin:
        instream.close()
    if outstream != sys.stdout:
        outstream.close()
