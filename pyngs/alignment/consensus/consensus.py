"""A module to calculate the consensus alignments from multiple alignments."""
import sys
from operator import attrgetter
from itertools import groupby
from .segment import segmenter


def consensus(alignments, quality_offset=32):
    """Get the consensus from a set of alignments.

    The consensus is determined by first segmenting the
    alignments based on their CIGAR operations. Each individual
    operation becomes a segment. These segments are then first
    solved per reference position. Afterwards, soft-clipped,
    insertions and deletions are solved.
    """
    # get the segments for the alignments
    segments = []
    for idx in range(len(alignments)):
        segments.extend(
            list(
                segmenter(alignments[idx],
                          quality_offset=quality_offset,
                          content=idx)))
    segments.sort(key=attrgetter("content", "refpos", "qpos"))
    for key, seg in groupby(segments, key=attrgetter("content", "refpos")):
        sys.stderr.write("{key}\t{n}\n".format(key=key, n=len(list(seg))))

