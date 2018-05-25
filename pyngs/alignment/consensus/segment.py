"""A module to represent segments in the alignment."""
import sys
from functools import total_ordering
from itertools import groupby
from operator import attrgetter, itemgetter
from collections import namedtuple

from .utils import decode_quality
import pyngs.alignment


def segmenter(alignment=None, quality_offset=32, content=None, default_quality=30):
    """
    Segment an alignment per base.

    Args:
        alignment - an alignment to segments
    Yields:
        a segment
    """
    cigar = pyngs.alignment.Cigar(alignment.cigar)
    rpos = alignment.position
    qpos = 0
    for nbase, cop in cigar.operations:
        seq = None
        qual = default_quality

        # return the bases in cigar operation separately
        for _ in range(nbase):
            # for operations on the query, calculate their quality
            if pyngs.alignment.Cigar.on_query(cop):
                seq = alignment.sequence[qpos]
                qual = decode_quality(alignment.quality[qpos], quality_offset)

            yield Segment(rpos, qpos, cop, seq, qual, content=content)

            if pyngs.alignment.Cigar.on_query(cop):
                qpos += 1
            if pyngs.alignment.Cigar.on_reference(cop):
                rpos += 1


@total_ordering
class Segment(object):
    """A class to represent a segment."""

    def __init__(self, refpos, qpos, operation, sequence, quality, content=None):
        """Initialize a segment."""
        self.refpos = refpos
        self.operation = operation
        self.sequence = sequence
        self.quality = quality
        self.qpos = qpos
        self.content = content
        self.anchored = []

    def __eq__(self, other):
        """2 segments are equal."""
        return (
            self.refpos,
            self.operation,
            self.qpos,
            self.sequence) == (
                other.refpos,
                other.operation,
                other.qpos,
                other.sequence)

    def __lt__(self, other):
        """One segment is less than the other."""
        return (
            self.refpos,
            self.operation,
            self.qpos,
            self.sequence) < (
                other.refpos,
                other.operation,
                other.qpos,
                other.sequence)

    def __repr__(self):
        return "{ref}:{op}:{qpos}:{seq}:{qual}:{content}:{anchored}".format(
            ref=self.refpos,
            op=self.operation,
            qpos=self.qpos,
            seq=self.sequence,
            qual=self.quality,
            content=repr(self.content),
            anchored=[repr(x) for x in self.anchored]
        )
