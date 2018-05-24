"""A module to represent segments in the alignment."""
import sys
from functools import total_ordering
from itertools import groupby
from operator import attrgetter, itemgetter
from collections import namedtuple
import pyngs.alignment


def segmenter(alignment=None, quality_offset=32, content=None):
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
        qual = None

        # return the bases in cigar operation separately
        for _ in range(nbase):
            # for operations on the query, calculate their quality
            if pyngs.alignment.Cigar.on_query(cop):
                seq = alignment.sequence[qpos]
                qual = ord(alignment.quality[qpos]) - quality_offset

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
        return "{ref}:{op}:{qpos}:{seq}:{qual}:{content}".format(
            ref=self.refpos,
            op=self.operation,
            qpos=self.qpos,
            seq=self.sequence,
            qual=self.quality,
            content=repr(self.content)
        )


class Anchor(object):
    """."""
    def __init__(self, segments):
        self.refpos = 

def combine_segments(segments):
    """Combine similar segments."""
    def compare(seg):
        return(
            seg.refpos,
            seg.operation,
            seg.qpos,
            seg.content,
            seg.sequence,
            seg.quality
        )

    # return object
    alignedbase = namedtuple(
        "alignedbase", ["refpos", "qpos", "sequence", "quality"])

    # sort the segments based on comparator above
    segments.sort(key=compare)
    backbone = []
    for key, grp in groupby(segments,
                            key=attrgetter(
                                "refpos",
                                "operation")):
        # skip cigar operations that we
        # currently do not support
        if key[1] in ("H", "P"):
            continue
        grp = list(grp)
        grp.sort(key=attrgetter("content", "qpos"))
        backbone.append((key, grp))

    #
    querypos = 0
    for key, grp in backbone:

        # process the M case
        if key[1] in ("M"):
            grp = list(grp)
            grp.sort(key=attrgetter("sequence"))
            seqs = []
            for seq, seg in groupby(grp, key=attrgetter("sequence")):
                seqs.append(seq, sum([x.quality for x in seg]))
            seqs.sort(key, itemgetter(1))
