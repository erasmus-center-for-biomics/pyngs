"""
Methods to generate consensus alignments from
multiple other alignments.
"""
import operator
from .cigar import operations, CigarOperation
from .cigar import CIGAR_OPERATIONS_ON_QUERY
from .cigar import CIGAR_OPERATIONS_ON_REFERENCE


def split_operations(alignment: list):
    """
    Get annotation CIGAR operations per increment.
    :param alignment: an alignment to process
    :yield: a CIGAR operation of length 1 for
            operations on the reference and the query.

    Note: in the case of the query only operations
          the query end position has been replaced
          with the relative coordinate to the next
          reference position. So:

          I  I  I  I  M
          0  1  2  3

          This should make the consensus calling easier
    """

    def dual(cigop):
        """A dual increment generator."""
        delta = cigop.reference[2] - cigop.reference[1]
        for offset in range(delta):
            yield CigarOperation(
                code=cigop.code,
                length=1,
                reference=(
                    cigop.reference[0],
                    cigop.reference[1] + offset,
                    cigop.reference[1] + offset + 1),
                query=(
                    cigop.query[0],
                    cigop.query[0], + offset + 1),
                sequence=cigop.sequence[offset],
                quality=cigop.quality[offset])

    def query(cigop):
        """A query increment generator."""
        delta = cigop.query[1] - cigop.query[0]
        for offset in range(delta):
            # the end position of the query has
            # been repurposed to the relative
            # offset in the cigar. This to make
            # a
            yield CigarOperation(
                code=cigop.code,
                length=1,
                reference=cigop.reference,
                query=(
                    cigop.query[0] + offset,
                    offset),
                sequence=cigop.sequence[offset],
                quality=cigop.quality[offset])

    def reference(cigop):
        """A reference increment generator."""
        delta = cigop.reference[2] - cigop.reference[1]
        for offset in range(delta):
            yield CigarOperation(
                code=cigop.code,
                length=1,
                reference=(
                    cigop.reference[0],
                    cigop.reference[1] + offset,
                    cigop.reference[1] + offset + 1),
                query=cigop.query,
                sequence=cigop.sequence,
                quality=cigop.quality)

    for cigarop in operations(alignment):
        if cigarop.code in CIGAR_OPERATIONS_ON_REFERENCE and \
                cigarop.code in CIGAR_OPERATIONS_ON_QUERY:
            for cop in dual(cigarop):
                yield cop
        elif cigarop.code in CIGAR_OPERATIONS_ON_REFERENCE:
            for cop in reference(cigarop):
                yield cop
        elif cigarop.code in CIGAR_OPERATIONS_ON_QUERY:
            for cop in query(cigarop):
                yield cop
        else:
            yield cigarop


def sorter(obj: tuple):
    """A function to get the attributes to sort on."""
    return (
        obj[0].reference[1],
        obj[0].code,
        obj[0].sequence,
        obj[1],
        obj[0].query[1])


class Consensus:
    """A class to generate consensus alignments."""

    def __init__(self, quality_offset: int=32, default_qual: int=30):
        self.quality_offset = quality_offset
        self.default_qual = default_qual

    def __call__(self, alignments: list):
        """Generate a new consensus alignment."""
        parts = []
        for idx, alignment in enumerate(alignments):
            for cigarop in split_operations(alignment):
                parts.append((cigarop, idx))
        parts.sort(key=sorter)

        for idx, (cigop, aidx) in enumerate(parts):
            pass
            # print(
            #     idx,
            #     cigop.reference[1],
            #     cigop.code,
            #     cigop.sequence,
            #     aidx,
            #     cigop.query[0],
            #     cigop.query[1])
