"""
Methods to generate consensus alignments from
multiple other alignments.
"""
import collections
from operator import itemgetter

from . import quality_to_score
from .cigar import operations, CigarOperation
from .cigar import CIGAR_OPERATIONS_ON_QUERY
from .cigar import CIGAR_OPERATIONS_ON_REFERENCE


# define a consensus object
ConsensusObject = collections.namedtuple("ConsensusObject", [
    "code", "length", "reference",
    "sequence", "quality", "alignments"])


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
        """Increment reference and query."""
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
        """Increment the query."""
        delta = cigop.query[1] - cigop.query[0]
        for offset in range(delta):
            # the end position of the query has
            # been repurposed to the relative
            # offset in the cigar.
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
        """Increment the reference."""
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


def consensus_sorter(obj: tuple):
    """Sort cigar operations for the consensus creation."""
    return (
        obj[0].reference[1],
        obj[0].code,
        obj[0].sequence,
        obj[1],
        obj[0].query[1])


def merge_non_reference(entries: list):
    """Merge non reference entries."""
    def group(entries: list):
        """Group non reference entries."""
        cidx = None
        operations = []
        for (cigop, aidx) in entries:
            if aidx != cidx:
                if len(operations):
                    yield operations
                cidx = aidx
                operations = []
            operations.append(cigop, aidx)

        if len(operations):
            yield operations


def create_consensus_operation(entries: list, qual_offset=32):
    """Create a consensus operation from cigar-operations."""
    quality = 0
    for cigop, _ in entries:
        if not len(cigop.quality):
            continue
        quality += quality_to_score(cigop.quality, qual_offset)

    base = entries[0][0]
    return ConsensusObject(
        code=base.code, length=base.length,
        reference=base.reference,
        sequence=base.sequence, quality=quality, alignments=len(entries))


def consensus_operations(operations):
    """Group consensus operations per reference position."""
    def emit(ref, qry):
        """Determine whether the return buffer should be emitted."""
        if ref.reference[1] != qry.reference[1]:
            return True
        if ref.code != qry.code:
            return True
        if ref.sequence != qry.sequence:
            return True
        return False

    def

    entries = []
    for (cigop, aidx) in operations:
        if entries and emit(entries[0][0], cigop):
            # special case for the query only operations
            if not entries[0][0].code in CIGAR_OPERATIONS_ON_REFERENCE:
                # TODO resort the retvals to the alignments, merge
                # operations per code and move it of to consensus_operation
                merge_cigar_operations(entries)
                pass
            else:
                yield create_consensus_operation(entries)
            entries = []
        entries.append((cigop, aidx))


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
        parts.sort(key=consensus_sorter)

        for idx, conop in enumerate(consensus_operations(parts)):
            print(idx, conop)
