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


SplitOp = collections.namedtuple(
    "SplitOp", 
    ["code", "length", "reference", "sequence", "quality"])


def split_operations(alignment):
    """."""
    
    for cigarop in operations(alignment):
        
        # return non reference operations as is
        if cigarop.code not in CIGAR_OPERATIONS_ON_REFERENCE:
            yield SplitOp(
                code=cigarop.code,
                length=cigarop.length,
                reference=cigarop.reference,
                sequence=cigarop.sequence,
                quality=cigarop.quality)
        else:
            # split reference operations per base
            for offset in range(0, cigarop.length):
                seq = ""
                qual = ""
                if cigarop.code in CIGAR_OPERATIONS_ON_QUERY:
                    seq = cigarop.sequence[offset],
                    qual = cigarop.quality[offset] 
                yield SplitOp(
                    code=cigarop.code,
                    length=1,
                    reference=[
                        cigarop.reference[0],
                        cigarop.reference[1] + offset,
                        cigarop.reference[1] + offset + 1],
                    sequence=seq,
                    quality=qual)


def splitobs_sort(obj):
    """Sort cigar operations for the consensus creation."""
    return (
        obj.reference[1],
        obj.reference[2],
        obj.code,
        obj.length,
        obj.sequence)


def differs(splitop_a, splitop_b):
    """Check whether 2 consensus alignments differ."""
    if splitop_a.reference != splitop_b.reference:
        return True
    if splitop_a.code != splitop_b.code:
        return True
    if splitop_a.length != splitop_b.length:
        return True
    if splitop_a.sequence != splitop_b.sequence:
        return True
    return False


def same_operation(splitops):
    """Batch the same operations in the consensus alignment."""
    batch = []
    for splitop in splitops:
        if batch:
            if differs(batch[0], splitop):
                yield batch
                batch = []
        batch.append(splitop)
    if batch:
        yield batch


def by_position(batches):
    """Group the consensus alignments by position."""
    sets = []
    for batch in batches:
        if sets:
            if sets[0][0].reference != batch[0].reference:
                yield sets
                sets = [] 
        sets.append(batch)

    if sets:
        yield sets


def performance_sort(obj):
    return (obj[0][0], obj[0][1])


class Consensus:
    """A class to generate consensus alignments."""

    def __init__(self, quality_offset: int=32, default_qual: int=30):
        self.quality_offset = quality_offset
        self.default_qual = default_qual

    def __call__(self, alignments: list):
        """Generate a new consensus alignment."""
        parts = []
        for alignment in alignments:
            for cigarop in split_operations(alignment):
                parts.append(cigarop)
        parts.sort(key=splitobs_sort)
        
        # determine the preliminary consensus
        preliminary = []
        for batches in by_position(same_operation(parts)):
            # get the performance of the possible entries per 
            # position in the consensus 
            to_choose = []
            for batch in batches:
                meas = self.performance(batch)
                to_choose.append((meas, batch[0]))
            
            # append the top hit to the preliminary consensus
            to_choose.sort(key=performance_sort, reverse=True)
            preliminary.append(to_choose[0])

        # TODO remove insertions with fewer than half the reads 
        # of the surrounding bases, internal soft-clipped and
        # hard-clipped bases.
        cleaned = []
        for idx, (meas, splitop) in enumerate(preliminary):
            
            # decide to keep or skip insertions
            if splitop.code in "I":
                if idx > 0 and meas[0] < preliminary[idx-1][0][0] / 2:
                    continue
                if idx < len(preliminary) - 1 and meas[0] < preliminary[idx+1][0][0] / 2:
                    continue

            # remove internal clipped bases
            if splitop.code in "SH":
                if idx > 0 and idx < len(preliminary) - 1:
                    continue
            cleaned.append((meas, splitop))
        
        # add N CIGAR operations for bases covered by the alignment
        # but absent in the reference. 
        consensus = []
        for idx, (meas, consop) in enumerate(cleaned):

            # only check after the first
            if idx > 0:
                prevop = cleaned[idx-1][1]

                # insert an N stretch if the previous operation
                # does not end at the current operations
                if prevop.reference[2] != consop.reference[1]:
                    size = consop.reference[1] - prevop.reference[2] 
                    insert = SplitOp(
                        code="N", length=size,
                        reference=[
                            consop.reference[0],
                            prevop.reference[2],
                            consop.reference[1]],
                        sequence="", quality="")
                    consensus.append((0, 0.0), insert)
            consensus.append((meas, consop))

        # return the consensus operation
        return consensus
     


    
    def qual_to_score(self, qual):
        vals = [ord(q) - self.quality_offset for q in qual]
        return sum(vals) / len(qual)

    def performance(self, batch):
        """."""
        quals = 0
        if batch[0].code in CIGAR_OPERATIONS_ON_QUERY:
            quals = sum([self.qual_to_score(b.quality) for b in batch])
        return len(batch), quals